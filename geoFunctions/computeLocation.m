function [posECEF, gnssVelECEF, dopMatrix] = computeLocation(p, pdot, satPos, satVel, initPosECEF, initVelECEF)

%COMPUTELOCATION Compute receiver position from GNSS measurements
%
%   Inputs
%       p              - Pseudoranges. S-element column vector. S is the 
%                        number of satellites.
%       pdot           - Pseudorange rates. S-element column vector. S is 
%                        the number of satellites.
%       satPos         - Satellite positions in ECEF (m). S-by-3 matrix. S 
%                        is the number of satellites.
%       satVel         - Satellite velocities in ECEF (m/s). S-by-3 matrix.
%                        S is the number of satellites.
%       initPosECEF    - Initial position estimate in ECEF (m). 3-element 
%                        row vector.
%       intVelECEF     - Initial velocity estimate in ECEF (m/s). 3-element
%                        row vector.
%   Outputs
%       posECEF        - Position estimate in ECEF (m). 3-element row 
%                        vector.
%       gnssVelECEF    - Velocity estimate in ECEF (m/s). 3-element row 
%                        vector.
%       dopMatrix      - 4x4 dilution of precision in ECEF. Diagonals are
%                        (x, y, z, c*tau) uncertainties (m^2).

if (numel(p) < 4)
    posECEF = NaN(1, 3);
    gnssVelECEF = NaN(1, 3);
    dopMatrix = NaN(4,4);
    return;
end

posPrev = [initPosECEF(:); 0];
velPrev = [initVelECEF(:); 0];

posEst = posPrev;
velEst = velPrev;

Hpos = ones(size(satPos, 1), 4);
Hvel = ones(size(satVel, 1), 4);
HposPrev = Hpos; % Used for DOP calculation.

resPos = Inf;
minResPos = 1e-4;
resVel = Inf;
minResVel = 1e-4;
maxIterations = 200;
iter = 0;
allResPos = NaN(maxIterations, 1, 'like', p);
allResVel = NaN(maxIterations, 1, 'like', p);
% Check if residuals are increasing, if so, save previous estimate that
% corresponds to smaller residual.
checkConverge = @(x) issorted(x, 'descend', 'MissingPlacement', 'last');
while (resPos > minResPos) && (iter < maxIterations)
    % Obtain current true range estimate and line-of-sight vector.
   % [pEst, ~, losVector] = nav.internal.gnss.calculatePseudoranges(satPos, satVel, posPrev(1:3).', velPrev(1:3).');
 recPos=  posPrev(1:3).' ;
 recVel= velPrev(1:3).' ;
 numSats = size(satPos, 1);
numAx = size(satPos, 2);

repRecPos = repmat(recPos, numSats, 1);
repRecVel = repmat(recVel, numSats, 1);

% Earth's rotation rate (rad/s).
[~, ~, ~, OmegaEDot] = fusion.internal.frames.wgs84ModelParams;
% Speed of light (m/s).
c = fusion.internal.ConstantValue.SpeedOfLight;

% Ranges without accounting for Earth's rotation.
posDiff = satPos - repRecPos;
rawRanges = vecnorm(posDiff, 2, 2);
losVector = posDiff ./ repmat(rawRanges, 1, numAx);
deltaTimes = permute(rawRanges ./ c, [3 2 1]);

rotECEF2ECI = repmat(eye(3), 1, 1, numel(deltaTimes));
rotECEF2ECI(1,1,:) = cos(OmegaEDot .* deltaTimes);
rotECEF2ECI(1,2,:) = sin(OmegaEDot .* deltaTimes);
rotECEF2ECI(2,1,:) = -sin(OmegaEDot .* deltaTimes);
rotECEF2ECI(2,2,:) = cos(OmegaEDot .* deltaTimes);

permSatPos = repmat( permute(satPos, [3 2 1]), numAx, 1 );
% Rotate each satellite position by each corresponding rotation matrix. This
% is the equivalent of a vectorized matrix multiply:
%   for s = 1:numSatellites
%       rotSatPos(s,:) = ( rotECEF2ECI(:,:,s) * satPos(s,:).' ).';
%   end
rotSatPos = permute(sum(rotECEF2ECI .* permSatPos, 2), [3 1 2]);

% Calculate pseudorange.
posDiff = rotSatPos - repRecPos;
p = vecnorm(posDiff, 2, 2);

skewSymOmegaEDot = [        0, -OmegaEDot, 0; 
                    OmegaEDot,          0, 0;
                            0,          0, 0];
rotSatVel = zeros(numSats, 3);
% Rotate each satellite velocity by each corresponding rotation matrix.
for s = 1:numSats
    rotSatVelECEF = satVel(s,:) + (skewSymOmegaEDot * satPos(s,:).').';
    rotSatVel(s,:) = ( rotECEF2ECI(:,:,s) * rotSatVelECEF.' ).';
end

% Calculate pseudorange rate.
pdot = dot(losVector, rotSatVel - (repRecVel + (skewSymOmegaEDot * repRecPos.').'), 2);
pEst = pdot;
    
    % Add previous clock bias error (m). This is the time offset of the 
    % receiver clock times the speed of light.
    pPred = pEst + posPrev(end);
    
    Hpos(:,1:3) = -losVector;
    
    posEst = posPrev + Hpos.' * Hpos \ Hpos.' * (p - pPred);

    resPos = norm(posEst - posPrev);
    
    iter = iter + 1;
    allResPos(iter) = resPos;
    if ~checkConverge(allResPos)
        posEst = posPrev;
        Hpos = HposPrev;
        break;
    end
    posPrev = posEst;
    HposPrev = Hpos;
end

iter = 0;
while (resVel > minResVel) && (iter < maxIterations) && checkConverge(allResVel)
    % Obtain current true range rate estimate and line-of-sight vector.
    [~, pdotEst, losVector] ...
        = nav.internal.gnss.calculatePseudoranges(satPos, satVel, ...
        posEst(1:3).', velPrev(1:3).');
    % Add previous clock drift error (m/s). This is the time drift of the 
    % receiver clock times the speed of light.
    pdotPred = pdotEst + velPrev(end);
    
    Hvel(:,1:3) = -losVector;
    
    velEst = velPrev + Hvel.' * Hvel \ Hvel.' * (pdot - pdotPred);
    
    resVel = norm(velEst - velPrev);
    
    iter = iter + 1;
    allResVel(iter) = resVel;
    if ~checkConverge(allResVel)
        velEst = velPrev;
        break;
    end
    velPrev = velEst;
end
posECEF = posEst(1:3).';

gnssVelECEF = velEst(1:3).';

dopMatrix = inv(Hpos.' * Hpos);
end
