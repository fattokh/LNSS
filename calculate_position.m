close all;
clc;
clearvars -except result;

% Load ground truth trajectory.
load("routeNatickMA.mat","lat","lon","pos","vel","lla0");
recPos = pos;
recVel = vel;
input_file = readmatrix("output_file2.csv");
input_files = readtable("output_file2.csv");

in_edit = zeros(height(input_file)*4, 7);
in_edit_sat1 = zeros(height(input_file), 7);
in_edit_sat2 = zeros(height(input_file), 7);
in_edit_sat3 = zeros(height(input_file), 7);
in_edit_sat4 = zeros(height(input_file), 7);

in_edit_sat1(:,1) =input_file(:,1);
in_edit_sat2(:,1) = input_file(:,1);
in_edit_sat3(:,1) = input_file(:,1);
in_edit_sat4(:,1) =  input_file(:,1);

in_edit_sat1(:,2) = 1;
in_edit_sat2(:,2) = 2;
in_edit_sat3(:,2) = 3;
in_edit_sat4(:,2) = 4;

in_edit_sat1(:,3) = input_file(:,3);  
in_edit_sat2(:,3) = input_file(:,4);
in_edit_sat3(:,3) = input_file(:,5);
in_edit_sat4(:,3) = input_file(:,6);
in_edit_sat1(:,4) = input_file(:,7);
in_edit_sat2(:,4) = input_file(:,8);
in_edit_sat3(:,4) = input_file(:,9);
in_edit_sat4(:,4) = input_file(:,10);
in_edit_sat1(:,5) = input_file(:,11);
in_edit_sat2(:,5) = input_file(:,12);
in_edit_sat3(:,5) = input_file(:,13);
in_edit_sat4(:,5) = input_file(:,14);
in_edit_sat1(:,6) = input_file(:,15);
in_edit_sat2(:,6) = input_file(:,16);
in_edit_sat3(:,6) = input_file(:,17);
in_edit_sat4(:,6) = input_file(:,18);
in_edit_sat1(:,7) = input_file(:,19);
in_edit_sat2(:,7) = input_file(:,20);
in_edit_sat3(:,7) = input_file(:,21);
in_edit_sat4(:,7) = input_file(:,22);

in_edit = [in_edit_sat1;in_edit_sat2; in_edit_sat3; in_edit_sat4];
input_file(:,26:33) = [];
input_file(:,3:22) = [];

in_edit = sortrows(in_edit,1);

 
% Specify simulation times.
startTime = datetime(2021,6,24,8,0,0,"TimeZone","America/New_York");
simulationSteps = size(pos,1);
dt = 20;
time = zeros(size(pos,1),1);
time(1) = 0;
stopTime = startTime + seconds((simulationSteps-1)*dt);
for temp = 2:simulationSteps
time(temp) = time(temp-1) + dt;

end
duration = seconds(stopTime - startTime);
% Specify mask angle.
maskAngle = 10; % degrees

% Convert receiver position from north-east-down (NED) to geodetic
% coordinates.
receiverLLA = ned2lla(recPos,lla0,"ellipsoid");
receiverECEF = lla2ecef(receiverLLA);
% Specify the RINEX file.
rinexFile = "GODS00USA_R_20211750000_01D_GN.rnx";

% Set RNG seed to allow for repeatable results. 
rng("default");



% Create scenario.
sc = satelliteScenario(startTime, stopTime, dt);

% Initialize satellites. 
navmsg = rinexread(rinexFile);
satellite(sc,navmsg);
satID = sc.Satellites.Name;

% Preallocate results.
numSats = numel(sc.Satellites);
allSatPos = zeros(numSats,3,simulationSteps);
allSatVel = zeros(numSats,3,simulationSteps);

% Save satellite states over entire simulation.
for i = 1:numel(sc.Satellites)
    [oneSatPos, oneSatVel] = states(sc.Satellites(i),"CoordinateFrame","ecef");
    allSatPos(i,:,:) = permute(oneSatPos,[3 1 2]);
    allSatVel(i,:,:) = permute(oneSatVel,[3 1 2]);
end
% Preallocate results.
allP = zeros(numSats,simulationSteps);
allPDot = zeros(numSats,simulationSteps);
allIsSatVisible = false(numSats,simulationSteps);

% Use the skyplot to visualize satellites in view.
sp = skyplot([],[],MaskElevation=maskAngle);

for idx = 1:simulationSteps
    satPos = allSatPos(:,:,idx);
    satVel = allSatVel(:,:,idx);
    
    % Calculate satellite visibilities from receiver position.
    [satAz,satEl,allIsSatVisible(:,idx)] = lookangles(receiverLLA(idx,:),satPos,maskAngle);
    
    % Calculate pseudoranges and pseudorange rates using satellite and
    % receiver positions and velocities.
    [allP(:,idx),allPDot(:,idx)] = pseudoranges(receiverLLA(idx,:),satPos,recVel(idx,:),satVel);
    
    set(sp,"AzimuthData",satAz(allIsSatVisible(:,idx)), ...
        "ElevationData",satEl(allIsSatVisible(:,idx)), ...
        "LabelData",satID(allIsSatVisible(:,idx)))
    drawnow limitrate
 %   pause(0.02);
end

% Preallocate results.
lla = zeros(simulationSteps,3);
gnssVel = zeros(simulationSteps,3);
hdop = zeros(simulationSteps,1);
vdop = zeros(simulationSteps,1);

for idx = 1:simulationSteps
    p = allP(:,idx);
    pdot = allPDot(:,idx);
    isSatVisible = allIsSatVisible(:,idx);
    satPos = allSatPos(:,:,idx);
    satVel = allSatVel(:,:,idx);
    
    % Estimate receiver position and velocity using pseudoranges,
    % pseudorange rates, and satellite positions and velocities.
  %  [posECEF(idx,:),gnssVelECEF(idx,:),lla(idx,:),hdop(idx,:),vdop(idx,:)] = receiver_pos(p(isSatVisible), satPos(isSatVisible,:),pdot(isSatVisible),satVel(isSatVisible,:));


  refFrame = fusion.internal.frames.NED;

  initPosECEF = [0 0 0]; %we can put here position of the first satellite
initVelECEF = [0 0 0]; %check if we can out here previous receiver postion

%[posECEF(idx,:),gnssVelECEF(idx,:), dopMatrix] = computeLocation( ...
%    p(isSatVisible), satPos(isSatVisible,:),pdot(isSatVisible),satVel(isSatVisible,:), initPosECEF, initVelECEF);



posPrev = [initPosECEF(:); 0];
velPrev = [initVelECEF(:); 0];

posEst = posPrev;
velEst = velPrev;

Hpos = ones(size(satPos(isSatVisible,:), 1), 4);
Hvel = ones(size(satVel(isSatVisible,:), 1), 4);
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
    [pEst, ~, losVector] ...
        = nav.internal.gnss.calculatePseudoranges(satPos(isSatVisible,:), satVel(isSatVisible,:), ...
        posPrev(1:3).', velPrev(1:3).');
    % Add previous clock bias error (m). This is the time offset of the 
    % receiver clock times the speed of light.
    pPred = pEst + posPrev(end);
    
    Hpos(:,1:3) = -losVector;
    
    posEst = posPrev + Hpos.' * Hpos \ Hpos.' * (p(isSatVisible) - pPred);

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
   [~, pdotEst, losVector] = calculate_Pseudoranges(satPos(isSatVisible,:), satVel(isSatVisible,:),  posEst(1:3).', velPrev(1:3).');
    % Add previous clock drift error (m/s). This is the time drift of the 
    % receiver clock times the speed of light.
    pdotPred = pdotEst + velPrev(end);
    
    Hvel(:,1:3) = -losVector;
    
    velEst = velPrev + Hvel.' * Hvel \ Hvel.' * (pdot(isSatVisible) - pdotPred);
    
    resVel = norm(velEst - velPrev);
    
    iter = iter + 1;
    allResVel(iter) = resVel;
    if ~checkConverge(allResVel)
        velEst = velPrev;
        break;
    end
    velPrev = velEst;
end
posECEF(idx,:) = posEst(1:3).';

gnssVelECEF(idx,:) = velEst(1:3).';

dopMatrix = inv(Hpos.' * Hpos);




lla(idx,:) = fusion.internal.frames.ecef2lla(posECEF(idx,:));
gnssVel = refFrame.ecef2framev(gnssVelECEF, ...
    lla(idx,1), lla(idx,2));

% Convert DOP matrix from ECEF to local NAV frame and extract horizontal
% and vertical dilutions of precision.
%[hdop, vdop] = nav.internal.gnss.calculateDOP(dopMatrix, refFrame, lla);


end


estPos = lla2ned(lla,lla0,"ellipsoid");
figure
[Pos_LS, Pos_KF] = kalman_filtering(p, satPos);

difference = posECEF - receiverECEF; 
displacement = zeros(length(estPos),1);
for temp =2:size(estPos)
displacement(temp) = sqrt((posECEF(temp,1)- posECEF(temp-1,1))^2 + (posECEF(temp,2)- posECEF(temp-1,2))^2 +(posECEF(temp,3)- posECEF(temp-1,3))^2);

end
disp(["Number of visible satellites: ",  string(nnz(isSatVisible))]);
result = ["Successful compiled: ", string(datetime) ];
disp(result);

timingg = 0:dt:duration;
for temp = 1:length(estPos)
displacement(temp) = displacement(temp) + rand(1,1) * (-1).^round(rand(1,1))*5;

end
velocity_diff = zeros(length(estPos),1);
accelerometerde = zeros(length(estPos),1);
velocity_grad= zeros(length(estPos),1);
velocity_integ= zeros(length(estPos),1);
displacement_2 = zeros(length(estPos),1);
displacement_integ = zeros(length(estPos),1);
velocity_vel = zeros(length(estPos),1);
velocity_grad = gradient(displacement,time);

accelerometerde = gradient(velocity_grad,time);
for temp= 2:length(estPos)
velocity_integ = integral(@(t)accelerometerde,time(temp),time(temp-1),'ArrayValued',true);
end



displacement_integ = cumtrapz(dt,velocity_integ);
displacement_2 = cumtrapz(time,velocity_grad);
displacement_2 = [displacement_integ, displacement_2, displacement];

velocity_integ = [velocity_integ, velocity_grad] ;
%%
 % position of the reference station posECEF = should be true position from
 % orbit team
 dGnss_radius = 3000;
dGnss_X = zeros(2,1);
dGnss_Y = zeros(2,1);
dGnss_Z = zeros(2,1);
dGnss_XYZ = ones(2,2);
dGnss_distance = ones(2,1);
dGnss_distance = dGnss_distance * dGnss_radius;
integer =1;
 for temp = 4:length(posECEF)-5:length(posECEF)
    
 dGnss_X(integer) = posECEF(temp,1);
 dGnss_Y(integer) =  posECEF(temp,2);
dGnss_Z(integer) =  posECEF (temp,3);
integer= integer+1;
     
 end

dGnss_XYZ = [dGnss_X, dGnss_Y];

dGnss_location = dGnss_XYZ \ dGnss_distance;
dGnss_location = [dGnss_location; (dGnss_Z(2)+  dGnss_Z(2))/2];  


% clock bias satellite will be given by orbit team, receiver clock bias so
% far is t=0; we have to find epsilon_m l
%

% for i = 1:length(pos)
% delta_rho(i) = R_m(i) - rho_m(i)  
% rho_u(i) = p(i) + c*t + delta_rho(i);
% delta_rho =
% rho_u_corr(i) = sqrt()
% end
 
%pseudorange measurement, i ρm, to the ith satellite




%ρ = R + c dt + em // Kaplan Chapter10 this measurement contains the range to the satellite along with the errors
%The reference station differences the computed geometric range, im R , with the
%pseudorange measurement to form the differential correction

%Δρ(i) = R_m(i) − ρ_m(i) = −c*dt_m − ε_m

% This correction, which may be a positive or negative quantity, is broadcast to the
% user receiver where it is added to the user receiver%s pseudorange measurement to
% the same satellite

%ρho_u(i) + Δρ_m(i) = R_u(i) + c*dt_u +   (−c*dt_m − ε_m)

%Δρ(t)(i) = Δρ t + Δρ(i)*m(t)(t − t_m)     (12.8)
%The corrected user receiver pseudorange, i ( ρ t , for time t is then calculated from
% ρ_m(t)(i) = ρ(t)(i) + Δρ_m(t)
%ρu cor = Ru + εum + c dtum

%εum = εu – εm
%Pseudorange correction
%rho_u_corr(i) = sqrt(x-x^2 + y-y^2 + z-z^2) + em + cdt_m



