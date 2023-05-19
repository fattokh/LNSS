close all ;
%clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Implement alternative PVT computation block to target one (or more) moon more realistic navigation
%   scenarios. Examples:
%   1) The receiver is equipped with an altimeter.
%   2) The receiver is a rover exploring the Mo% on, equipped with an inertial measurement unit (IMU).
%   (Kalman Filter).
%   3) There is a Moon station transmitting a ranging signal to the receiver/rover (Differential
%       positioning).
%   To be described in the report
%       Discuss the performances of the developed PVT techniques, comparing the measured with the true
%       data, estimating
%       1) solution accuracy @ Kaplan 11.2.3
%       2) availability,    (@ Kaplan Chapter 11.3)
%       3) dilution of precision (DOP)  @.Kaplan 11.2.1
%
%
%   Also we can improve our code if we detect anomalies: 
%    - e.g. 92% cases are clock jumps and other clock anomalies. We can plot the result
%    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

settings = initSettings();
addpath include            
addpath geoFunctions 

gnssdata = rinexread("BRDC00WRD_R_20221260000_01D_EN.rnx").Galileo;

%% We get unique satellites (however, in the final version we include all dublicates as well)
[~,satIdx] = unique(gnssdata.SatelliteID);
gnssdata = gnssdata(satIdx,:);

%% Ephemeris preparation
eph.t_oc      	=	gnssdata.TransmissionTime;
eph.a_f2	=	gnssdata.SVClockDriftRate;
eph.a_f1 	=	gnssdata.SVClockDrift;
eph.a_f0 	=	gnssdata.SVClockBias;
eph.T_GD 	=	gnssdata.BRDCOrbit5Spare4;
eph.sqrtA 	=	gnssdata.sqrtA;
eph.t_oe 	=	gnssdata.Toe;
eph.deltan 	=	gnssdata.Delta_n;
eph.M_0 	=	gnssdata.M0;
eph.e	=	gnssdata.Eccentricity;
eph.omega  	=	gnssdata.omega;
eph.C_uc	=	gnssdata.Cuc;
eph.C_us 	=	gnssdata.Cus;
eph.C_rc	=	gnssdata.Crc;
eph.C_rs	=	gnssdata.Crs;
eph.i_0	=	gnssdata.i0;
eph.C_ic	=	gnssdata.Cic;
eph.C_is	=	gnssdata.Cis;
eph.omega_0	=	gnssdata.OMEGA0;
eph.omegaDot	=	gnssdata.OMEGA_DOT;
eph.iDot = gnssdata.IDOT;
eph.transmitTime = gnssdata.GALWeek;


%% Until now we do not have pseudorange and pseudorange rate. Therefore, we have dummy input data.
prnList = zeros(height(gnssdata),1);
prnList(:) = 24701805;
for i = 1:height(gnssdata)
prnList(i) = prnList(i) + rand(1,1)*51145709;
end
prnRate = zeros(height(gnssdata),3);
prnRate(:,:) = -500;
for i = 1:height(gnssdata)
prnRate(i,1) = prnRate(i,1) + rand(1,1)*1000;
prnRate(i,2) = prnRate(i,2) + rand(1,1)*1000;
prnRate(i,3) = prnRate(i,3) + rand(1,1)*1000;
end

%ls w/o alyitude


%% Satellite position and satellite clock correction
[satPositions, satClkCorr] = satposition(eph.transmitTime, prnList, eph, settings);

satpos = satPositions;
obs = prnRate;
%[pos, el, az,dop] = leastSquarePosition(satPositions, prnRate, settings);

nmbOfIterations = 7;


nmbOfSatellites = height(prnRate);
dtr     = pi/180;
pos     = zeros(4, 1);
X       = satpos;
A       = zeros(nmbOfSatellites, 4);
omc     = zeros(nmbOfSatellites, 1);
az      = zeros(1, nmbOfSatellites);
el      = az;

%=== Iteratively find receiver position ===================================
for iter = 1:nmbOfIterations

    for i = 1:nmbOfSatellites
        if iter == 1
            %--- Initialize variables at the first iteration --------------
            Rot_X = X(:, i);

            %% Variance covariance matrix including elevation angles for WLS @ Kaplan
            R = ones(25);
        else
            %--- Update equations -----------------------------------------
            rho2 = (X(1, i) - pos(1))^2 + (X(2, i) - pos(2))^2 + ...
                   (X(3, i) - pos(3))^2;
            traveltime = sqrt(rho2) / settings.c;

            %--- Correct satellite position (do to earth rotation) --------
            Rot_X = e_r_corr(traveltime, X(:, i));

            %--- Find the elevation angel of the satellite ----------------
            [az(i), el(i), dist] = topocent(pos(1:3, :), Rot_X - pos(1:3, :));
            
            for ii = 1:nmbOfSatellites
                 for j = 1:nmbOfSatellites
                      R (ii,j) = el(ii)*el(j);
                       
                 end
            end
            
        end % if iter == 1 ... ... else 

        %--- Apply the corrections ----------------------------------------
        omc(i) = (obs(i) - norm(Rot_X - pos(1:3), 'fro') + pos(4) );

        %--- Construct the A matrix ---------------------------------------
        A(i, :) =  [ ((Rot_X(1) - pos(1))) / obs(i), ...
                     ((Rot_X(2) - pos(2))) / obs(i), ...
                     ((Rot_X(3) - pos(3))) / obs(i), ...
                     1 ];
        
    end % for i = 1:nmbOfSatellites

    % These lines allow the code to exit gracefully in case of any errors
    if rank(A) ~= 4
        pos     = zeros(1, 4);
        
    end
Q       = inv(A'*A);
    %--- Find position update ---------------------------------------------


%% Full formula for WLS
    x   =A \ omc;%(A.' * inv(R) * A) \ A.' * inv(R)  * omc;
   %fd(iter,:) =x;
    %--- Apply position update --------------------------------------------
    pos = pos + x;
    
end % for iter = 1:nmbOfIterations

pos = pos';

%=== Calculate Dilution Of Precision ======================================

    %--- Initialize output ------------------------------------------------
dop     = zeros(1, 5);
    
    %--- Calculate DOP ----------------------------------------------------
    
    
    dop(1)  = sqrt(trace(Q));                       % GDOP    
    dop(2)  = sqrt(Q(1,1) + Q(2,2) + Q(3,3));       % PDOP
    dop(3)  = sqrt(Q(1,1) + Q(2,2));                % HDOP
    dop(4)  = sqrt(Q(3,3));                         % VDOP
    dop(5)  = sqrt(Q(4,4));                         % TDOP



skyplot(az,abs(el),satIdx);
poss = pos;
poss(3) = [];
LLA = xyz2lla(poss(1:3));


%{
p = input(6,:);
pdot = input(10,:);
satVel = input(11:13, :).';
satpos = satpos.';

%{
validateattributes(p, {'double', 'single'}, {'vector', 'real', 'finite'});

N = numel(p);
validateattributes(satpos, {'double', 'single'}, ...
    {'2d', 'nrows', N, 'ncols', 3, 'real', 'finite'});

validateattributes(pdot, {'double', 'single'}, ...
    {'vector', 'numel', N, 'real', 'finite'});

validateattributes(satVel, {'double', 'single'}, ...
    {'2d', 'nrows', N, 'ncols', 3, 'real', 'finite'});
%}
refFrame = fusion.internal.frames.NED;

initPosECEF = [0 0 0];
initVelECEF = [0 0 0];
[posECEF, gnssVelECEF, dopMatrix] = computeLocation(p, pdot, satpos, satVel, initPosECEF, initVelECEF);
%[Az, El, D] = topocent(X, dx);
lla = ecef2lla(posECEF);

%}
