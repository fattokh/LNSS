close all ;
%clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Implement alternative PVT computation block to target one (or more) moon more realistic navigation
%   scenarios. Examples:
%   1) The receiver is equipped with an altimeter.
%         [altimeter] = We have input data
%   2) The receiver is a rover exploring the Moon, equipped with an inertial measurement unit (IMU).
%   (Kalman Filter). 
%
%        
%            [accelerometer] = calculate distance, we need double integration of acceleration,
%            but we will have noise
%            [magnetometer] = to identify the north direction 
%            [gyroscope] = to identify change in orientation or rotation
%            [Kalman] = to estimate receiver position, smoothing the
%            tracking (or estimate when there is temporary GNSS
%            unavailability
%              
%
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
prnList(i) = prnList(i) + rand(1,1)*5002709;
end
satVel = zeros(height(gnssdata),3);
satVel(:,:) = -3500;
prnRate = zeros(height(gnssdata),1);
prnRate(:,:) = -500+rand(1,1)*1000;
for i = 1:height(gnssdata)
satVel(i,1) = satVel(i,1) + rand(1,1)*7000;
satVel(i,2) = satVel(i,2) + rand(1,1)*7000;
satVel(i,3) = satVel(i,3) + rand(1,1)*7000;
end

%ls w/o alyitude    


%% Satellite position and satellite clock correction
[satPositions, satClkCorr] = satposition(eph.transmitTime, prnList, eph, settings);

satpos = satPositions;
satpos1 = satPositions.';
obs = prnRate;
%[pos, el, az,dop] = leastSquarePosition(satPositions, prnRate, settings);
nmbOfIterations = 7;

dtr     = pi/180;
pos     = zeros(3, 1);
X       = satpos;
nmbOfSatellites = size(satpos, 2);

A       = zeros(nmbOfSatellites, 3);
omc     = zeros(nmbOfSatellites, 1);
az      = zeros(1, nmbOfSatellites);
el      = az;

%=== Iteratively find receiver position ===================================
for iter = 1:nmbOfIterations

    for i = 1:nmbOfSatellites
        if iter == 1
            %--- Initialize variables at the first iteration --------------
            Rot_X = X(:, i);
            trop = 2;
        else
            %--- Update equations -----------------------------------------
            rho2 = (X(1, i) - pos(1))^2 + (X(2, i) - pos(2))^2 + ...
                   (X(3, i) - pos(3))^2;
            traveltime = sqrt(rho2) / settings.c ;

            %--- Correct satellite position (do to earth rotation) --------
            Rot_X = e_r_corr(traveltime, X(:, i));

            %--- Find the elevation angel of the satellite ----------------
            [az(i), el(i), dist] = topocent(pos(1:3, :), Rot_X - pos(1:3, :));

            if (settings.useTropCorr == 1)
                %--- Calculate tropospheric correction --------------------
                trop = tropo(sin(el(i) * dtr), ...
                             0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
            else
                % Do not calculate or apply the tropospheric corrections
                trop = 0;
            end
        end % if iter == 1 ... ... else 

        %--- Apply the corrections ----------------------------------------
        omc(i) = (obs(i) - norm(Rot_X - pos(1:3), 'fro') );

        %--- Construct the A matrix ---------------------------------------
        A(i, :) =  [ (-(Rot_X(1) - pos(1))) / obs(i); ...
                     (-(Rot_X(2) - pos(2))) / obs(i); ...
                     (-(Rot_X(3) - pos(3))) / obs(i)   ];
    end % for i = 1:nmbOfSatellites

    % These lines allow the code to exit gracefully in case of any errors


    %--- Find position update ---------------------------------------------
    x   = A \ omc;
    
    %--- Apply position update --------------------------------------------
    pos = pos + x;
    
end % for iter = 1:nmbOfIterations

pos = pos';

%=== Calculate Dilution Of Precision ======================================

    %--- Initialize output ------------------------------------------------
    dop     = zeros(1, 5);
    
    %--- Calculate DOP ----------------------------------------------------
    Q       = inv(A'*A);
    
    dop(1)  = sqrt(trace(Q));                       % GDOP    
    dop(2)  = sqrt(Q(1,1) + Q(2,2) + Q(3,3));       % PDOP
    dop(3)  = sqrt(Q(1,1) + Q(2,2));                % HDOP
    dop(4)  = sqrt(Q(3,3));                         % VDOP
 %   dop(5)  = sqrt(Q(4,4));                         % TDOP



[lla1, gnssVel1, hdop1, vdop1] = receiverposition(prnList,satpos1,prnRate,satVel);
%geoplot(lla1(1),lla1(2));
%geobasemap topographic
%skyplot(az,abs(el),satIdx);
poss = pos;

LLA = xyz2lla(poss(1:3));


    c=imread("Lunar_map.jpg");
    center = [0 0 0]; % center of the sphere
  radius = 1737;
    % Initalize longitude and colatitude range
    longitude = deg2rad(linspace(0,360,size(c,2)));
    colatitude = deg2rad(linspace(90,-90,size(c,1)));
    
    circle_center = [0.5 0.5 0.5]; % center of the circle
circle_radius = 0.3; % radius of the circle
theta_circle = linspace(0, 2*pi, 100); % azimuthal angle for the circle
phi_circle = linspace(0, pi, 50); % polar angle for the circle
[theta_circle, phi_circle] = meshgrid(theta_circle, phi_circle);

x_circle = circle_center(1) + circle_radius*sin(phi_circle).*cos(theta_circle);
y_circle = circle_center(2) + circle_radius*sin(phi_circle).*sin(theta_circle);
z_circle = circle_center(3) + circle_radius*cos(phi_circle);


    % Initialize empty variables
    X = zeros(size(c,1),size(c,2));
    Y = zeros(size(c,1),size(c,2));
    Z = zeros(size(c,1),size(c,2));
    % Convert spherical to cartesian coordinates
    for j = 1:length(colatitude)
        for i = 1:length(longitude)
            [X(j,i), Y(j,i), Z(j,i)] = sph2cart( ...
                                            longitude(i),...
                                            colatitude(j),...
                                            radius);
        end
    end
n=30;
r = 10; % radius of the sphere

thetha = 0:pi/(n/2):2*pi; 
phi    = -pi:2*pi/n:pi;
xp     = r.*sin(phi).*cos(thetha);
yp     = r.*sin(thetha).*sin(phi);
zp     = r.*cos(phi);

    figure

    hold on 

    scatter3(X(100),Y(100),Z(100));
    hold on
        surf(X,Y,Z,c,'EdgeColor','none','FaceColor','texturemap');

    axis equal


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
