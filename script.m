close all;
clear all;
clc;
data = rinexread("BRDC00WRD_R_20231370000_01D_CN.rnx");
file = "input.csv";
input = readmatrix(file);
settings = initSettings();
addpath include            
addpath geoFunctions 
input(:,1) = [];
ephemeris.prn = input(1,:);
ephemeris.tow = input(2,:);
ephemeris.week = input(3,:);
ephemeris.rx_time = input(4,:) ;
ephemeris.user_clk_offset = input(5,:);
ephemeris.pseudorange= input(6,:);
ephemeris.pos_x = input(7,:);
ephemeris.pos_y = input(8,:);
ephemeris.pos_z = input(9,:);
ephemeris.pseudorange_rate = input(10,:);
ephemeris.vel_x = input(11,:);
ephemeris.vel_y = input(12,:);
ephemeris.vel_z = input(13,:);



satpos = zeros(3,size(input,2));
obs = zeros(1, size(input,2));
receiver_position = [0;0;0];


for k = 1:size(input,2)
satpos(:,k) = input(7:9,k);
obs(:,k) = input(6,k);
end



c_0 = 299792458;

[pos, el, az,dop] = leastSquarePosition(satpos, obs, settings);

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

lla = ecef2lla(posECEF);
xo =zeros(8,1);
prs = p;
gpsEph = pdot;

%[xHat,z,svPos,H,Wpr,Wrr] = WeightedLeastSquares(p,pdot,xo);
jWk=1; jSec=2; jSv=3; jPr=4; jPrSig=5; jPrr=6; jPrrSig=7;%index of columns


xHat=[]; z=[]; H=[]; svPos=[];
xyz0 = xo(1:3);
bc = xo(4);
%error = "error";
numVal=4;
if numVal<4
  error("error");
end
ttxWeek = prs(:,jWk); %week of tx. Note - we could get a rollover, when ttx_sv
%goes negative, and it is handled in GpsEph2Pvt, where we work with fct
ttxSeconds =  prs(:,jSec) - prs(:,jPr)/settings.c; %ttx by sv clock 
% this is accurate satellite time of tx, because we use actual pseudo-ranges 
% here, not corrected ranges
% write the equation for pseudorange to see the rx clock error exactly cancel
% to get precise GPS time: we subtract the satellite clock error from sv time, 
% as done next:
 %make into a column for compatibility with other time vectors
ttx = ttxSeconds - p/settings.c; %subtract dtsv from sv time to get true gps time

%calculate satellite position at ttx
%[svXyzTtx,dtsv,svXyzDot,dtsvDot]=GpsEph2Pvt(ephemeris,[ttxWeek,ttx]);

xyzM=[]; dtsvS=[]; 

% [bOk,gpsEph,gpsWeek,ttxSec] = CheckGpsEphInputs(gpsEph,[ttxWeek,ttx]);

gpsTime = [ttxWeek,ttx];
gpsWeek = gpsTime(:,1);
ttxSec  = gpsTime(:,2);



p=length(gpsEph);
%Now we are done checking and manipulating the inputs
% the time vectors: gpsWeek, ttxSec are the same length as gpsEph
%-------------------------------------------------------------------------------

%set fitIntervalSeconds
fitIntervalHours = [gpsEph.Fit_interval]';
%Rinex says "Zero if not known", so adjust for zeros
fitIntervalHours(fitIntervalHours == 0) = 2;
fitIntervalSeconds = fitIntervalHours*3600; 

%Extract variables from gpsEph, into column vectors
%orbit variable names follow RINEX 2.1 Table A4
%clock variables af0, af1, af2 follow IS GPS 200
TGD     = [gpsEph.TGD]';
Toc     = [gpsEph.Toc]';
af2     = [gpsEph.af2]';
af1     = [gpsEph.af1]';
af0     = [gpsEph.af0]';
Crs     = [gpsEph.Crs]';
Delta_n = [gpsEph.Delta_n]';
M0      = [gpsEph.M0]';
Cuc     = [gpsEph.Cuc]';
e       = [gpsEph.e]';
Cus     = [gpsEph.Cus]';
Asqrt   = [gpsEph.Asqrt]';
Toe     = [gpsEph.Toe]';
Cic     = [gpsEph.Cic]';
OMEGA   = [gpsEph.OMEGA]';
Cis     = [gpsEph.Cis]';
i0      = [gpsEph.i0]';
Crc     = [gpsEph.Crc]';
omega   = [gpsEph.omega]';
OMEGA_DOT=[gpsEph.OMEGA_DOT]';  
IDOT    = [gpsEph.IDOT]';
ephGpsWeek = [gpsEph.GPS_Week]';
    
%Calculate dependent variables -------------------------------------------------

%Time since time of applicability accounting for weeks and therefore rollovers
%subtract weeks first, to avoid precision errors:
tk =(gpsWeek-ephGpsWeek)*GpsConstants.WEEKSEC + (ttxSec-Toe);

I = find(abs(tk)>fitIntervalSeconds);
if ~isempty(I)
    numTimes = length(I);
    fprintf(sprintf('WARNING in GpsEph2Xyz.m, %d times outside fit interval.',...
        numTimes));
end

A = Asqrt.^2;    %semi-major axis of orbit
n0=sqrt(GpsConstants.mu./(A.^3));    %Computed mean motion (rad/sec)
n=n0+Delta_n;      %Corrected Mean Motion
h = sqrt(A.*(1-e.^2).*GpsConstants.mu); 
Mk=M0+n.*tk;     %Mean Anomaly
Ek=Kepler(Mk,e);  %Solve Kepler's equation for eccentric anomaly

%Calculate satellite clock bias (See ICD-GPS-200 20.3.3.3.3.1)
%subtract weeks first, to avoid precision errors:
dt =(gpsWeek-ephGpsWeek)*GpsConstants.WEEKSEC + (ttxSec-Toc);

%Calculate satellite clock bias
sin_Ek=sin(Ek);
cos_Ek=cos(Ek);
dtsvS = af0 + af1.*dt + af2.*(dt.^2)  + ...
    GpsConstants.FREL.*e.*Asqrt.*sin_Ek -TGD;

%true anomaly:
vk=atan2(sqrt(1-e.^2).*sin_Ek./(1-e.*cos_Ek),(cos_Ek-e)./(1-e.*cos_Ek));
Phik=vk + omega;   %Argument of latitude          

sin_2Phik=sin(2*Phik);
cos_2Phik=cos(2*Phik);
% The next three terms are the second harmonic perturbations
duk = Cus.*sin_2Phik + Cuc.*cos_2Phik;   %Argument of latitude correction
drk = Crc.*cos_2Phik + Crs.*sin_2Phik;   %Radius Correction
dik = Cic.*cos_2Phik + Cis.*sin_2Phik;   %Correction to Inclination

uk = Phik + duk;  %Corrected argument of latitude
rk = A.*((1-e.^2)./(1+e.*cos(vk))) + drk; %Corrected radius
ik = i0 + IDOT.*tk + dik; %Corrected inclination

sin_uk=sin(uk);
cos_uk=cos(uk);
xkp = rk.*cos_uk; % Position in orbital plane
ykp = rk.*sin_uk; % Position in orbital plane

% Wk = corrected longitude of ascending node, 
Wk  = OMEGA + (OMEGA_DOT - GpsConstants.WE).*tk - GpsConstants.WE*[gpsEph.Toe]';

%for dtflight, see FlightTimeCorrection.m
sin_Wk = sin(Wk);
cos_Wk = cos(Wk);

xyzM=zeros(p,3);
% The matrix xyzM contains the ECEF coordinates of position
sin_ik=sin(ik);
cos_ik=cos(ik);
xyzM(:,1)=xkp.*cos_Wk-ykp.*cos_ik.*sin_Wk;
xyzM(:,2)=xkp.*sin_Wk+ykp.*cos_ik.*cos_Wk;
xyzM(:,3)=ykp.*sin_ik; 



svXyzTrx = svXyzTtx; %initialize svXyz at time of reception

%Compute weights ---------------------------------------------------
Wpr = diag(1./prs(:,jPrSig));
Wrr = diag(1./prs(:,jPrrSig));

%iterate on this next part tilL change in pos & line of sight vectors converge
xHat=zeros(4,1);
dx=xHat+inf;
whileCount=0; maxWhileCount=100; 
%we expect the while loop to converge in < 10 iterations, even with initial
%position on other side of the Earth (see Stanford course AA272C "Intro to GPS")
while norm(dx) > GnssThresholds.MAXDELPOSFORNAVM
    whileCount=whileCount+1;
    assert(whileCount < maxWhileCount,...
        'while loop did not converge after %d iterations',whileCount);
    for i=1:length(gpsEph)
        % calculate tflight from, bc and dtsv
        dtflight = (prs(i,jPr)-bc)/GpsConstants.LIGHTSPEED + dtsv(i);
        % Use of bc: bc>0 <=> pr too big <=> tflight too big.
        %   i.e. trx = trxu - bc/GpsConstants.LIGHTSPEED
        % Use of dtsv: dtsv>0 <=> pr too small <=> tflight too small.
        %   i.e ttx = ttxsv - dtsv
        svXyzTrx(i,:) = FlightTimeCorrection(svXyzTtx(i,:), dtflight);
    end

  %calculate line of sight vectors and ranges from satellite to xo
  v = xyz0(:)*ones(1,numVal,1) - svXyzTrx';%v(:,i) = vector from sv(i) to xyz0
  range = sqrt( sum(v.^2) );
  v = v./(ones(3,1)*range); % line of sight unit vectors from sv to xo

  svPos=[prs(:,3),svXyzTrx,dtsv(:)];

  %calculate the a-priori range residual
  prHat = range(:) + bc -GpsConstants.LIGHTSPEED*dtsv;
  % Use of bc: bc>0 <=> pr too big <=> rangehat too big
  % Use of dtsv: dtsv>0 <=> pr too small
    
  zPr = prs(:,jPr)-prHat; 
  H = [v', ones(numVal,1)]; % H matrix = [unit vector,1]
  
  %z = Hx, premultiply by W: Wz = WHx, and solve for x:
  dx = pinv(Wpr*H)*Wpr*zPr;

  % update xo, xhat and bc
  xHat=xHat+dx;
  xyz0=xyz0(:)+dx(1:3);
  bc=bc+dx(4);

  %Now calculate the a-posteriori range residual
  zPr = zPr-H*dx;
end

% Compute velocities ---------------------------------------------------------
rrMps = zeros(numVal,1);
for i=1:numVal
    %range rate = [satellite velocity] dot product [los from xo to sv]
    rrMps(i) = -svXyzDot(i,:)*v(:,i);
end
prrHat = rrMps + xo(8) - GpsConstants.LIGHTSPEED*dtsvDot;
zPrr = prs(:,jPrr)-prrHat;
%z = Hx, premultiply by W: Wz = WHx, and solve for x:
vHat = pinv(Wrr*H)*Wrr*zPrr;
xHat = [xHat;vHat]; 

z = [zPr;zPrr];