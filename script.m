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


