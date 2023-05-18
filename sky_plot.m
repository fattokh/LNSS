close all;
clear all;
clc;

file = "input.csv";
input = readmatrix(file);
input(:,1) = [];
recPos = [42 -71 50];

maskAngle =5;
fileID = fopen('gnss_log_2023_05_18_00_56_32.txt');
fileIDD = rinexread("BRDC00WRD_R_20221260000_01D_EN.rnx");


gnssdata = fileIDD.Galileo;
[~,satIdx] = unique(gnssdata.SatelliteID);
gnssdata = gnssdata(satIdx,:);
% Read the text file
gpsData = fscanf(fileID,'%c');
parserObj = nmeaParser('MessageId','GSV');
% Parse the NMEA Data
gsvData = parserObj(gpsData);
az = []; el = []; satID = []; prevSentenceNumber = 0;
% Create an empty skyplot
sp = skyplot([], [],[]);
t = datetime('now','TimeZone','Local');
%[satPos,satVel,satID] = gnssconstellation(t );
%[satAz,satEl,slantRange] = lookangles(recPos, satPos,maskAngle); %ecef2aer(x,y,z,lat0,lon0,h0,wgs84);

t = gnssdata.Time(1);
[satPos,satVel,satID] = gnssconstellation(t,gnssdata);
prn =  timetable2table(gnssdata(:,1));
prn(:,1) = [];
prn = table2array(prn);
[satAz,satEl,slantRange] = lookangles(recPos, satPos,maskAngle);
skyplot(satAz,abs(satEl),prn) 


