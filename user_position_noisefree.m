close all
clc
clear

P1 = 1000;  % Pseudorange measurement to satellite 1
P2 = 1200;  % Pseudorange measurement to satellite 2
P3 = 900;  % Pseudorange measurement to satellite 3
P4 = 1100;  % Pseudorange measurement to satellite 4

% Satellite positions in ECEF coordinates
satellite_1_pos = [1000; 2000; 3000]; %[x1; y1; z1]; ECEF coordinates of satellite 1
satellite_2_pos = [-1500; 2500; -3500];%[x2; y2; z2];  % ECEF coordinates of satellite 2
satellite_3_pos = [2000; -3000; 4000];%[x3; y3; z3];  % ECEF coordinates of satellite 3
satellite_4_pos = [-2500; 1000; -1500];%[x4; y4; z4];  % ECEF coordinates of satellite 4

% Construct the matrix A and vector b for the linear system of equations
A = [
    2*(satellite_1_pos - satellite_2_pos)';
    2*(satellite_1_pos - satellite_3_pos)';
    2*(satellite_1_pos - satellite_4_pos)';
];
b = [
    P2^2 - P1^2 + norm(satellite_2_pos)^2 - norm(satellite_1_pos)^2;
    P3^2 - P1^2 + norm(satellite_3_pos)^2 - norm(satellite_1_pos)^2;
    P4^2 - P1^2 + norm(satellite_4_pos)^2 - norm(satellite_1_pos)^2;
];

% Solve the linear system of equations using least squares method
user_pos = pinv(A) * b;

% Extract the user position coordinates
x = user_pos(1);
y = user_pos(2);
z = user_pos(3);

disp(user_pos);
settings = initSettings();

[pos, el, az, dop] = leastSquarePos(satpos, obs, settings);