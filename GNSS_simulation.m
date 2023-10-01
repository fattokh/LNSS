%% Readme

% this script simulates a gps receiver on the Earth. We use the last version
% of the file "data_Rx_group.csv" we received from the receiver gruop and
% a real gps data file that we have downloaded on the Internet because we
% were not provided with the navigation message, which is necessary to
% compute the orbits.

% we have made this script to test our implementation in an easier enviroment like
% Earth, before adapting it to the Moon.

% As output we peovide:
%   - figure 1: Satellite positions
%   - figure 2: Position of the receiver
%   - figure 3: Position (ECEF) Error fo each axis
%   - figure 4: Accelerometer Readings
%   - figure 5: Position (ECEF) Error
%   - figure 6: Dilution of Precision (DOP)
%   - figure 7: Number of visible satellites
%   - figure 8: Satellite visibility

% related to figure 1, if you do not want to observe the satellite movements
% comment the line 172 (pause(0.05);).


%% Reading the initial data

close all;
clearvars;
clc;

% % simulated user trajectory (from receiver group)
% % format: latitude, longitude, position in xyz, velocity along x, y, z
initial_data = readmatrix("data_Rx_group.csv");
% extract the data
% receiver position data
rec_pos = initial_data(:,3:5);
% receiver velocity data
rec_vel = initial_data(:,6:8);

% this file contains the satellite parameters: orbits, eccentricity,
% positions, pseudoranges and so on (from orbital team).
orbital_data = readmatrix("data_Orbit_group.csv");
% remove first row with the columns' names
orbital_data(1,:) = [];

% vectors with information about each satellite
% format: sat ID, simulation time, X, Y, Z, distance, pseudoranges.
sat1 = [ones(height(orbital_data),1)*1, orbital_data(:,1),orbital_data(:,3), orbital_data(:,7),orbital_data(:,11), orbital_data(:,15),  orbital_data(:,19)];
sat2 = [ones(height(orbital_data),1)*2, orbital_data(:,1),orbital_data(:,4), orbital_data(:,8),orbital_data(:,12), orbital_data(:,16),  orbital_data(:,20)];
sat3 = [ones(height(orbital_data),1)*3, orbital_data(:,1),orbital_data(:,5), orbital_data(:,9),orbital_data(:,13), orbital_data(:,17), orbital_data(:,21)];
sat4 = [ones(height(orbital_data),1)*4, orbital_data(:,1),orbital_data(:,6),orbital_data(:,10),orbital_data(:,14), orbital_data(:,18) , orbital_data(:,22)];

% specify simulation times as indicated on the orbit csv
start_time = datetime(2022,01,06,8,0,0,"TimeZone","America/New_York");
% number of samples
steps = size(rec_pos,1);
% time interval
dt = 25;
stop_time = start_time + seconds((steps-1)*dt);
duration = seconds(stop_time - start_time);
% mask angle [degrees], specified by orbit team
mask_angle = 10;

% initial guess in NED geodetic coordinates
lla_guess = [42.31451, -71.38236, 0];
% receiver position conversion from north-east-down coordinates to
% geodetic coordinates and ECEF
receiver_LLA = ned2lla(rec_pos,lla_guess,"ellipsoid");
receiver_ECEF = lla2ecef(receiver_LLA);


%% Calculate Recevier parameters

% calculate the distance between the receiver positions in two
% successive time instants
dist_receiver_x = zeros(steps-1,1);
dist_receiver_y=  zeros(steps-1,1);
dist_receiver_z = zeros(steps-1,1);
for temp = 1:steps-1
    dist_receiver_x(temp)= receiver_ECEF(temp+1,1) - receiver_ECEF(temp,1);
    dist_receiver_y(temp) = receiver_ECEF(temp+1,2) - receiver_ECEF(temp,2);
    dist_receiver_z(temp) = receiver_ECEF(temp+1,3) - receiver_ECEF(temp,3);
end

% time vector 1 to steps-1
time_input = linspace(1,steps-1,steps-1);

% velocity as derivative of space wrt the time
velocity_x = gradient(dist_receiver_x,time_input);
velocity_y = gradient(dist_receiver_y,time_input);
velocity_z = gradient(dist_receiver_z,time_input);

% calculate acceleration as derivative of velocity wrt time
acceleration_x = gradient(velocity_x,time_input);
acceleration_y = gradient(velocity_y,time_input);
acceleration_z = gradient(velocity_z,time_input);

% create matrices of distance and acceleration of receiver along
% the three axis
distance = [dist_receiver_x, dist_receiver_y, dist_receiver_z];
acceleration = [acceleration_x,acceleration_y,acceleration_z];


%% Satellite simulation

% this is an example of GPS navigation message file
gps_data = "real_gps_data.rnx";
navmsg = rinexread(gps_data);
% simulate the satellites in the sky
scenario = satelliteScenario(start_time, stop_time, dt);
% adding satellites from the rinex file
satellite(scenario,navmsg);
sat_ID = scenario.Satellites.Name;
num_sats = numel(scenario.Satellites);

% initializing satellites' positions and velocities
sat_positions = zeros(num_sats,3,steps);
sat_velocities = zeros(num_sats,3,steps);
for i = 1:numel(scenario.Satellites)
    [single_sat_pos, single_sat_vel] = states(scenario.Satellites(i),"CoordinateFrame","ecef");
    % fixing the order of parameters
    sat_positions(i,:,:) = permute(single_sat_pos, [3 1 2]);
    sat_velocities(i,:,:) = permute(single_sat_vel, [3 1 2]);
end

% vectors with all pseudoranges and pseudorange rates
all_pseudoranges = zeros(num_sats,steps);
all_pseudorange_rates = zeros(num_sats,steps);
% pseudoranges differential (station point of view)
all_pseudoranges_diff = zeros(num_sats,steps);
% vector with visibility information
all_are_visible = false(num_sats,steps);

% plot the position of satellites in the sky, removing the ones that are
% located in the lower 10° because we assume that they are too low
figure(1)
sp = skyplot([],[],MaskElevation=mask_angle);
title('Satellites Positions')


%% Simulating the ground station - DGNSS

pseudo_sat_dgnss = zeros(num_sats,steps);
% let's place a reference station near the receiver
station_position = [mean(receiver_LLA(:,1:2)) , mean(receiver_LLA(:,3)+30)];
for idx = 1:steps
    % take postion and velocity of each satellite at istant idx
    single_sat_pos = sat_positions(:,:,idx);
    single_sat_vel = sat_velocities(:,:,idx);

    % returns the look angles azimuth and elevation in degrees,
    % of the satellites using the satellite positions single_sat_pos in the
    % ECEF coordinate system in meters and the receiver position in
    % geodetic coordinates.
    % the output is an array of azimuth and elevation angles and a logical
    % array specifying the visibility of each satellite.
    [satAz,satEl,all_are_visible(:,idx)] = lookangles(receiver_LLA(idx,:),single_sat_pos,mask_angle);

    % returns the pseudorange and pseudorange rates, pdot, in meters per second,
    % between the receiver and satellites.
    [all_pseudoranges(:,idx),all_pseudorange_rates(:,idx)] = pseudoranges(receiver_LLA(idx,:),single_sat_pos,rec_vel(idx,:),single_sat_vel);
    % computing the pseudoranges but from station point of view
    all_pseudoranges_diff(:,idx) = pseudoranges(station_position,single_sat_pos);
    % manually computing the distance for each satellite to the station
    for ii =1:num_sats
        pseudo_sat_dgnss(ii,idx) = sqrt((station_position(1)-single_sat_pos(ii,1))^2+(station_position(2)-single_sat_pos(ii,2))^2+(station_position(3)-single_sat_pos(ii,3))^2);
    end
    set(sp,"AzimuthData",satAz(all_are_visible(:,idx)), "ElevationData",satEl(all_are_visible(:,idx)), "LabelData",sat_ID(all_are_visible(:,idx)))

    % if you do not want to observe the satellite movements comment this line
    pause(0.05);
end


%% DOP and position

% initializing latitude, longitude, and altitude coordinates
lla = zeros(steps,3);
% latitude, longitude, and altitude coordinates with altimeter
lla_alt = zeros(steps,3);
% initializing DOP parameters
hdop = zeros(steps,1);
vdop = zeros(steps,1);
% initializing DOP parameters for differential
hdop_diff = zeros(steps,1);
vdop_diff = zeros(steps,1);
% initializing DOP parameters with altimeter
hdop_alt = zeros(steps,1);
vdop_alt = zeros(steps,1);

% taking position x,y,z at first time instant for a satellite (ECEF)
init_sat_pos = [sat_positions(1,1,1)/10 sat_positions(1,2,1)/10 sat_positions(1,3,1)/10] ;
% previous positions (with time index)
pos_prev = [init_sat_pos(:); 0];
pos_prev_diff = [init_sat_pos(:); 0];
% set the velocity to zero in x, y, z at instant zero
vel_prev = [0; 0; 0; 0];
pos_ECEF = zeros(steps, 3);
pos_ECEF_diff = zeros(steps, 3);
% set the desired frame
refFrame = fusion.internal.frames.NED;
% preallocation
lla_diff = zeros(steps,3);

for idx = 1:steps
    % copy pseudoranges at step idx
    pseudoranges_idx = all_pseudoranges(:,idx);
    % take only the current values at instant idx
    pdot = all_pseudorange_rates(:,idx);
    is_visible = all_are_visible(:,idx);
    single_sat_pos = sat_positions(:,:,idx);
    single_sat_vel = sat_velocities(:,:,idx);

    % position estimates for the loop
    pos_est = pos_prev;
    pos_est_diff = pos_prev_diff;

    % H matrices
    H_pos = ones(size(single_sat_pos(is_visible,:), 1), 4);
    H_pos_diff = ones(size(single_sat_pos(is_visible,:), 1), 4);
    H_pos_prev = H_pos;

    % we begin the computation of the position: initially, the resolution
    % is infinite (max uncertainty) and we want to achieve the resolution
    % specified in the variable target_res with the successive iterations.
    res_pos = Inf;
    res_pos_diff = Inf;
    target_res = 1;
    max_iter = 200;
    iter = 0;
    all_res_pos = NaN(max_iter, 1);
    all_res_pos_diff = NaN(max_iter, 1);

    % loop until either we pass max iterations count or arrive to the
    % desired resolution
    while (res_pos > target_res) && (iter < max_iter)
        % p_est: pseudorange for all visible satellites
        % los_vector: line-of-sight unit vectors from receiver to satellite
        [p_est, ~, los_vector] = nav.internal.gnss.calculatePseudoranges(single_sat_pos(is_visible,:), single_sat_vel(is_visible,:), pos_prev(1:3).', vel_prev(1:3).');
        % computing the prediction
        p_pred = p_est + pos_prev(end);
        H_pos(:,1:3) = - los_vector;
        % calculate the next position estimation (state update equation)
        pos_est = pos_prev + H_pos.' * H_pos \ H_pos.' * (pseudoranges_idx(is_visible) - p_pred);
        res_pos = norm(pos_est - pos_prev);
        % increment iteration count
        iter = iter + 1;
        % save the current result
        all_res_pos(iter) = res_pos;
        % set the current results as previous so we can proceed to the
        % next iteration
        pos_prev = pos_est;
        H_pos_prev = H_pos;
    end

    % second loop for differential computation
    % loop until either we pass max iterations count or arrive to the
    % desired resolution
    iter = 0;
    while (res_pos_diff > target_res) && (iter < max_iter)
        % dgnss position accounting for the ground station
        [pEst_dgnss, ~,~] = nav.internal.gnss.calculatePseudoranges(single_sat_pos(is_visible,:), single_sat_vel(is_visible,:), station_position, [0,0,0]);
        delta_rho = pEst_dgnss-pseudo_sat_dgnss(is_visible,idx);
        % pseudoranges with previous estimated position
        [p_est_diff, ~, los_vector_diff] = nav.internal.gnss.calculatePseudoranges(single_sat_pos(is_visible,:), single_sat_vel(is_visible,:), pos_prev_diff(1:3).', [0,0,0]);
        % new orbit estimate adding Least Square deviation (lecture 21, slide 27, point 3)
        pPred_1 = p_est_diff + delta_rho;
        H_pos_diff(:,1:3) = -los_vector_diff;
        % calculate the next position estimation (state update equation)
        pos_est_diff = pos_prev_diff + H_pos_diff.' * H_pos_diff \ H_pos_diff.' * (pseudoranges_idx(is_visible) - pPred_1);
        % compute the resolution of the current position estimate
        res_pos_diff = norm(pos_est_diff - pos_prev_diff);
        % update iteration and save the result
        iter = iter + 1;
        all_res_pos_diff(iter) = res_pos_diff;
        % set the current results as previous so we can proceed to the next iteration
        pos_prev_diff = pos_est_diff;
        HposPrev_1 = H_pos_diff;
    end

    % calculate the dop parameters and formatting the results for both the
    % normal measurements and differential
    pos_ECEF(idx,:) = pos_est(1:3).';
    pos_ECEF_diff(idx,:) = pos_est_diff(1:3).';

    dop_matrix = inv(H_pos.' * H_pos);
    dop_matrix_diff = inv(H_pos_diff.' * H_pos_diff);

    % convert the results in lat long alt for normal, differential and altimeter
    lla(idx,:) = ecef2lla(pos_ECEF(idx,:));
    lla_diff(idx,:) = ecef2lla(pos_ECEF_diff(idx,:));
    lla_alt(idx,:)= [lla_diff(idx,1:2),receiver_LLA(idx,3) + (randn(1,1)*0.2-0.1)];

    [hdop(idx), vdop(idx)] = calculate_DOP(dop_matrix, refFrame, lla);
    [hdop_diff(idx), vdop_diff(idx)] = calculate_DOP(dop_matrix_diff, refFrame, lla_diff);
    [hdop_alt(idx), vdop_alt(idx)] = calculate_DOP(dop_matrix_diff, refFrame, lla_alt);
end


%% Results plotting

% 3D plot of receiver position
figure(2)
plot3(pos_ECEF(:,1)/1000,pos_ECEF(:,2)/1000,pos_ECEF(:,3)/1000);
legend("user trajectory")
xlabel("X - axis position [km]")
ylabel("Y - axis position [km]")
zlabel("Z - axis position [km]")
title("Position of the receiver")
grid on

% ECEF position errors with and without altimeter
figure(3)
% convert position to north east down system
estPos = lla2ned(lla,lla_guess,"ellipsoid");
winSize = floor(size(estPos,1)/10);
alt_ECEF= lla2ecef(lla_alt);
subplot(3,1,1)
plot(smoothdata(abs(pos_ECEF(:,1)-receiver_ECEF(:,1)),"movmedian",winSize))
hold on
plot(smoothdata(abs(alt_ECEF(:,1)-receiver_ECEF(:,1)),"movmedian",winSize))
legend("without altimeter","with altimeter")
xlabel("Time [s]")
ylabel("Error X [m]")
xlim([1,steps])
title("Position (ECEF) Error")

subplot(3,1,2)
plot(smoothdata(abs(pos_ECEF(:,2)-receiver_ECEF(:,2)),"movmedian",winSize))
hold on
plot(smoothdata(abs(alt_ECEF(:,2)-receiver_ECEF(:,2)),"movmedian",winSize))
legend("without altimeter","with altimeter")
xlabel("Time [s]")
ylabel("Error Y [m]")
xlim([1,steps])

subplot(3,1,3)
plot(smoothdata(abs(pos_ECEF(:,3)-receiver_ECEF(:,3)),"movmedian",winSize))
hold on
plot(smoothdata(abs(alt_ECEF(:,3)-receiver_ECEF(:,3)),"movmedian",winSize))
legend("without altimeter","with altimeter")
xlabel("Time [s]")
ylabel("Error Z [m]")
xlim([1,steps])

% Accelerometer Readings
figure(4)
t = 1:steps-1;
subplot(3,1,1)

plot(t,acceleration(:,1)+randn(steps-1, 1)*(-1).^(round(randn(1,1))),'Color',[0 0.4470 0.7410])
hold on
legend('X-axis')
xlabel("Time [s]")
ylabel('Acceleration [m/s^2]')
xlim([1,steps])
title('Accelerometer Readings')

subplot(3,1,2)
plot(t,acceleration(:,2)+randn(steps-1, 1)*(-1).^(round(randn(1,1))),'Color',[0.8500 0.3250 0.0980])
hold on
legend('Y-axis')
xlim([1,steps])
xlabel("Time [s]")
ylabel('Acceleration [m/s^2]')

subplot(3,1,3)
plot(t,acceleration(:,3)+randn(steps-1, 1)*(-1).^(round(randn(1,1))),'Color',[0.9290 0.6940 0.1250])
hold off
legend('Z-axis')
xlim([1,steps])
xlabel("Time [s]")
ylabel('Acceleration [m/s^2]')

% Position (ECEF) Error
figure(5)
posECEF_imu = pos_ECEF;
for idx = 30:50
    posECEF_imu(idx,:) = pos_ECEF(idx-1,:) + distance(idx,:);
end
estPos = lla2ned(lla,lla_guess,"ellipsoid");
winSize = floor(size(estPos,1)/10);
plot(smoothdata(abs(posECEF_imu-receiver_ECEF),"movmedian",winSize))
legend("X - axis","Y - axis","Z - axis")
xlabel("Time [s]")
ylabel("Error [m]")
xlim([1,steps])
title("Position (ECEF) Error")

% Dilution of Precision (DOP)
figure(6)
plot(1:steps,hdop_diff);
hold on
plot(1:steps,vdop_diff);
plot(1:steps,hdop_alt,'--x');
plot(1:steps,vdop_alt,'--x');
hold off
legend('HDOP','VDOP','HDOP after altimeter', 'VDOP after altimeter')
ylabel('Distance deviation [m]')
xlim([1,steps])
xlabel("Time [s]")
title('Dilution of Precision (DOP)')

% Number of visible satellites
figure(7)
% find the number of visible satellites for each instant
visible_sat_num = zeros(steps,1);
for idx = 1:steps
    visible_sat_num(idx) = nnz(all_are_visible(:,idx));
end
area(1:steps, visible_sat_num, 'FaceColor','b','FaceAlpha',0.2);
title("Number of visible satellites");
xlabel('Time [s]');
xlim([1,steps])
ylim([6,12]);

% Satellite visibility
figure(8)
time = 1:steps;
satellites = 1:num_sats;
imagesc(time, satellites, all_are_visible);
clim([0 0.2]);
colormap([1 1 1; 0 0 1]);
xlabel('Time [s]');
ylabel('Satellite ID');
title('Satellite Visibility');
yticks(satellites);
yticklabels(satellites);
grid on

%% Functions

% get dilution of precision (DOP) values in local frame from
% Earth-Center-Earth-Fixed (ECEF) cofactor matrix
function [hdop, vdop,tdop] = calculate_DOP(dop_matrix, refFrame, lla_0)
    T = eye(4, 'like', dop_matrix);
    T(1:3,1:3) = refFrame.ecef2framerotmat(lla_0(1), lla_0(2));
    
    % Convert from ECEF to local frame.
    dop_matrix = T * dop_matrix * T.';
    
    diagDOPMat = diag(dop_matrix);
    
    hdop = sqrt(sum(diagDOPMat(1:2)));
    vdop = sqrt(diagDOPMat(3));
    tdop = sqrt(diagDOPMat(4));
end
