%% Readme

% this script simulates a gps receiver on the Moon. We use the file
% "data_Orbit_group.csv" we received from the orbital gruop because
% we were not provided with the navigation message, which is
% necessary to compute the orbits. Furthermore, we have discovered
% that the satellite consellation is designed in such a way that we almost
% never have four or more visible satellites. For this reason we are not
% able to provide a consistent position fix for each istant.

% we are not able to provide the same output plots as for GNSS_simulation.m
% beacuse we would obtain meaningless results since for the most part of
% the time we not able to see at least four satellites.

% the code structure is very similar to GNSS_simulation.m but it is adapted
% for Moon enviroment and different orbital parameters.

% As output we peovide:
%   - figure 1: Position of the receiver
%   - figure 2: Satellite orbits
%   - figure 3: Accelerometer Readings along X, Y, Z axis
%   - figure 4: Number of visible satellites
%   - figure 5: Satellite visibility


%% Reading the initial data

close all;
clearvars;
clc;

% simulated user trajectory (from orbitals group)
% format: latitude, longitude, position in xyz, velocity along x, y, z
initial_data = readmatrix("data_Orbit_group.csv");

% we need access to the user position, which defeats the whole purpose of
% GPS, because the receiver team did not provide us with anything
% meaningful: we do not know which satellites are visible and we do not
% have access to any navigation message
% receiver position (X, Y, Z)
rec_pos = initial_data(:,23:25);
% time interval as decided by orbital team
dt = 300;
% calculating the velocity of the receiver
rec_vel_x = calculateVel(rec_pos(:, 1), dt);
rec_vel_y = calculateVel(rec_pos(:, 2), dt);
rec_vel_z = calculateVel(rec_pos(:, 3), dt);
rec_vel = [rec_vel_x, rec_vel_y, rec_vel_z];

% vectors with information about each satellite
% format: sat ID, simulation time, X, Y, Z, distance, pseudorange.
sat1 = [ones(height(initial_data),1)*1, initial_data(:,1),initial_data(:,3), initial_data(:,7),initial_data(:,11), initial_data(:,15),  initial_data(:,19)];
sat2 = [ones(height(initial_data),1)*2, initial_data(:,1),initial_data(:,4), initial_data(:,8),initial_data(:,12), initial_data(:,16),  initial_data(:,20)];
sat3 = [ones(height(initial_data),1)*3, initial_data(:,1),initial_data(:,5), initial_data(:,9),initial_data(:,13), initial_data(:,17), initial_data(:,21)];
sat4 = [ones(height(initial_data),1)*4, initial_data(:,1),initial_data(:,6),initial_data(:,10),initial_data(:,14), initial_data(:,18) , initial_data(:,22)];

% specify simulation times as indicated on the orbit csv
start_time = datetime(2022,01,06,8,0,0,"TimeZone","America/New_York");
steps = size(rec_pos,1);

stop_time = start_time + seconds((steps-1)*dt);
duration = seconds(stop_time - start_time);
% mask angle [degrees], specified by orbit team - smaller than GNSS since
% on the Moon we assume to have less obstacles than down on Earth
maskAngle = 3;

% receiver position conversion from north-east-down coordinates to
% geodetic coordinates
receiver_ECEF = rec_pos;
receiver_LLA = ecef2lla(receiver_ECEF);


%% Calculate Recevier parameters

% calculate the distance between the receiver positions in two
% successive time instants
dist_receiver_x = zeros(steps-1,1);
dist_receiver_y=  zeros(steps-1,1);
dist_receiver_z = zeros(steps-1,1);
for ii = 1:steps-1
    dist_receiver_x(ii)= receiver_ECEF(ii+1,1) - receiver_ECEF(ii,1);
    dist_receiver_y(ii) = receiver_ECEF(ii+1,2) - receiver_ECEF(ii,2);
    dist_receiver_z(ii) = receiver_ECEF(ii+1,3) - receiver_ECEF(ii,3);
end
distance = [dist_receiver_x, dist_receiver_y, dist_receiver_z];

% time vector 1 to steps
time = 1:steps;

% calculate the accelerations along x, y, z axis
acceleration_x = rec_vel_x/dt;
acceleration_y = rec_vel_y/dt;
acceleration_z = rec_vel_z/dt;
acceleration = [acceleration_x,acceleration_y,acceleration_z];

% number of satellites in the constellation (decided by orbital team)
num_sats = 4;

% initialize satellites positions and velocities
sat_positions = zeros(num_sats,3,steps);
sat_velocities = zeros(num_sats,3,steps);
for idx = 1:steps
    sat_1_xyz = [initial_data(idx,3), initial_data(idx,7), initial_data(idx,11)];
    sat_2_xyz = [initial_data(idx,4), initial_data(idx,8), initial_data(idx,12)];
    sat_3_xyz = [initial_data(idx,5), initial_data(idx,9), initial_data(idx,13)];
    sat_4_xyz = [initial_data(idx,6), initial_data(idx,10), initial_data(idx,14)];
    sat_positions(:, :, idx) = [sat_1_xyz; sat_2_xyz; sat_3_xyz; sat_4_xyz];
end

% calculating the velocities along x, y, z for each one of the four
% satellites. Since we do not have the navigation message, we need to do
% this by hand by looking at their positions in time
for idx = 1:steps-1
    sat1_difference = [initial_data(idx+1,3) - initial_data(idx,3), initial_data(idx+1,7) - initial_data(idx,7), initial_data(idx+1,11) - initial_data(idx,11)];
    sat2_difference = [initial_data(idx+1,4) - initial_data(idx,4), initial_data(idx+1,8) - initial_data(idx,8), initial_data(idx+1,12) - initial_data(idx,12)];
    sat3_difference = [initial_data(idx+1,5) - initial_data(idx,5), initial_data(idx+1,9) - initial_data(idx,9), initial_data(idx+1,13) - initial_data(idx,13)];
    sat4_difference = [initial_data(idx+1,6) - initial_data(idx,6), initial_data(idx+1,10) - initial_data(idx,10), initial_data(idx+1,14) - initial_data(idx,14)];
    sat1_vel = sat1_difference/dt;
    sat2_vel = sat2_difference/dt;
    sat3_vel = sat3_difference/dt;
    sat4_vel = sat4_difference/dt;
    sat_velocities(:, :, idx) = [sat1_vel; sat2_vel; sat3_vel; sat4_vel];
end

% initialize vectors with all pseudoranges and pseudorange rates
all_pseudoranges = zeros(num_sats,steps);
all_pseudorange_rates = zeros(num_sats,steps);
% pseudoranges for differential LNSS (station point of view)
all_pseudoranges_diff = zeros(num_sats,steps);
% initialize vector with visibility information
all_are_visible = false(num_sats,steps);


%% Simulating the ground station - DLNSS

pseudo_sat_dgnss = zeros(num_sats, steps);
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
    [satAz, satEl, all_are_visible(:,idx)] = lookangles(receiver_LLA(idx,:),single_sat_pos,maskAngle);

    % returns the pseudorange and pseudorange rates, pdot, in meters per second,
    % between the receiver and satellites.
    [all_pseudoranges(:,idx), all_pseudorange_rates(:,idx)] = pseudoranges(receiver_LLA(idx,:),single_sat_pos,rec_vel(idx,:),single_sat_vel);
    % computing the pseudoranges but from station point of view
    all_pseudoranges_diff(:,idx) = pseudoranges(station_position,single_sat_pos);
    % manually computing the distance for each satellite
    for ii =1:num_sats
        pseudo_sat_dgnss(ii,idx) = sqrt((station_position(1)-single_sat_pos(ii,1))^2+(station_position(2)-single_sat_pos(ii,2))^2+(station_position(3)-single_sat_pos(ii,3))^2);
    end
end


%% DOP and position

% initializing latitude, longitude, and altitude coordinates
lla = zeros(steps,3);
% latitude, longitude, and altitude coordinates with altimeter and station
lla_alt = zeros(steps,3);
lla_diff = zeros(steps,3);
% initializing DOP parameters
hdop = zeros(steps,1);
vdop = zeros(steps,1);
% initializing DOP parameters for differential
hdop_diff = zeros(steps,1);
vdop_diff = zeros(steps,1);
% initializing DOP parameters for altimeter
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

for idx = 1:steps
    % copy pseudorange at step idx
    pseudoranges_idx = all_pseudoranges(:,idx);
    % take only the current values at instant idx
    pdot = all_pseudorange_rates(:,idx);
    is_visible = all_are_visible(:,idx);
    single_sat_pos = sat_positions(:,:,idx);
    single_sat_vel = sat_velocities(:,:,idx);
    % we had to put a condition on the number of visible satellites,
    % because otherwise the code throws many warnings of "matrix is singular"
    % since, because of the way the satellite constellation is designed, we almost
    % never have 4 or more visible satellites, which is the minimum number to be
    % able to perform the position fix. For that reason we take in consideration
    % data only when we can see at least four satellites.
    if sum(is_visible(:)) >= 4
        % position estimates
        pos_est = pos_prev;
        pos_est_diff = pos_prev_diff;

        % H matrices
        H_pos = ones(size(single_sat_pos(is_visible,:), 1), 4);
        H_pos_diff = ones(size(single_sat_pos(is_visible,:), 1), 4);
        H_pos_prev = H_pos;

        % we begin the computation of the position: initially, the resolution
        % is infinite (max uncertainty) and we want to achieve the resolution
        % specified in the variable target_res with the successive iterations
        res_pos = Inf;
        res_pos_diff = Inf;
        target_res = 1e-4;
        max_iter= 200;
        iter = 0;
        all_res_Pos = NaN(max_iter, 1);
        all_res_pos_diff = NaN(max_iter, 1);

        % loop until either we pass max iterations count or arrive to the
        % desired resolution
        while (res_pos > target_res) && (iter < max_iter)
            % p_est: pseudorange for all visible satellites 
            % los_vector: line-of-sight unit vectors from receiver to satellite
            [pEst, ~, los_vector] = nav.internal.gnss.calculatePseudoranges(single_sat_pos(is_visible,:), single_sat_vel(is_visible,:), pos_prev(1:3).', vel_prev(1:3).');
            % computing the prediction 
            p_pred = pEst + pos_prev(end);
            H_pos(:,1:3) = - los_vector;
            % calculate the next position estimation (state update equation)
            pos_est = pos_prev + H_pos.' * H_pos \ H_pos.' * (pseudoranges_idx(is_visible) - p_pred);
            res_pos = norm(pos_est - pos_prev);
            % increment iteration count
            iter = iter + 1;
            % save the current result
            all_res_Pos(iter) = res_pos;
            % set the current results as previous so we can proceed to the
            % next iteration
            pos_prev = pos_est;
            H_pos_prev = H_pos;
        end

        % second loop for altimeter computation
        % loop until either we pass max iterations count or arrive to the
        % desired resolution
        iter = 0;
        while (res_pos_diff > target_res) && (iter < max_iter)
            % dgnss position accounting for the ground station
            [pEst_dlnss, ~,~] = nav.internal.gnss.calculatePseudoranges(single_sat_pos(is_visible,:), single_sat_vel(is_visible,:), station_position, [0,0,0]);
            delta_rho = pEst_dlnss-pseudo_sat_dgnss(is_visible,idx);
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
        % normal measurements and altimeter
        pos_ECEF(idx,:) = pos_est(1:3).';
        pos_ECEF_diff(idx,:) = pos_est_diff(1:3).';

        dop_matrix = inv(H_pos.' * H_pos);
        dop_matrix_diff = inv(H_pos_diff.' * H_pos_diff);

        % convert the results in lat long alt for normal, differential and altimeter
        lla(idx,:) = ecef2lla(pos_ECEF(idx,:));
        lla_diff(idx,:) = ecef2lla(pos_ECEF_diff(idx,:));
        lla_alt(idx,:)= [lla_diff(idx,1:2),receiver_LLA(idx,3)+(randn(1,1)*0.2-0.1)];

        [hdop(idx), vdop(idx)] = calculate_DOP(dop_matrix, refFrame, lla);
        [hdop_diff(idx), vdop_diff(idx)] = calculate_DOP(dop_matrix_diff, refFrame, lla_diff);
        [hdop_alt(idx), vdop_alt(idx)] = calculate_DOP(dop_matrix_diff, refFrame, lla_alt);
    end
end


%% Results plotting

% note: there is an outlier, which is the initial guess of the user
% position
% 3D plot of receiver position
figure(1)
plot_pos_ECEF = pos_ECEF;
plot_pos_ECEF(plot_pos_ECEF==0) = nan;
plot3(plot_pos_ECEF(:,1)/1000,plot_pos_ECEF(:,2)/1000,plot_pos_ECEF(:,3)/1000, 'linestyle','none','marker','o');
legend("user trajectory")
xlabel("X - axis position [km]")
ylabel("Y - axis position [km]")
zlabel("Z - axis position [km]")
title("Position of the receiver")
grid on

% plotting the satellite orbits around the Moon
figure(2)
% moon sphere
[X,Y,Z] = sphere;
r = 1737;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;
plot3(sat1(:,3)/1000, sat1(:,4)/1000, sat1(:,5)/1000 , "g");
hold on
grid on
plot3(sat2(:,3)/1000, sat2(:,4)/1000, sat2(:,5)/1000 , "c");
plot3(sat3(:,3)/1000, sat3(:,4)/1000, sat3(:,5)/1000, "b" );
plot3(sat4(:,3)/1000, sat4(:,4)/1000, sat4(:,5)/1000, "y" );
moon = surf(X2,Y2,Z2, 'FaceColor', [0.66274,0.66274,0.66274]);
xlim([-8000, 8000]);
ylim([-8000, 8000]);
zlim([-8000, 8000]);
legend("sat 1", "sat 2", "sat 3", "sat 4");
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
title("Satellite orbits")

% plots with the accelerometer readings
figure(3)
% x axis acceleration
subplot(3,1,1)
plot(time,acceleration_x+randn(steps, 1)*(-1).^(round(randn(1,1))),'Color',[0 0.4470 0.7410])
hold on
legend('X-axis')
xlabel("Time [s]")
ylabel('Acceleration [m/s^2]')
xlim([7,steps])
ylim([-8, 8]);
title('Accelerometer Readings')

% y axis acceleration
subplot(3,1,2)
plot(time,acceleration_y+randn(steps, 1)*(-1).^(round(randn(1,1))),'Color',[0.8500 0.3250 0.0980])
hold on
legend('Y-axis')
xlim([5,steps])
ylim([-8, 8]);
xlabel("Time [s]")
ylabel('Acceleration [m/s^2]')

% z axis acceleration
subplot(3,1,3)
plot(time,acceleration_z+randn(steps, 1)*(-1).^(round(randn(1,1))),'Color',[0.8500 0.3250 0.0980])
hold off
legend('Z-axis')
xlim([5,steps])
ylim([-8, 8]);
xlabel("Time [s]")
ylabel('Acceleration [m/s^2]')

% plot the number of visible satellites
figure(4)
% find the number of visible satellites for each instant
visible_sat_num = zeros(1, 7869);
for idx = 1:steps
    visible_sat_num(idx) = nnz(all_are_visible(:,idx));
end
area(1:1:steps, visible_sat_num, 'FaceColor','b','FaceAlpha',0.2);
title("Number of visible satellites");
xlabel('Time [s]');
xlim([1,steps])
ylim([0,5]);
ay = gca;
ay.YTick = unique(round(ay.YTick));

% plot the visibility for each satellite
figure(5)
satellites = 1:num_sats;
imagesc(time, satellites, all_are_visible);
clim([0 0.2]);
colormap([1 1 1; 0 0 1]);
xlabel('Time [s]');
ylabel('Satellite ID');
title('Satellite Visibility');
yticks(satellites);
yticklabels(satellites);

%% Other plots

% in the file GNSS_simulation.m we were able to extract and plot more results,
% but here we can't provide them since we can't get a stable position fix
% because we almost never have 4 satellites visible at the same time.
% Please refer to the file for the GNSS to have the working plots with
% data that actually make sense.


%% Functions

% calculate the instantaneous velocity starting from a vector of positions
% and time interval between positions
function out = calculateVel(positions, dt)
    diffPositions = diff(positions);
    out = [0; diffPositions/dt];
    out(2) = 0;
end

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