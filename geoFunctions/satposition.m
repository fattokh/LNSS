function [satPositions, satClkCorr] = satpos(transmitTime, prnList, eph, settings) 
numOfSatellites = height(prnList);

% GPS constatns


gpsPi          = 3.1415926535898;  % Pi used in the GPS coordinate 
                                   % system

%--- Constants for satellite position calculation -------------------------
Omegae_dot     = 7.2921151467e-5;  % Earth rotation rate, [rad/s]
GM             = 3.986005e14;      % Universal gravitational constant times
                                   % the mass of the Earth, [m^3/s^2]
F              = -4.442807633e-10; % Constant, [sec/(meter)^(1/2)]

%% Initialize results =====================================================
satClkCorr   = zeros(1, numOfSatellites);
satPositions = zeros(3, numOfSatellites);

%% Process each satellite =================================================

for satNr = 1 : numOfSatellites
    
    prn = prnList(satNr);
    dtr = zeros(satNr,1);
%% Find initial satellite clock correction --------------------------------

    %--- Find time difference ---------------------------------------------
    dt(satNr) = check_t(transmitTime(satNr) - eph.t_oc(satNr));

    %--- Calculate clock correction ---------------------------------------
    satClkCorr(satNr) = (eph.a_f2(satNr) * dt(satNr) + eph.a_f1(satNr)) * dt(satNr) + ...
                         eph.a_f0(satNr) - ...
                         eph.T_GD(satNr);

    time(satNr) = transmitTime(satNr) - satClkCorr(satNr);

%% Find satellite's position ----------------------------------------------

    %Restore semi-major axis
    a(satNr)   = eph.sqrtA(satNr) * eph.sqrtA(satNr);

    %Time correction
    tk(satNr)  = check_t(time(satNr) - eph.t_oe(satNr));

    %Initial mean motion
    n0(satNr)  = sqrt(GM / a(satNr)^3);
    %Mean motion
    n(satNr)   = n0(satNr) + eph.deltan(satNr);

    %Mean anomaly
    M(satNr)   = eph.M_0(satNr) + n(satNr) * tk(satNr);
    %Reduce mean anomaly to between 0 and 360 deg
    M(satNr)   = rem(M(satNr) + 2*gpsPi, 2*gpsPi);

    %Initial guess of eccentric anomaly
    E(satNr)   = M(satNr);

    %--- Iteratively compute eccentric anomaly ----------------------------
    for ii = 1:10
        E_old   = E;
        E       = M + eph.e(satNr) * sin(E);
        dE      = rem(E - E_old, 2*gpsPi);

        if abs(dE) < 1.e-12
            % Necessary precision is reached, exit from the loop
            break;
        end
    end

    %Reduce eccentric anomaly to between 0 and 360 deg
    E(satNr)   = rem(E(satNr) + 2*gpsPi, 2*gpsPi);

    %Compute relativistic correction term
    dtr(satNr) = F * eph.e(satNr) * eph.sqrtA(satNr) * sin(E(satNr));

    %Calculate the true anomaly
    nu(satNr)   = atan2(sqrt(1 - eph.e(satNr)^2) * sin(E(satNr)), cos(E(satNr))-eph.e(satNr));

    %Compute angle phi
    phi(satNr) = nu(satNr) + eph.omega(satNr);
    %Reduce phi to between 0 and 360 deg
    phi(satNr) = rem(phi(satNr), 2*gpsPi);

    %Correct argument of latitude
    u(satNr) = phi(satNr) + eph.C_uc(satNr) * cos(2*phi(satNr)) + ...
        eph.C_us(satNr) * sin(2*phi(satNr));
    %Correct radius
    r(satNr) = a(satNr) * (1 - eph.e(satNr)*cos(E(satNr))) + ...
        eph.C_rc(satNr) * cos(2*phi(satNr)) + ...
        eph.C_rs(satNr) * sin(2*phi(satNr));
    %Correct inclination
    i(satNr) = eph.i_0(satNr) + eph.iDot(satNr) * tk(satNr) + ...
        eph.C_ic(satNr) * cos(2*phi(satNr)) + ...
        eph.C_is(satNr) * sin(2*phi(satNr));

    %Compute the angle between the ascending node and the Greenwich meridian
    Omega(satNr) = eph.omega_0(satNr) + (eph.omegaDot(satNr) - Omegae_dot)*tk(satNr) - ...
            Omegae_dot * eph.t_oe(satNr);
    %Reduce to between 0 and 360 deg
    Omega(satNr) = rem(Omega(satNr) + 2*gpsPi, 2*gpsPi);

    %--- Compute satellite coordinates ------------------------------------
    satPositions(1, satNr) = cos(u(satNr))*r(satNr) * cos(Omega(satNr)) - sin(u(satNr))*r(satNr) * cos(i(satNr))*sin(Omega(satNr));
    satPositions(2, satNr) = cos(u(satNr))*r(satNr) * sin(Omega(satNr)) + sin(u(satNr))*r(satNr) * cos(i(satNr))*cos(Omega(satNr));
    satPositions(3, satNr) = sin(u(satNr))*r(satNr) * sin(i(satNr));


%% Include relativistic correction in clock correction --------------------
    satClkCorr(satNr) = (eph.a_f2(satNr) * dt(satNr) + eph.a_f1(satNr)) * dt(satNr) + ...
                         eph.a_f0(satNr) - ...
                         eph.T_GD(satNr) + dtr(satNr);   
end



 % for satNr = 1 : numOfSatellites
end