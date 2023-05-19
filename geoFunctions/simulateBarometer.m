function pressure = simulateBarometer()
    % Constants for the simulation
    seaLevelPressure = 1013.25; % Standard sea level pressure in hPa
    altitude = 0; % Initial altitude in meters
    lapseRate = 0.0065; % Lapse rate in K/m
    temperature = 288.15; % Initial temperature in K

    % Simulation parameters
    timeStep = 1; % Time step in seconds
    totalTime = 3600; % Total simulation time in seconds

    % Generate time vector
    time = 0:timeStep:totalTime;

    
    % Preallocate pressure vector
    pressure = zeros(size(time));

    % Simulate barometric pressure over time
    for i = 1:numel(time)
        % Calculate current pressure based on altitude and temperature
        pressure(i) = seaLevelPressure * (1 - (lapseRate * altitude) / temperature)^(5.255);

        % Update altitude using a dummy model (e.g., constant ascent rate)
        ascentRate = 1; % Assumed ascent rate in meters per second
        altitude = altitude + ascentRate * timeStep;

        % Update temperature using a dummy model (e.g., constant lapse rate)
        temperature = temperature - lapseRate * altitude * timeStep;
    end

    % Add random noise to the pressure signal (optional)
    noiseMagnitude = 0.5; % Adjust the magnitude of noise as desired
    pressure = pressure + noiseMagnitude * randn(size(pressure));
end
