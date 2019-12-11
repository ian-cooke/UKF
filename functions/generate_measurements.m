function [measurements, stations_state_eci] = generate_measurements(config)
% GENERATE_MEASUREMENTS

% params
lat_lon = config.lat_lon;
R_meas = config.R_meas;
dt = config.dt;
mask = config.mask;
theta_0 = config.theta_0;
w_E = config.w_E;
R_E = config.R_E;
rot_E = config.rot_E;
x_true = config.x_true;
n = config.n;

% Calculate station ecef coordinates
n_stations = size(lat_lon, 1);
stations_ecef = latlon2ecef(n_stations, lat_lon, R_E);

% Loop over x_base to create measurements and add noise
N = size(x_true,2);
measurements = NaN(2*n_stations, N);
stations_state_eci = zeros(n, N, n_stations);
theta = theta_0;
for k = 1:N
    % spacecraft pos and velocity
    R = x_true(1:3, k);
    V = x_true(4:6, k);
    % calculate dcm to translate ecef to eci
    Q_ECEF_ECI = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
    % loop over each station
    for j = 1:n_stations
        % Station position and velocity
        R_s = Q_ECEF_ECI*stations_ecef(:, j);
        V_s = cross(rot_E, R_s);
        stations_state_eci(:, k, j) = [R_s;V_s];
        % range and range rate
        rho = norm(R - R_s);
        rhod = (R - R_s)'*(V - V_s)/rho;
        % line of sight and elevation angle
        LOS = (R - R_s)/rho;
        elev = rad2deg(acos(dot(R_s, LOS)/(norm(R_s))));
        % check if elevation angle above mask angle, if so add noise and 
        % place in array
        if elev < 90 - mask
            measurements(2*j-1:2*j,k) = [rho;rhod] + mvnrnd([0,0], R_meas)';
        end
    end
    % update angle
    theta = theta + w_E*dt;
end
        
    
end