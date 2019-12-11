function [xhat, Phat, resid, mat_cond] = ekf(tspan, xhat_0, Phat_0, measurements, config)
%CKF classic Kalman filter, using Joseph formulation for update
%
% Inputs
%   tspan - vector of times
%   dxhat_0 - initial estimated state deviation
%   xhat_0 - initial estimated state
%   Phat_0 - initial estimated state error covariance
%   measurements - range and range rate measurements
%   config - every other thing the ckf needs to run, including
%       data - string that tells which data to use 'rho', 'rhod', or 'both'
%       m - number of variables in state vector
%       w_E - [rad/s] rotation rate of the Earth
%       R_E - [m] mean radius of Earth
%       rot_E - [rad/s] vector of angular velocity difference between eci
%           and ecef frames
%       R_meas - [m^2, m^2/s^2] measurement errror covariance
%       A_sat - [m^2] drag area of satellite
%       m_sat - [kg] mass of satellite
%       rho_0 - [kg/m^3] initial density
%       H - [m] constant in exponential drag model
%       r_0 - [m] constant in exponential drag model
%       options - ode tolderances
% Outputs
%   dxhat - esimtated state deviation
%   Phat - esimated state error covariance
%   traj_ref - reference trajecrory (includes STM)
%   resid - post-fit measurement residuals
%   mat_cond - condition number of posterior state error covariance matrix

% Extract params
data = config.data;
n = config.n;
R_meas = config.R_meas;
options = config.options;
sigma_u = config.sigma_u;
use_snc = config.use_snc;
use_ric = config.use_ric;
stations_state_eci = config.stations_state_eci;
Q_CT = sigma_u^2*eye(3);
Gamma = config.Gamma;

% Loop setup
N = length(tspan);
xhat = zeros(n, N);
Phat = zeros(n, n, N);
mat_cond = zeros(1, N);
mat_cond(1) = cond(Phat_0);
% if you only want one type of data
switch data
    case 'both'
        resid = NaN(6, N);
    otherwise
        resid = NaN(3, N);
        if strcmp(data, 'rho')
            R_meas = R_meas(1,1);
        else
            R_meas = R_meas(end,end);
        end
end
% Populate first index
xhat(:, 1) = xhat_0;
Phat(:, :, 1) = Phat_0;

% Get first post-fit residual
current_meas = measurements(:, 1);
station = get_station(current_meas);
current_stations_state_eci = [stations_state_eci(:, 1, 1), stations_state_eci(:, 1, 2), stations_state_eci(:, 1, 3)]; 
[x_obs, y_true] = get_x_obs_meas(station, current_stations_state_eci, current_meas, data);
x_sc = xhat(1:6, 1);
y_predicted = get_rho_rhod(x_sc, x_obs, data);
res = y_true - y_predicted;
% in case you want only one type of data
switch data
    case 'both'
        resid(2*station-1:2*station, 1) = res;
    otherwise
        resid(station, 1) = res;
end

% Initialize loop variables
x_kP = xhat_0;
P_kP = Phat_0;

for k = 1:N-1
    % 1) Propagate reference trajectory, update, and get STM
    [~, x_temp] = ode113(@(t,x) keplerJ2_wPhi_ODE(t, x, config), [tspan(k), tspan(k+1)],...
        [x_kP(1:n); reshape(eye(n), n^2, 1)], options);
    x_temp = x_temp(end, :)';
    x_kp1M = x_temp(1:n);
    Phi_k = reshape(x_temp(n+1:end), n, n);
    
    % 2) Prediction step
    % SNC or not SNC
    if use_snc
        % w/ SNC
        % if using ric frame
        if use_ric
            Q_ECI_RIC = eci2ric_dcm(traj_ref(1:3, k) + dx_kP(1:3), traj_ref(4:6, k) + dx_kP(4:6));
            Q_CT = Q_ECI_RIC*Q_CT*Q_ECI_RIC';
        end
        P_kp1M = Phi_k*P_kP*Phi_k' + Gamma*Q_CT*Gamma';
    else
        % w/o SNC
        P_kp1M = Phi_k*P_kP*Phi_k';
    end
    
    % 3) Observation deviation, H_tilde, gain matrix - only calcualte if
    % there is actually a measurement
    % get station, current measurement & trajectory based on station,
    % predicted measurement based on trjectory, and observation deviation r
    % see which station is available & get measurement
    current_meas = measurements(:, k+1);
    station = get_station(current_meas);
    if station ~= 0
        x_sc = x_kp1M;
        current_stations_state_eci = [stations_state_eci(:, k+1, 1), stations_state_eci(:, k+1, 2), stations_state_eci(:, k+1, 3)]; 
        [x_obs, y_true] = get_x_obs_meas(station, current_stations_state_eci, current_meas, data);
        y_predicted = get_rho_rhod(x_sc, x_obs, data);
        r_k = y_true - y_predicted;

        % get H_tilde based on current spacecraft state and station states
        Htilde = get_Htilde(x_sc, x_obs, data);
        % ensure symmetric
        %P_kp1M = (P_kp1M + P_kp1M')/2;
        % Kalman Gain - inverse is / for b*inv(A)
        K_k = (P_kp1M*Htilde')*inv(Htilde*P_kp1M*Htilde' + R_meas);
        
        % 4) Measurement update - if there is a measurement
        x_kp1P = x_kp1M + K_k*r_k;
        P_kp1P = (eye(n) - K_k*Htilde)*P_kp1M*(eye(n) - K_k*Htilde)' + K_k*R_meas*K_k';
        % ensure symmetric
        %P_kp1P = (P_kp1P + P_kp1P')/2;
        
        % calculate the condition number
        mat_cond(k+1) = cond(P_kp1P);
        
        % 5) Calculate the post-fit measurement residual by getting new
        % y_predicted
        x_sc = x_kp1P;
        current_stations_state_eci = [stations_state_eci(:, k+1, 1), stations_state_eci(:, k+1, 2), stations_state_eci(:, k+1, 3)]; 
        [x_obs, y_true] = get_x_obs_meas(station, current_stations_state_eci, current_meas, data);
        y_predicted = get_rho_rhod(x_sc, x_obs, data);
        res = y_true - y_predicted;
        % in case you want only one type of data
        switch data
            case 'both'
                resid(2*station-1:2*station, k+1) = res;
            otherwise
                resid(station, k+1) = res;
        end
        
        % update for next iteration
        x_kP = x_kp1P;
        P_kP = P_kp1P;

    % If no measurement, just pure prediction
    else
        % update for next iteration
        x_kP = x_kp1M;
        P_kP = P_kp1M;
        % calculate the condition number
        mat_cond(k+1) = cond(P_kP);
    end
    
    
    % store
    xhat(:, k+1) = x_kP;
    Phat(:, :, k+1) = P_kP;
    
    
end
    


