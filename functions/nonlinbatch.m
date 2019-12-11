function [dxhat_0_batch, Phat_0_batch, Phat, xhat, dxhat, resid] = nonlinbatch(tspan, dxhat_0, xhat_0, Phat_0, measurements, config)
%NONLINBATCH Nonlinear Batch Filter, iterated three times
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
%   dxhat_0_batch - vector of batch estimated intitial state deviations
%   Phat_0_batch - 3d matrix of batch estimated initial covariances
%   Phat - 3d matrix of estimated error covariances
%   xhat_batch - reference trajectory computed after final batch iteration
%   resid - post-fit measurement residuals from reference trajectory

% Extract params
data = config.data;
n = config.n;
w_E = config.w_E;
R_meas = config.R_meas;
inv_R_meas = inv(R_meas);
options = config.options;
n_iter = config.n_iter;
stations_state_eci = config.stations_state_eci;

% Batch loop setup
dxhat_0_batch = zeros(n, n_iter + 1);
dxhat_0_batch(:, 1) = dxhat_0;
Phat_0_batch = zeros(n, n, n_iter + 1);
Phat_0_batch(:, :, 1) = Phat_0;
N = length(tspan);
traj_ref = zeros(n^2+n, N);
xstar_0 = xhat_0;
dx_0_prior = dxhat_0;

% main batch loop
iter = 1;
while iter <= n_iter
    % Setup for inner loop a priori estimates
    Lambda = inv(Phat_0);
    N_mat = Lambda*dx_0_prior;
    traj_ref(:, 1) = xstar_0;
    
    % Inner loop
    for k = 1:N-1
        
        % 1) Integrate reference trajectory & get Phi
        [~, x_temp] = ode113(@(t,x) keplerJ2_wPhi_ODE(t, x, config), [tspan(k), tspan(k+1)],...
            traj_ref(:, k), options);
        traj_ref(:, k+1) = x_temp(end, :)';
        Phi_k = reshape(traj_ref(n+1:end, k+1), n, n);
        
        % 2) Get current measurement and station
        current_meas = measurements(:, k+1);
        station = get_station(current_meas);
        if station ~= 0 % if no observation, Lambda and N remain unchanged
            % 3) get current predicted measurement
            current_stations_state_eci = [stations_state_eci(:, k+1, 1), stations_state_eci(:, k+1, 2), stations_state_eci(:, k+1, 3)]; 
            [x_obs, y_true] = get_x_obs_meas(station, current_stations_state_eci, current_meas, data);
            x_sc = traj_ref(1:6, k+1);
            y_predicted = get_rho_rhod(x_sc, x_obs, data);
            % 4) compute residual
            r_k = y_true - y_predicted;
            
            % 5) get H_tilde based on current spacecraft state and station states
            Htilde = get_Htilde(x_sc, x_obs, data);
            % calculate modified H matrix
            H = Htilde*Phi_k;
            
            % 6) sum up
            switch data
                case 'both'
                    Lambda = Lambda + H'*inv_R_meas*H;
                    N_mat = N_mat + H'*inv_R_meas*r_k;
                otherwise
                    Lambda = Lambda + H'*1/R_meas(1,1)*H;
                    N_mat = N_mat + H'*1/R_meas(2,2)*r_k;
            end
        end
    end
      
    % solve normal equations & store
    P_0 = inv(Lambda);
    dx_0_post = P_0*N_mat;
    dxhat_0_batch(:, iter+1) = dx_0_post;
    Phat_0_batch(:, :, iter+1) = P_0;
    
    % update initial reference state
    xstar_0(1:n) = xstar_0(1:n) + dx_0_post;
    xstar_0(n+1:end) = reshape(eye(n), n^2, 1);
    dx_0_prior = dx_0_prior - dx_0_post;
    
        
    % iterate loop variable
    iter=iter+1;
end

% Trajectory Calculation Loop Setup
% if you only want one type of data
switch data
    case 'both'
        resid = NaN(6, N);
    otherwise
        resid = NaN(3, N);
end
% calculate new trajectory
[~, xhat] = ode113(@(t,x) keplerJ2_wPhi_ODE(t, x, config), tspan,...
            xstar_0, options);
xhat = xhat';

% loop through to calculate the covariance and residuals
Phat = zeros(n, n, N);
P_0 = Phat_0_batch(:, :, end);
for k = 1:N
    Phi_1_k = reshape(xhat(n+1:end, k), n, n);
    Phat(:, :, k) = Phi_1_k*P_0*Phi_1_k';
    % Calculate residual
    current_meas = measurements(:, k);
    station = get_station(current_meas);
    if station ~= 0 
        current_traj = xhat(1:n, k);
        current_stations_state_eci = [stations_state_eci(:, k, 1), stations_state_eci(:, k, 2), stations_state_eci(:, k, 3)];
        [x_obs, y_true] = get_x_obs_meas(station, current_stations_state_eci, current_meas, data);
        x_sc = current_traj(1:6);
        y_predicted = get_rho_rhod(x_sc, x_obs, data);
        
        % Calculate residual
        res = y_true - y_predicted;
        % in case you want only one type of data
        switch data
            case 'both'
                resid(2*station-1:2*station, k) = res;
            otherwise
                resid(station, k) = res;
        end
    end
end

% Generate reference trajectory to get dxhat
[~, x_ref] = ode113(@(t,x) keplerJ2_wPhi_ODE(t, x, config), tspan,...
            [xstar_0(1:n) + dx_0_post; reshape(eye(n), n^2, 1)], options);
x_ref = x_ref';
dxhat = xhat(1:n, :) - x_ref(1:n, :);
xhat = xhat(1:n, :);

end

