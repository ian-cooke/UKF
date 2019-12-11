function [xhat, Phat, resid, mat_cond] = ukf(tspan, xhat_0, Phat_0, measurements, config)
%UKF Unscented Kalman Filter

% extract params
% Extract params
data = config.data;
n = config.n;
alpha = config.alpha;
beta = config.beta;
kappa = config.kappa;
R_meas = config.R_meas;
options = config.options;
stations_state_eci = config.stations_state_eci;
use_J3_filter = config.use_J3_filter;

% process noise stuff
use_snc = config.use_snc;
sigma_u = config.sigma_u;
Gamma = config.Gamma;
Q_CT = sigma_u^2*eye(3);

% Loop setup
N = length(tspan);
xhat = zeros(n, N);
Phat = zeros(n, n, N);
mat_cond = zeros(1, N);
resid = NaN(6, N);
% Populate first index
xhat(:, 1) = xhat_0;
Phat(:, :, 1) = Phat_0;
mat_cond(1) = cond(Phat_0);

% ukf specific parameter computations
lambda = alpha^2*(kappa + n) - n;
gamma = sqrt(n + lambda);
% compute weights
w_m = zeros(1, 2*n+1);
w_c = w_m;
for j = 1:2*n+1
    if j == 1
        w_m(j) = lambda/(lambda + n);
        w_c(j) = w_m(j) + (1 - alpha^2 + beta);
    else
        w_m(j) = 0.5/(lambda+n);
        w_c(j) = w_m(j);
    end
end

% Initialize loop variables
x_kP = xhat_0;
P_kP = Phat_0;

for k = 1:N-1
    
    % 1) compute sigma points for k
    % compute principal square root of P_kP (getting sigma_j's)
    S_k = sqrtm(P_kP);
    capX_kP = zeros(n, 2*n+1);
    capX_kp1M = zeros(n, 2*n+1);
    for j = 1:2*n+1
        if j == 1
            capX_kP(:, j) = x_kP;
        elseif (j > 1) && (j < n + 2)
            capX_kP(:, j) = x_kP + gamma*S_k(:, j-1);
        else
            capX_kP(:, j) = x_kP - gamma*S_k(:, j-(n+1));
        end
        % propagate forward with full nonlinear dynamics
        if use_J3_filter
            [~, x_temp] = ode113(@(t,x) keplerJ2J3_wPhi_ODE(t, x, config), [tspan(k), tspan(k+1)],...
                [capX_kP(:, j); reshape(eye(n), n^2, 1)], options);
        else
            [~, x_temp] = ode113(@(t,x) keplerJ2_wPhi_ODE(t, x, config), [tspan(k), tspan(k+1)],...
                [capX_kP(:, j); reshape(eye(n), n^2, 1)], options);
        end
        capX_kp1M(:, j) = x_temp(end, 1:n)';
    end
    
    % 2) calculate x_kp1M and P_kp1M for kp1
    sum_x = 0;
    sum_P = 0;
    % compute x_kp1M
    for j = 1:2*n+1
        sum_x = sum_x + w_m(j)*capX_kp1M(:, j);
    end
    x_kp1M = sum_x;
    % compute P_kp1M
    for j = 1:2*n+1
        sum_P = sum_P + w_c(j)*(capX_kp1M(:, j) - x_kp1M)*(capX_kp1M(:, j) - x_kp1M)';
    end
    P_kp1M = sum_P;
    
    % 3) if using process noise, need to recompute sigma points
    if use_snc
        P_kp1M = P_kp1M + Gamma*Q_CT*Gamma';
        % recompute sigma points based on time update mean and covariance
        S_k = sqrtm(P_kp1M);
        % compute sigma points, kp1M
        capX_kp1M = zeros(n, 2*n+1);
        for j = 1:2*n+1
            if j == 1
                capX_kp1M(:, j) = x_kp1M;
            elseif (j > 1) && (j < n + 2)
                capX_kp1M(:, j) = x_kp1M + gamma*S_k(:, j-1);
            else
                capX_kp1M(:, j) = x_kp1M - gamma*S_k(:, j-(n+1));
            end
        end
    end
    
    % 4) measurement update equations
    current_meas = measurements(:, k+1);
    station = get_station(current_meas);
    % only do measurement update if measurement is available
    if station ~= 0
        current_stations_state_eci = [stations_state_eci(:, k+1, 1), stations_state_eci(:, k+1, 2), stations_state_eci(:, k+1, 3)]; 
        [x_obs, y_true] = get_x_obs_meas(station, current_stations_state_eci, current_meas, data);
        
        % compute capY for each sigma point, and yhat_kM
        y_predicted = 0;
        capY = zeros(2, 2*n+1);
        for j = 1:2*n+1
            capY(:, j) = get_rho_rhod(capX_kp1M(:, j), x_obs, data);
            y_predicted = y_predicted + w_m(j)*capY(:, j);
        end
        
        % compute P_xy and P_yy
        P_xy = 0;
        P_yy = 0;
        for j = 1:2*n+1
            P_yy = P_yy + w_c(j)*(capY(:, j) - y_predicted)*(capY(:, j) - y_predicted)';
            P_xy = P_xy + w_c(j)*(capX_kp1M(:, j) - x_kp1M)*(capY(:, j) - y_predicted)';
        end
        
        % compute the kalman gain
        K_k = P_xy*inv(P_yy + R_meas);
        
        % Kalman update equations - the same (almost)
        x_kp1P = x_kp1M + K_k*(y_true - y_predicted);
        P_kp1P = P_kp1M - K_k*P_yy*K_k';
        
        % 5) Calculate the post-fit measurement residual
        x_sc = x_kp1P;
        current_stations_state_eci = [stations_state_eci(:, k+1, 1), stations_state_eci(:, k+1, 2), stations_state_eci(:, k+1, 3)]; 
        [x_obs, y_true] = get_x_obs_meas(station, current_stations_state_eci, current_meas, data);
        y_predicted = get_rho_rhod(x_sc, x_obs, data);
        resid(2*station-1:2*station, k+1) = y_true - y_predicted;
        
    % if no measurement available
    else
        x_kp1P = x_kp1M;
        P_kp1P = P_kp1M;
    end
    
    % get post-update matrix condition number
    mat_cond(k+1) = cond(P_kp1P);
    
    % store
    xhat(:, k+1) = x_kp1P;
    Phat(:, :, k+1) = P_kp1P;
    % update
    x_kP = x_kp1P;
    P_kP = P_kp1P;
       
end

end