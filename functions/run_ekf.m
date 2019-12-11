function dataset = run_ekf(config, measurements)

% Extract variables
tspan = config.tspan;
xhat_0 = config.xhat_0_ekf;
Phat_0 = config.Phat_0;

% Run EKF
[xhat, Phat, resid, mat_cond] = ekf(tspan, xhat_0, Phat_0, measurements, config);

% Get Dataset
dataset = get_dataset_ekf(xhat, Phat, resid, mat_cond, config);

end

