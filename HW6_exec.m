%HW6 script for execution
%% Setup
clear; close all; restoredefaultpath;
addpath('functions');
addpath('macros');
rng(100);

%% Inputs
% filter (ekf or ukf)
filter = 'ukf';
filter_2 = 'none';

% truth ('J2, J2drag', 'J2J3', etc.) - tells the program which integrator
% to use for truth values
truth = 'J2J3';

% use snc (0 or 1)
use_snc = 1;

% for doing part d, where you compare RMS values for various initial
% conditions
compare_ekf_ukf = 1;
n_comparisons = 15;
scale_interval = 20;

% use the ric frame
use_ric = 0;

% use J3 filter dynamics
use_J3_filter = 1;

% which types of data to use
data = 'both';

% number of vars in state
n = 6;

% orbit elems
a = 10000e3;
e = 0.001;
inc = 40;
Omega = 80;
omega = 40;
nu_0 = 0;

% Dynamics
J2 = 1.082626925638815e-03; % [none]
J3 = -2.33936e-3*J2;
mu = 3.986004415e14; % [km^3/s^2]
R_E = 6378136.3; % [km]
w_E = 7.2921158553*10^-5; % [rad/s] 
rot_E = [0;0;w_E];

% convert to r,v for initial spacecraft state
[R_0, V_0] = kep2eci(a,e,inc,Omega,omega,nu_0,mu);

% time step, orbital period, and time span
dt = 20;
T = 2*pi*a^(3/2)/sqrt(mu);
n_orbits = 15;
if compare_ekf_ukf
    tspan = 0:dt:5*T;
else
    tspan = 0:dt:n_orbits*T;
end

% station parameters
lat_lon = [-35.398333, 148.981944;...
           40.427222, 355.749444;...
           35.247164, 243.205];
mask = 10; % [deg] elevation mask
theta_0 = deg2rad(122);

% ode tolerances
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Initial State and Covariance
x_0 = [R_0; V_0];
sigma_r = 1e3; % [m]
sigma_v = 1; % [m/s]
Phat_0 = diag([ones(1,3)*sigma_r, ones(1,3)*sigma_v]).^2;

% Measurement covariance
sigma_rho = 1; % [m]
sigma_rhod = 1e-3; % [m/s]
R_meas = diag([sigma_rho, sigma_rhod]).^2;

% State noise sigma
sigma_u = 1e-5;
Gamma = [eye(3)*1/2*dt^2; eye(3)*dt];

% ukf parameters
alpha = 1;
beta = 2;
kappa = 0;

% Initial state deviation - sampled from initial covariance - will use
% the same deviation for every run depending on random seed
dxhat_0 = mvnrnd(zeros(1,n), Phat_0); dxhat_0 = dxhat_0';
%dxhat_0 = zeros(n,1);
%dxhat_0 = [-.5,-0.2,0.5,0.5,-0.5,0.2]';


% Assemble config structure
config.filter = filter; config.filter_2 = filter_2;
config.use_snc = use_snc;
config.use_ric = use_ric;
config.data = data;
config.use_J3_filter = use_J3_filter;
config.n = n;
config.a = a; config.e = e; config.inc = inc; config.Omega = Omega; config.omega = omega; config.nu_0 = nu_0;
config.J2 = J2; config.J3 = J3; config.mu = mu; config.R_E = R_E; config.w_E = w_E; config.rot_E = rot_E;
config.R_meas = R_meas;
config.dt = dt; config.tspan = tspan;
config.lat_lon = lat_lon; config.mask = mask; config.theta_0 = theta_0;
config.options = options;
config.x_0 =  x_0; config.Phat_0 = Phat_0;
config.xhat_0 = [x_0; reshape(eye(n), n^2, 1)];
config.xhat_0_ekf = x_0 + dxhat_0;
config.dxhat_0 = dxhat_0;
config.Gamma = Gamma;
config.sigma_u = sigma_u;
config.alpha = alpha; config.beta = beta; config.kappa = kappa;

% calculate 'base' truth data and append to config structure
switch truth
    case 'J2'
        [~, x_true] = ode113(@(t,x) keplerJ2_wPhi_ODE(t, x, config), tspan,...
            [x_0; reshape(eye(n), n^2, 1)], options);
    case 'J2J3'
        [~, x_true] = ode113(@(t,x) keplerJ2J3_wPhi_ODE(t, x, config), tspan,...
            [x_0; reshape(eye(n), n^2, 1)], options);
end
x_true=x_true';
config.x_true = x_true;

% Get data
[measurements, stations_state_eci] = generate_measurements(config);

% add station states to config
config.stations_state_eci = stations_state_eci;

%% Run Filters
% comparing multiple runs
if compare_ekf_ukf
    % preallocate and set new tspan
    RMS_pos_ekf = zeros(1, n_comparisons);
    RMS_vel_ekf = RMS_pos_ekf;
    RMS_pos_ukf = RMS_pos_ekf;
    RMS_vel_ukf = RMS_pos_ekf;
    % loop over number of comparisons
    for i = 1:n_comparisons
        config.xhat_0_ekf = x_0 + dxhat_0*(i/(1/scale_interval)+(1-scale_interval));
        dataset = run_ekf(config, measurements);
        RMS_pos_ekf(i) = dataset.RMS_3d_pos;
        RMS_vel_ekf(i) = dataset.RMS_3d_vel;
        dataset = [];
        dataset = run_ukf(config, measurements);
        RMS_pos_ukf(i) = dataset.RMS_3d_pos;
        RMS_vel_ukf(i) = dataset.RMS_3d_vel;
        dataset = [];
    end
    f3 = figure(3);
    plotter(f3, [1:15]./(1/scale_interval) + (1-scale_interval), 'compare_ekf_ukf', data, filter, filter_2, [RMS_pos_ekf; RMS_pos_ukf], [RMS_vel_ekf; RMS_vel_ukf]);
% non comparing, single runs
else
    switch filter
        case 'ekf'
            dataset = run_ekf(config, measurements);
        case 'ukf'
            dataset = run_ukf(config, measurements);
        otherwise
            error('Not a valid filter name')
    end
    if strcmp(filter_2, 'ekf')
        dataset_2 = run_ekf(config, measurements);
    elseif strcmp(filter_2, 'ukf')
        dataset_2 = run_ukf(config, measurements);
    else
        dataset_2 = [];
    end
end

%% Plotting & Results

% plot dstate from filter
f1 = figure(1);
plotter(f1, tspan, 'dstate', data, filter, filter_2, dataset, dataset_2)

% plot post-fit measurement residuals
f2 = figure(2);
plotter(f2, tspan, 'resid', data, filter, filter_2, dataset, dataset_2)

% Condition number of posteriori covariance vs time
if strcmp(filter, 'ckf') || strcmp(filter, 'potter') || strcmp(filter, 'srif')
    f13 = figure(13);
    plotter(f13, tspan, 'cond', data, filter, filter_2, dataset, dataset_2)
end

% plot trace of state and velocity covariance
f14 = figure(14);
plotter(f14, tspan, 'trace', data, filter, filter_2, dataset, dataset_2)

% plot ckf vs. srif
if strcmp(filter_2, 'ukf')
    f15 = figure(15);
    f16 = figure(16);
    f17 = figure(17);
    f18 = figure(18);
    plotter([f15, f16, f17], tspan, 'compare_ckf_other', data, filter, filter_2, dataset, dataset_2)
    plotter(f18, tspan, 'compare_cond', data, filter, filter_2, dataset, dataset_2)
end


