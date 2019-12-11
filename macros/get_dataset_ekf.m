function dataset = get_dataset_ekf(xhat, Phat, resid, mat_cond, config)
%GET_DATASET_CKF
% Creates dataset for the ckf

tspan = config.tspan;
n = config.n;
x_true = config.x_true(1:n, :);

% Assemble vector of vectors that only record the state/ covariance when there are
% measurements
% also calculate trace of position and velocity covariance, also translate
% station positions to ecef
% estimated state & nan
xhat_nan = NaN(n, length(tspan));
% trace of position and velocity covariance
trace_pos = zeros(1,length(tspan));
trace_vel = zeros(1,length(tspan));
trace_pos_nan = NaN(1,length(tspan));
trace_vel_nan = NaN(1,length(tspan));
% covariance nan
Phat_nan = NaN(n, n, length(tspan));
% x_true_nan for deltax
x_true_nan = NaN(n, length(tspan));
% condition number nan
mat_cond_nan = NaN(1,length(tspan));
for k = 1:length(tspan)
    % traces
    trace_pos(k) = trace(Phat(1:3, 1:3, k));
    trace_vel(k) = trace(Phat(4:6, 4:6, k));
    % NaN vectors are defined as spots where measurements are unavailable,
    % it makes it so that the plots are very clear where measurements are
    % unavailable.
    if all(isnan(resid(:,k)))
        Phat_nan(:, :, k) = Phat(:, :, k);
        trace_pos_nan(k) = trace_pos(k);
        trace_vel_nan(k) = trace_vel(k);
        xhat_nan(:, k) = xhat(:, k);
        x_true_nan(:, k) = x_true(:, k);
        mat_cond_nan(k) = mat_cond(k);
    end
end

% Calculate covariance ellipse coordinates at the final time
xyz_ellipse_r = cov_ellipsoid(Phat(1:3, 1:3, end));
xyz_ellipse_v = cov_ellipsoid(Phat(4:6, 4:6, end));

% post-fit residual RMS values
RMS_resid = nanrms(resid);
RMS_3d_pos = rms3d(x_true(1:3,:) - xhat(1:3,:));
RMS_3d_vel = rms3d(x_true(4:6,:) - xhat(4:6,:));
[RMS_comp_pos, RMS_comp_vel] = rmscomp(x_true - xhat);

% Assemble dataset from return ckf vars
dataset.xhat = xhat;
dataset.xhat_nan = xhat_nan;
dataset.dxhat = x_true - xhat;
dataset.dxhat_nan = x_true_nan - xhat_nan;
dataset.Phat = Phat;
dataset.Phat_nan = Phat_nan;
dataset.resid = resid;
dataset.trace_pos = trace_pos;
dataset.trace_vel = trace_vel;
dataset.trace_pos_nan = trace_pos_nan;
dataset.trace_vel_nan = trace_vel_nan;
dataset.mat_cond = mat_cond;
dataset.mat_cond_nan = mat_cond_nan;
dataset.xyz_ellipse_r = xyz_ellipse_r;
dataset.xyz_ellipse_v = xyz_ellipse_v;
dataset.RMS_resid = RMS_resid;
dataset.x_true = x_true;
dataset.RMS_3d_pos = RMS_3d_pos;
dataset.RMS_3d_vel = RMS_3d_vel;
dataset.RMS_comp_pos = RMS_comp_pos;
dataset.RMS_comp_vel = RMS_comp_vel;

end