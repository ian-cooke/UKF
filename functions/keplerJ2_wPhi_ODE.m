function dzdt = keplerJ2_wPhi_ODE(t, z, params)
%KEPLER_J2_DRAG_WPHI_ODE Primary ode integration function for project 1
% Uses standard gravity w/ J2 and drag perturbations
%
% Inputs
%   t - current time
%   z - current state vector and phi matrix
%   params - struct of necessary params, including
%       R_E - radius of earth
%       w_E - spin rate of earth
%       J2 - first spherical harmonic
%       J3 - second spherical harmonic
% Outputs
%   dxdt - time rate of change of state vector
%-------------------------------------------------------------------------%

%% constants

% from params
R_E = params.R_E;
n = params.n;
mu = params.mu;
J2 = params.J2;

% Extract from z
x = z(1:n);
R = x(1:3);
V = x(4:6);
Phi_flat = z(n+1:end);
Phi = reshape(Phi_flat, n, n);

%% State vector change
% change in state due to gravity
r = norm(R);
rd = V;
rdd = -mu*R/r^3;

% transform to ECEF for J2 perturbation
X = R(1);
Y = R(2);
Z = R(3);

% Perturbation due to J2
p_J2 = 3/2*J2*mu/r^2*(R_E/r)^2*[X/r*(5*Z^2/r^2 - 1), Y/r*(5*Z^2/r^2 - 1), Z/r*(5*Z^2/r^2 - 3)]';

% state vector change
dxdt = [rd;rdd + p_J2];


%% STM Change
A = dfdx_wJ2(x(1:3), params);
Phid = A*Phi;
Phid_flat = reshape(Phid, n^2, 1);

%% Assemble
dzdt = [dxdt; Phid_flat];

end

