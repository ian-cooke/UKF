function dzdt = keplerJ2J3_wPhi_ODE(t, z, params)
%KEPLER_J2_DRAG_WPHI_ODE Primary ode integration function for project 1
% Uses standard gravity w/ J2 & J3 perturbations
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
J3 = params.J3;

% Extract from z
x = z(1:n);
R = x(1:3);
V = x(4:6);
Phi_flat = z(n+1:end);
Phi = reshape(Phi_flat, n, n);

%% State vector change
% change in state due to gravity
r = norm(x(1:3));
rd = V;
rdd = -mu*R/r^3;

r = norm(R);
X = R(1);
Y = R(2);
Z = R(3);

% Perturbation due to J2
p_J2 = 3/2*J2*mu/r^2*(R_E/r)^2*[X/r*(5*Z^2/r^2 - 1), Y/r*(5*Z^2/r^2 - 1), Z/r*(5*Z^2/r^2 - 3)]';
p_J3 = [-(5*J3*R_E^3*X*Z*mu*(3*X^2 + 3*Y^2 - 4*Z^2))/(2*r^9);...
        -(5*J3*R_E^3*Y*Z*mu*(3*X^2 + 3*Y^2 - 4*Z^2))/(2*r^9);...
         (J3*R_E^3*mu*(3*X^4 + 6*X^2*Y^2 - 24*X^2*Z^2 + 3*Y^4 - 24*Y^2*Z^2 + 8*Z^4))/(2*r^9)];
 

% state vector change
dxdt = [rd;rdd + p_J2 + p_J3];


%% STM Change
A = dfdx_wJ2J3(x(1:3), params);
Phid = A*Phi;
Phid_flat = reshape(Phid, n^2, 1);

%% Assemble
dzdt = [dxdt; Phid_flat];

end

