function [R,V] = kep2eci(a,e,i,Omega,omega,f,mu)
%KEP2ECI Convert Keplerian orbit Elements to ECI state vector
%
% Inputs
%   a - semi-major axis [m]
%   e - eccentricity
%   i - eccentricity
%   Omega - right ascension of the ascending node [deg]
%   omega - argument of perigee [deg]
%   f - true anomaly [deg]
% Outputs
%   R - position vector [X,Y,Z]' in ECI
%   V - velocity vector [Xd,Yd,Zd]' in ECI

% params & conversions
i = deg2rad(i);
Omega = deg2rad(Omega);
omega = deg2rad(omega);
f = deg2rad(f);

% Define a few values
p = a*(1-e^2);
h = sqrt(mu*p);
V_r = h*e/p*sin(f);
r = p/(1+e*cos(f));
V_theta = h/r;

% Q matrix
Q_11 = cos(omega)*cos(Omega) - sin(omega)*sin(Omega)*cos(i);
Q_12 = -sin(omega)*cos(Omega) - cos(omega)*sin(Omega)*cos(i);
Q_21 = cos(omega)*sin(Omega) + sin(omega)*cos(Omega)*cos(i);
Q_22 = -sin(omega)*sin(Omega) + cos(omega)*cos(Omega)*cos(i);
Q_31 = sin(omega)*sin(i);
Q_32 = cos(omega)*sin(i);
Q = [Q_11, Q_12; Q_21, Q_22; Q_31, Q_32];

% Star Values
Xdstar = V_r*cos(f) - V_theta*sin(f);
Ydstar = V_r*sin(f) + V_theta*cos(f);

% Calculate the state & return
R = r*Q*[cos(f);sin(f)];
V = Q*[Xdstar;Ydstar];
end

