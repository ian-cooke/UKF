function dfdz_sym = dfdx_wJ2_sym()
%DFDX_WJ2J3
% state vector is [x, y, z, vx, vy, vz, mue, J2, J3]' in ECI
% xdot = [vx, vy, vz, ax, ay, az, 0, 0, 0] in ECI
syms x y z xd yd zd J2 mue r phi R_E real;
r = sqrt(x^2 + y^2 + z^2);
phi = asin(z/r);

% potentials in gravity and J2/J3 components
f_grav = -mue/r;
f_J2 = mue/r*(J2/2*(R_E/r)^2*(3*sin(phi)^2-1));
% Compute gravity acceleration in ECI
rddot_grav = -gradient(f_grav, [x,y,z]);
% Compute J2/J3 acceleration in ECI
rddot_J2 = -gradient(f_J2, [x,y,z]);
% sum the two to get the acceleration in ECI
rddot = rddot_grav + rddot_J2;

% put together xdot
f = [xd, yd, zd, rddot']'; 
% symbolic state vector
z_sym = [x, y, z, xd, yd, zd];
% symbolic jacobian
dfdz_sym = jacobian(f, z_sym);

end
