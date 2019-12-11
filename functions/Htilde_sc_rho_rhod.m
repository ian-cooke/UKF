function Htilde = Htilde_sc_rho_rhod(x_sc, x_obs, data)
%Htilde_sc_rho_rhod Summary of this function goes here

% Vectors
R = x_sc(1:3);
V = x_sc(4:6);
R_s = x_obs(1:3);
V_s = x_obs(4:6);

% derivatives w range
drho_dR = (R - R_s)/norm(R - R_s);
drho_dV = zeros(1,3);

% derivatives w range rate
drhod_dR = -((R - R_s)'*(V - V_s))*(R - R_s)/norm(R - R_s)^3 + (V - V_s)/norm(R - R_s);
drhod_dV = (R - R_s)/norm(R - R_s);

% assemble
Htilde = [drho_dR', drho_dV; drhod_dR', drhod_dV'];

% If you only want one kind of data
switch data
    case 'rho'
        Htilde = Htilde(1,:);
    case 'rhod'
        Htilde = Htilde(2,:);
end
end

