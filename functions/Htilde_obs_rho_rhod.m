function Htilde = Htilde_obs_rho_rhod(x_sc, x_obs, data, w_e)
%UNTITLED10 Summary of this function goes here
x = x_sc(1);
y = x_sc(2);
z = x_sc(3);
xd = x_sc(4);
yd = x_sc(5);
zd = x_sc(6);
xs = x_obs(1);
ys = x_obs(2);
zs = x_obs(3);
xds = x_obs(4);
yds = x_obs(5);
zds = x_obs(6);
Htilde = [                                                                                                                                                     -(abs(x - xs)*sign(x - xs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2),                                                                                                                                                     -(abs(y - ys)*sign(y - ys))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2),                                                                                                                           -(abs(z - zs)*sign(z - zs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2);...
           (abs(x - xs)*sign(x - xs)*((xd + w_e*ys)*(x - xs) + (yd - w_e*xs)*(y - ys) + zd*(z - zs)))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - (xd + w_e*ys + w_e*(y - ys))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2), (w_e*xs - yd + w_e*(x - xs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2) + (abs(y - ys)*sign(y - ys)*((xd + w_e*ys)*(x - xs) + (yd - w_e*xs)*(y - ys) + zd*(z - zs)))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2), (abs(z - zs)*sign(z - zs)*((xd + w_e*ys)*(x - xs) + (yd - w_e*xs)*(y - ys) + zd*(z - zs)))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - zd/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2)];       
% if you only want one type of data
switch data
    case 'rho'
        Htilde = Htilde(1,:);
    case 'rhod'
        Htilde = Htilde(2,:);
end
end

