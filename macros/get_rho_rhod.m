function y_predicted = get_rho_rhod(x_sc, x_obs, data)
%GET_RHO_RHOD get expected range and range rate based on current spacecraft
%state vector and station vector

R = x_sc(1:3,1);
V = x_sc(4:6,1);
Rs = x_obs(1:3,1);
Vs = x_obs(4:6,1);
rho = norm(R - Rs);
rhod = (R - Rs)'*(V - Vs)/rho;

% if you want to ignore one type of data
switch data
    case 'both'
        y_predicted = [rho;rhod];
    case 'rho'
        y_predicted = rho;
    case 'rhod'
        y_predicted = rhod;
end

end