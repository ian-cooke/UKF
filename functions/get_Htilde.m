function Htilde = get_Htilde(x_sc, x_obs, data)
%GET_HTILDE get 2x18 Htilde matrix
%
% Inputs
%   x_sc - current spacecraft state
%   x_obs - current station state
%   params - necessary paramater struct
%   Q - eci to ecef transformation matrix
%   data - if you want to only use one type of data
Htilde = Htilde_sc_rho_rhod(x_sc, x_obs, data); 
switch data
    case 'rho'
        Htilde = Htilde(1,:);
    case 'rhod'
        Htilde = Htilde(2,:);
end

end

