function [x_obs, meas] = get_x_obs_meas(station, current_stations_state_eci, current_meas, data)
%GET_X_OBS get the current x_obs vector based on the station number,
%requires also a transformation matrix from eci to ecef in order to
%calculate the station velocity because it is not kept track of in the
%state
% 
% Inputs
%   station - current station # (1,2, or 3)
%   current_traj - 6x1 urrent trajectory being kept track of in kf
%   current_meas - 6x1 current true measurement
%   data - string both, rho, or rhod
%   t - scalar current time
%   w_E - scalar rotation rate of the Earth
% Ouputs
%   x_obs - 6x1 current station state
%   meas - 2x1 current measurement at specified station
%-------------------------------------------------------------------------%

switch station
    case 1
        x_obs = [current_stations_state_eci(1:3,1); current_stations_state_eci(4:6,1)];
        meas = current_meas(1:2);
    case 2
        x_obs = [current_stations_state_eci(1:3,2); current_stations_state_eci(4:6,2)];
        meas = current_meas(3:4);
    case 3
        x_obs = [current_stations_state_eci(1:3,3); current_stations_state_eci(4:6,3)];
        meas = current_meas(5:6);
end
% if only want one kind of data
switch data
    case 'rho'
        meas = meas(1);
    case 'rhod'
        meas = meas(2);
end
end