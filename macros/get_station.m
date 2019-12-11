function station = get_station(meas)
%GET_STATION get the integer station number (1,2,3) of the measurement
%currently being observed
% 
% Inputs
%   meas - current measurement
% Outputs
%   station - station number

station = 0;
for i = 1:3
    if ~isnan(meas(2*i))
        station = i;
    end
end

end