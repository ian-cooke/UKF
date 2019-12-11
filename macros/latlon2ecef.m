function stations_ecef = latlon2ecef(n_stations, lat_lon, R)
%LATLON2ECEF Convert Latitude and longitude of array of station coordinates
%(assumed to be at sea level) to ecef
% 
% Inputs
%   n_stations - number of stations in array
%   lat_lon - n_stations x 2 array of latitude/longitude
%   R - radius of body
% Outputs
%   stations_ecef - ecef coordinates of stations
stations_ecef = zeros(3, n_stations);
for i = 1:n_stations
    lat = deg2rad(lat_lon(i, 1));
    lon = deg2rad(lat_lon(i, 2));
    s_x = R*cos(lat)*cos(lon);
    s_y = R*cos(lat)*sin(lon);
    s_z = R*sin(lat);
    s_i_ecef = [s_x; s_y; s_z];
    stations_ecef(:, i) = s_i_ecef;
end
end

