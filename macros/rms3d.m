function RMS = rms3d(error)
%RMS3D calculate 3d rms error (rms of element-wise norm of error)
%
% Inputs
%   errror - 3xN vector of 3d error
% Outputs
%   RMS - scalar RMS value of norms
%

% loop over N
error_norm = zeros(1, size(error, 2));
for i = 1:size(error, 2)
    error_norm(i) = norm(error(:, i));
end
% calculate rms
RMS = rms(error_norm);
end