function [RMS_comp_pos, RMS_comp_vel] = rmscomp(vec)
%RMSCOMP Get component wise state vector error residuals
RMS_comp_pos = zeros(3, 1);
RMS_comp_vel = RMS_comp_pos;
for i = 1:6
    if i < 4
        RMS_comp_pos(i) = rms(vec(i,:));
    else
        RMS_comp_vel(i) = rms(vec(i-3,:));
    end
end
end
        