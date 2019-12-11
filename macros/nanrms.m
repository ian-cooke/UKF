function RMS = nanrms(A)
%NANRMS calculate the rms value row-wise of a matrix ignoring NaN values
m1 = 0;
m2 = 0;
RMS_1 = 0;
RMS_2 = 0;
[r,c] = size(A);
for i = 1:r
    for j = 1:c
        if ~isnan(A(i,j))
            switch i
                case {1,3,5}
                    RMS_1 = RMS_1 + A(i,j)^2;
                    m1 = m1 + 1;
                case {2,4,6}
                    RMS_2 = RMS_2 + A(i,j)^2;
                    m2 = m2 + 1;
            end
        end
    end
end
RMS_1 = sqrt(1/m1*RMS_1);
RMS_2 = sqrt(1/m2*RMS_2);
RMS = [RMS_1; RMS_2];
end

