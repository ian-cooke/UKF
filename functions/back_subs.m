function x = back_subs(b, R)
%BACK_SUBS Backward substitution algorithm

x = zeros(length(b), 1);
for i = n:-1:1
    sum = 0;
    for j = i+1:n
        sum = sum + R(i,j)*x(j);
    end
    x(i) = (b(i) - sum)/R(i,i);
end