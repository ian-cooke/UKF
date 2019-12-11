function S = inv_back_recursion(R)
%INV_BACK_RECURSION Compute the inverse of a square upper triangular matrix
%efficiently using backwards recursion
%
% Inputs
%   R - upper square triangular matrix to be inverted
% Outputs
%   S - inverse of R

n = size(R, 1);
S = zeros(n,n);
% diag elements
for i = 1:n
    % diag elements
    S(i,i) = 1/R(i,i); 
end
% off-diag elements
for i = 1:n
    for j = (i+1):n
        sum = 0;
        for k = i:(j-1)
            sum = sum + R(k,j)*S(i,k);
        end
        S(i,j) = -S(j,j)*sum;
    end
end

end

