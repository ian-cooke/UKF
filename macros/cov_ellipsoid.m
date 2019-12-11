function [xyz] = cov_ellipsoid(P)
%COV_ELLIPSOID Calcualtes xyz coordinates for surf() function to plot
%covariance ellipsoid centered about mu = [0,0,0]

% eigenvectors & values
[vecs, vals] = eig(P);

% # of standard deviations
N = 3;
radii = N*sqrt(diag(vals));

% data for unrotated ellipsoid
[xu, yu, zu] = ellipsoid(0,0,0, radii(1), radii(2), radii(3));

% rotate data using kronecker tensor product 
a = kron(vecs(:, 1), xu);
b = kron(vecs(:, 2), yu);
c = kron(vecs(:, 3), zu);

% assemble data
data = a + b + c;
n = size(data, 2);

% compute new rotated coordinates (still centered at zero)
x = data(1:n,  :);
y = data(n+1:2*n, :);
z = data(2*n+1:end, :);

xyz.x = x;
xyz.y = y;
xyz.z = z;



end

