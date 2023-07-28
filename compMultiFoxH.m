function result = compMultiFoxH(params, nsubdivisions)
% This module estimates a multivariate integral using simple rectangule quadrature. 
% In most practical applications, 20 points per dimension provide sufficient accuracy.

% Inputs:
%  'params': list containing z, mn, pq, c, d, a, b.
%  'nsubdivisions': the number of divisions taken along each dimension. 
%   Note that the total number of points will be nsubdivisions**dim.
%  'boundaryTol': tolerance used for determining the boundaries
% Output:
%  'result': the estimated value of the multivariate Fox H function.

boundaryTol = 0.0001;   boundaries = detBoundaries(params, boundaryTol);
dim = size(boundaries,2);
signs = dec2bin(2^dim-1:-1:0)-'0';  signs =  2*signs-1;
X = zeros(round(nsubdivisions / 2)^dim,dim);    A = (1:round(nsubdivisions / 2));
for i = 1:size(X,1)
    ixVect = ind2subVect(repmat(round(nsubdivisions / 2),1,dim),i); 
    for j = 1:dim
        X(i,j) = A(ixVect(j));
    end
end
code = sortrows(X); coder.varsize('quad');
quad = 0i;  res = zeros((0));
for k = 1: size(signs,1)
    points = signs(k,:) .* (code - 0.5) .* boundaries .* 2 ./ nsubdivisions;
    res = [res, real(compMultiFoxHIntegrand(points, params))];
    quad = quad + sum(compMultiFoxHIntegrand(points, params));
end
volume = prod(2 * boundaries / nsubdivisions);
result = quad * volume;
end