function [p,w] = GaussLegendre(np);
% Name      : GaussLegendre(np)
%               Compute the Gauss Legendre quadrature rule for the interval
%               [-1, 1]
%   np= number of quadrature
%   p = quadrature points
%   w = quadrature weigth 
%
% Author    : Carlos A. Dorao cadorao@nt.ntnu.no   2005
% Location  : <directory in path>/HighOrderLib

A = zeros(np,np);
p = zeros(np,1);
w = zeros(np,1);

% This loop finds the A-matrix
A(1,2) = 1;
if np>2
  for i = 2:np-1
    A(i,i-1) = (i-1)/(2*i-1);
    A(i,i+1) = i/(2*i-1);
  end
end  
A(np,np-1) = (np-1)/(2*np-1);


% The array of the sorted eigenvalues/zeros

p=sort(eig(A));

% This loop finds the associated weights
for j=1:np
  w(j) = 2/((1-p((j))^2)*(LegendreDerivative(np,p(j)))^2);
end
