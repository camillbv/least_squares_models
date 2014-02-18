function [p,w] = GLL_(N, x_min, x_max)
% Name      : GLL
%               Compute the Gauss Lobatto Legendre quadrature rule for the interval [x_min, x_max]
% Author    : Carlos A. Dorao cadorao@nt.ntnu.no   2005
% Location  : <directory in path>/@LSQ_0D_T

% Syntaxis  : [p,w]= GLL(N, x_min, x_max)
%               N             : Polynomial order
%               x_min         : Minimum size 
%               x_max         : Maximum size 
%              
% Post      : Return the quadrature points and weigth  
%               p            : vector of the quadrature points
%               w            : vector of the quadrature weigth
%
% See also  :  GaussLegendre

[p,w]=GaussLobattoLegendre(N);
delta=x_max-x_min;
for i=1:N
    p(i)=delta/2*(p(i)+1)+x_min;   % mapping from (-1,1) -> (x_min, x_max)
    w(i)=delta/2*w(i);
end
