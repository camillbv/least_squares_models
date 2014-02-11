function [p,w] = GaussLobattoLegendre(np)
% [p,w]=GaussLobattoLegendre(np)
% GaussLobattoLegendre quadrature
% np= number of points 
% p =quadrature points 
% q =quadrature weigth
%
% The tolerance in the Newton-iteration
% tol = 1E-14;

tol = 1E-14;
p = zeros(np,1);
p(1) = -1;
p(np) = 1;

w = zeros(np,1);

if np<3
  for i = 1:np
    L = LegendreP(np-1,p(i));
    w(i) = 2/((np-1)*np*L^2);
  end  
  return
end

% These points are needed for the startvalues in the Newton-iteration

[GLpoints,GLweights] = GaussLegendre(np-1);

startvalues = zeros(np,1); 
startvalues(2:np-1) = (GLpoints(1:np-2)+GLpoints(2:np-1))/2;

%This loop executes the Newton-iteration to find the GLL-points
for i = 2:np-1
  p(i) = startvalues(i);
  p_old = 0;
  while (abs(p(i)-p_old)) > tol
    p_old = p(i);
    L = LegendreP(np-1,p_old);
    Ld = LegendreDerivative(np-1,p_old);
    p(i) = p_old + ((1-p_old^2)*Ld)/((np-1)*np*L);
  end
end

% This loop finds the associated weights
for i = 1:np
  L = LegendreP(np-1,p(i));
  w(i) = 2/((np-1)*np*L^2);
end
