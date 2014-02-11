function L = Legendre(n,x);

% This function returns the value of the n'th Legendre polynomial
% evaluated in the point, x.
%
% authors: Kay Hansen-Zahl
%          Tormod Bjøntegaard
% 
%          Fall 2002

Ln = zeros(n+1,1);

Ln(1) = 1;
Ln(2) = x;
if n>1
  for i = 1:n-1
    Ln(i+2) = (2*i+1)/(i+1)*x*Ln(i+1) - i/(i+1)*Ln(i);
  end
end
L = Ln(n+1);