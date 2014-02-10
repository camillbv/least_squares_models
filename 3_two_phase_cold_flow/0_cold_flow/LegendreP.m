function L = LegendreP(n,x);
% Legendre(n,x)
% n= polynomial order 0,1,...
% x= coordinate -1=< x =< 1
%
% Dorao, C.A.  2004

Ln = zeros(n+1,length(x));
Ln(1,:) = ones(1,length(x));
Ln(2,:) = x;
if n>1
  for i = 1:n-1
    Ln(i+2,:) = (2*i+1)/(i+1)*x.*Ln(i+1,:) - i/(i+1)*Ln(i,:);
  end
end
L = Ln(n+1,:);

