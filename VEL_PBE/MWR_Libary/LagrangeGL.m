function y=LagrangeGL(n,x,N)
% Lagrange Interpolant Polynomial based on GL points 
% y=LagrangeGL(n,x,N)
% n= polynomial 1,...,N
% x= coordinate -1=< x =< 1
% N= polynomial order 
%
% Dorao, C.A.  2004

y=ones(length(x),1);

%[xp,w]=GaussLegendre(N);
[xp,w]=GL(N,-1,1);

ini=1;
iend=length(x);

for i=1:length(x)
                aux =( (x(i)-xp(n))*LegendreDerivative(N,xp(n)) );
                aux1= LegendreP(N,x(i));
                if abs(aux*10^-2) < eps  & abs(aux1*10^-2)<eps
                    y(i)=1;
                    
                else
                    y(i)=aux1/ aux;
                end
            
end

if N==2 & n==2
    for j=1:length (x)
        if x(j)==0    
            y(j)=1;    
        end
    end
end
    