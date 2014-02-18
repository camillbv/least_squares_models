function y=LagrangeGLL2(n,x,N,xp)
% Lagrange(n,x,N)
% n= polynomial 0,1,...,N+1
% x= coordinate -1=< x =< 1
% N= polynomial order 
%
% Dorao, C.A.  2004

y=ones(length(x),1);

% [xp,w]=GaussLobattoLegendre(N);
ini=1;
iend=length(x);
%if x(1)==-1 & n==1
%    y(1)=1;
%    ini=2;
%    iend=length(x);
%end
%if x(length(x))==1 & n==N+1
%    y(length(x))=1;
%    ini=1;
%    iend=length(x)-1;
%end   
for i=ini:iend
        aux=((N-1)*(N).*(x(i)-xp(n))*LegendrePoly(N-1,xp(n)));
        aux1=(x(i)-1).*(x(i)+1)*LegendreDerivative(N-1,x(i));

        if abs(aux*10^-2) < eps  
            if abs(aux1*10^-2)<eps
                y(i)=1;
                %else
                %fprintf('error !!')
            end
         else
                y(i)=aux1./...
                aux;
                %               (N*(N+1).*(x(i)-xp(n))*Legendre(N,xp(n)));
                
         end
            
end

%only for N=2
%if N==2 & n==2
%    for j=1:length (x)
%        if x(j)==0    
%            y(j)=1;    
%        end
%    end
%end
    