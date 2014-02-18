function y=Lagrange(n,x,N)
% Lagrange(n,x,N)
% n= polynomial 0,1,...,N+1
% x= coordinate -1=< x =< 1
% N= polynomial order 
%
% Dorao, C.A.  2004

y=ones(length(x),1);

[xp,w]=GaussLobattoLegendre(N+1);
ini=1;
iend=length(x);
if x(1)==-1 & n==1
    y(1)=1;
    ini=2;
    iend=length(x);
end
if x(length(x))==1 & n==N+1
    y(length(x))=1;
    ini=1;
    iend=length(x)-1;
end   
for i=ini:iend
        aux=(N*(N+1).*(x(i)-xp(n))*LegendreP(N,xp(n)));
        if aux==0
            aux1=(x(i)-1).*(x(i)+1)*LegendreDerivative(N,x(i));
            if aux1==0
                y(i)=1;
                %else
                %fprintf('error !!')
            end
         else
                y(i)=(x(i)-1).*(x(i)+1)*LegendreDerivative(N,x(i))./...
                aux;
                %               (N*(N+1).*(x(i)-xp(n))*LegendreP(N,xp(n)));
                
         end
            
end

%only for N=2
if N==2 & n==2
    for j=1:length (x)
        if x(j)==0    
            y(j)=1;    
        end
    end
end
    