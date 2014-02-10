function Ld = LegendreDerivative(n,x);
% LegendreDerivative(n,x)
% This function returns the value of the n'th derivative of the 
% Legendre polynomial evaluated in the point, x.
% n= derivative order 
% x= coordinate 
% Dorao, C.A. 
%
Ln = zeros(n+1,1);

%if length(x)==1
    Ln(1) = 1;
    Ln(2) = x;

    Ld(1) = 0;
    Ld(2) = 1;

    % Have to treat the endpoints separatly
    if abs(x)==1
      Ld = x^(n-1)*(1/2)*n*(n+1);
    else
      for i = 1:n-1
        Ln(i+2) = (2*i+1)/(i+1)*x*Ln(i+1) - i/(i+1)*Ln(i);
      end
      Ld = n/(1-x^2)*Ln(n) - n*x/(1-x^2)*Ln(n+1);
    end
%else
%    
%     der=x*0.0;
%     for idx=1:length(x)
%         Ln = zeros(n+1,1);
%         Ln(1) = 1;
%         Ln(2) = x(idx);
% 
%         Ld(1) = 0;
%         Ld(2) = 1;
% 
%         % Have to treat the endpoints separatly
%         if abs(x(idx))==1
%           Ld = x(idx)^(n-1)*(1/2)*n*(n+1);
%         else
%           for i = 1:n-1
%             Ln(i+2) = (2*i+1)/(i+1)*x(idx)*Ln(i+1) - i/(i+1)*Ln(i);
%           end
%           Ld = n/(1-x(idx)^2)*Ln(n) - n*x(idx)/(1-x(idx)^2)*Ln(n+1);
%         end
%         
%         der(idx)=Ld;
%     end%endfor idx
% 
%     Ln=zeros(length(x),n+1);
%     Ln(:,1)=ones(length(x),1);
%     Ln(:,2)=x;
%     
%     Ld=zeros(length(x),n+1);
%     %Ld(:,1)=0 0 0 0 ... 0;
%     Ld(:,2)=ones(length(x),1);
%     
%     for i=1:n-1
%         Ln(:,i+2)=(2*i+1)/(i+1)*x.*Ln(:,i+1) - i/(i+1)*Ln(:,i);
%     end%endfor i
%     Ld= n./(1-x.^2).*Ln(:,n) - n*x./(1-x.^2).*Ln(:,n+1);
%     
%     I=find(abs(x==1));
%     if I>0
%        Ld(I) = x(I).^(n-1)*(1/2)*n*(n+1);
%     end
end %endfor length