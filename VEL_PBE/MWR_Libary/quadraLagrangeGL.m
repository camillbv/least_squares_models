function quadraLagrangeGL(N)

%used for ploting the exact solution
vecX=-1:.01:1;

%p_sample are the points to sample the function
[p_sample,q]=GLL_Int(N,-1,1);

%functionY contain the nodal points
%f(x)=sum( functionY(i)*Li(x) ) from i=1 to N
functionY=feval('func',p_sample);

%example of the integration by nodal approximation
[p,q]=GLL_Int(N,p_sample(2),1); %points and weigth in the interval (0,1)
L=lagrangeMatrix(N,p_sample(2),1); %Lagrange matrix

figure(1)
plot(vecX,feval('func',vecX),'--',p_sample,functionY,'o',p,functionY'*L,'<:')
legend('exact','sample points','Interpolation function')

%integral using the Interpolation
integralNystrom=functionY'*L*q;

%exact Integration
integralExact=-cos(1)+cos(p_sample(2));

abs(integralNystrom-integralExact);
error=zeros(length(p_sample),1);

for i=2:length(p_sample)-1
    fprintf('i: %d \n',i)
    [p,q]=GLL_Int(N,p_sample(i),1); %points and weigth in the interval (0,1)
    L=lagrangeMatrix(N,p_sample(i),1); %Lagrange matrix
    L(i,1)=1;
    
%    %integral using the Interpolation
    integralNystrom=functionY'*L*q;

%    %exact Integration
    integralExact=-cos(1)+cos(p_sample(i));
    
    er=100*abs(integralNystrom-integralExact)/integralExact;
    error(i,1)=er;
end
p_sample
error(1)
figure(2)
semilogy(p_sample,error,'ro:')
title(['quadrature points N= ' num2str(N)])
xlabel('x')
ylabel('percentual error')
axis([-1 1 10E-19 10E1])


function y=h(xtr,x)
    y=x;
    
    for i=1:length(y(:,1))
        
        if (x(i)>xtr)
            y(i)=1;
        else
            y(i)=0;
        end
        
    end
    %y
    
    

function y=func(x)
    y=sin(x);
    
function [p,q]=GLL_Int(N,a,b)
    [p,q]=GaussLobattoLegendre(N);
    p=(b-a)/2.*(p+1)+a;
    q=q*(b-a)/2;

    
function L=lagrangeMatrix(N,a,b)
% p vector of points to construct the Lagrange Matrix
% N the order of the matrix 

    if (a==-1 & b==1)
        id=ones(N,1);
        L=diag(id);

    else
        [p,q]=GLL_Int(N,a,b);
        L=zeros(length(p),N);
        for i=1:N
            L(i,:)=Lagrange(i,p,N-1)';
        end

        for i=1:N
            for j=1:N
                if  isinf(L(i,j))==1
                    L(i,j)=1;
                end
            end
        end 
        L(N,N)=1;
    end
    
    %figure(3)
%plot (p,L(1,:))
%figure(4)
%p=-1:.01:1;
%plot(p,Lagrange(1,p,N-1))
