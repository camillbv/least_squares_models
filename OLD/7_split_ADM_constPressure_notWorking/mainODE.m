%% header 
% purpose: solve gas plug flow model 
% for Fischer-Tropsch reactor
% camilla.berge.vik@ntnu.no
% 28.01.2014

%% clear memory
clc
clear all
close all

%% declare global variables
global nCompGas nCompLiq 
       
%% set parameters and initial conditions
parameters

%% assemble initial conditions vector
y0 = [wtFracGasInit' wtFracLiqInit' supVelGasInit supVelLiqInit];

%% set height span
heightSpan = [0 totalHeight];

%% call ode15s
[z,y] = ode15s('model_equations',heightSpan,y0);

%% unpack results
wRES_G = y(:,1:nCompGas);               
wRES_L = y(:,nCompGas+1:nCompGas+nCompLiq); 
vRES_G = y(:,nCompGas+nCompLiq + 1);    
vRES_L = y(:,nCompGas+nCompLiq + 2); 

%% make plots
figure(1)
plot(z,wRES_G(:,1),...
     z,wRES_G(:,2),...      
     z,wRES_G(:,3),...
     z,wRES_G(:,4),'g--',...
     z,wRES_G(:,5),'r-.',...
     z,wRES_G(:,6),'b-' )
title('gas mass fractions')
legend('CO','H_2','H2O','C1-C20','C21-C30','C31+')
ylabel('gas phase mass fractions')
xlabel('reactor length')
grid on

figure(2)
plot(z,wRES_L(:,1),...
     z,wRES_L(:,2),...      
     z,wRES_L(:,3),...
     z,wRES_L(:,4),'g--',...
     z,wRES_L(:,5),'r-.',...
     z,wRES_L(:,6),'b-' )
title('liquid mass fractions')
legend('CO','H_2','H2O','C1-C20','C21-C30','C31+')
ylabel('liquid phase mass fractions')
xlabel('reactor length')
grid on

%% calculate conversion of CO (%) on mass basis
wconversion     = 100*wRES_G(:,1);
wcum_conversion = 100*(wtFracGasInit(1) - wRES_G(:,1))./wtFracGasInit(1);

figure(3)
plot(z,wconversion,z,wcum_conversion,'r')
title('Conversion on mass basis (%)')
ylabel('conversion (%)')
xlabel('reactor length')
legend('Conversion','Cumulative Conversion')

figure(4)
plot(z,vRES_G)
title('Superficial gas velocity')

figure(5)
plot(z,vRES_L)
title('Superficial liquid velocity')
xlabel('reactor length')
ylabel('liquid velocity')
