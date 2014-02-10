%% header 
% purpose: solve two-phase plug flow model
% for Fischer-Tropsch reactor
% camilla.berge.vik@ntnu.no, 19.09.2013
% camilla.berge.vik@ntnu.no, 19.12.2013
% camilla.berge.vik@ntnu.no, 07.01.2014

%% clear memory
clc
clear all
close all

%% declare global variables
global nCompGas nCompLiq gasConst ...
       gasDensityInit liqDensityInit ...
       tempGasInit tempSluInit...
       
%% get parameters
parameters
        
%% assemble initial conditions vector
y0 = [wtFracGasInit' wtFracLiqInit' supVelGasInit supVelLiqInit ...
     tempGasInit tempSluInit];

%% set height 
heightValues    = [0:0.1:totalHeight];  
nHeightValues   = length(heightValues);   

%% call ode15s
[z,y] = ode15s('model_equations',heightValues,y0);

%% unpack results
wRES_G      = y(:,1:nCompGas);               
wRES_L      = y(:,nCompGas+1:nCompGas+nCompLiq); 
vRES_G      = y(:,nCompGas+nCompLiq  + 1);    
vRES_L      = y(:,nCompGas+nCompLiq  + 2); 
tempGas     = y(:,nCompGas+nCompLiq  + 3); 
tempSlu     = y(:,nCompGas+nCompLiq  + 4);

%% outlet properties
molFracGasOut   = wRES_G(end,:)'./Mw/(sum(wRES_G(end,:)'./Mw));
avMolMassGasOut = molFracGasOut'*Mw;
gasDensityOut   = pTot.*avMolMassGasOut./gasConst./tempGas(end);

%% mass conservation calculations
fluxIn = gasDensityInit*supVelGasInit + ...
          liqDensityInit*supVelLiqInit % mass flux in feed

fluxOut = gasDensityOut*vRES_G(end) + ...
           liqDensityInit*vRES_L(end) % mass flux out of the reactor

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
     z,wRES_G(:,6),'b-', ...
     z,wRES_G(:,7),'k:' )
title('gas mass fractions')
legend('CO','H_2','H2O','C1-C10','C11-C20','C21-C30','C31+')
ylabel('gas phase mass fractions')
xlabel('reactor length')
grid on

figure(2)
plot(z,wRES_L(:,1),...
     z,wRES_L(:,2),...      
     z,wRES_L(:,3),...
     z,wRES_L(:,4),'g--',...
     z,wRES_L(:,5),'r-.',...
     z,wRES_L(:,6),'b-', ...
     z,wRES_L(:,7),'k:' )
title('liquid mass fractions')
legend('CO','H_2','H2O','C1-C10','C11-C20','C21-C30','C31+')
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


figure(9)
plot(z,tempGas)
title('Gas Temperature')
xlabel('Reactor Length [m]')
ylabel('Temperature [K]')

figure(10)
plot(z,tempSlu)
title('Slurry Temperature')
xlabel('Reactor Length [m]')
ylabel('Temperature [K]')
