%% header 
% purpose: solve gas plug flow model 
% for Fischer-Tropsch reactor
% camilla.berge.vik@ntnu.no
% 23.09.2013

%% clear memory
clc
clear all
close all

%% declare global variables
global Cmax nCompGas nCompLiq gasConst T pTot...
       gasDensityInit liqDensityInit nu...
       
%% get parameters
parameters

%% inlet gas phase concentrations
molFracGasCOInit    = 0.3;   
molFracGasH2Init    = 0.7; 
molFracGasH2OInit   = 0.0;     
molFracGasAlkInit   = zeros(Cmax+1,1)';
molFracGasInit      = [molFracGasCOInit molFracGasH2Init molFracGasH2OInit molFracGasAlkInit]';
molFracGasInit      = molFracGasInit./sum(molFracGasInit);  % scale it to sum to 1

%% inlet liquid phase concentrations
molFracLiqCOInit    = 0;
molFracLiqH2Init    = 0;
molFracLiqH2OInit   = 0;
molFracLiqAlkInit   = nu(4:length(nu)).*ones(Cmax+1,1);    % prod. gen. by one kmole reacted CO
molFracLiqAlkInit(1:20)=0.00001; 
molFracLiqInit      = [molFracLiqCOInit molFracLiqH2Init molFracLiqH2OInit molFracLiqAlkInit']';     % add water
molFracLiqInit      = molFracLiqInit./sum(molFracLiqInit);  % scale it to sum to 1

%% convert to weight fractions
wtFracGasInit   = molFracGasInit.*Mw./(molFracGasInit'*Mw); %  kg/kg (unitless)  gas feed weight fractions
wtFracLiqInit   = molFracLiqInit.*Mw./(molFracLiqInit'*Mw); %  kg/kg (unitless)  gas feed weight fractions

%% number of components
nCompGas = length(wtFracGasInit);
nCompLiq = length(wtFracLiqInit);

%% feed average molar mass
avMolMassGasInit   = molFracGasInit'*Mw; % kg/kmol gas feed average molar mass

%% feed gas density
gasDensityInit    = pTot*avMolMassGasInit/gasConst/T; % kg/m^3  gas feed total density
liqDensityInit    = 600; % estimate based on HYSYS flash simulations

%% assemble initial conditions vector
y0 = [wtFracGasInit' wtFracLiqInit' supVelGasInit supVelLiqInit];

%% set height 
heightValues    = [0:0.1:totalHeight];  
nHeightValues   = length(heightValues);   
heightSpan = [0 totalHeight];

%% call ode15s
[z,y] = ode15s('model_equations',heightSpan,y0);

%% unpack results
wRES_G = y(:,1:nCompGas);               
wRES_L = y(:,nCompGas+1:nCompGas+nCompLiq); 
vRES_G = y(:,nCompGas+nCompLiq + 1);    
vRES_L = y(:,nCompGas+nCompLiq + 2); 

for iz = 1:size(wRES_G,1)
    % gas
    wRES_C1toC10_G(iz) =  sum(wRES_G(iz,4:13));
    wRES_C11toC20_G(iz) = sum(wRES_G(iz,14:23));
    wRES_C21toC30_G(iz) = sum(wRES_G(iz,24:33));
    wRES_C31plus_G(iz)  = sum(wRES_G(iz,34));
    
    % liquid
    wRES_C1toC10_L(iz)  = sum(wRES_L(iz,4:13));
    wRES_C11toC20_L(iz) = sum(wRES_L(iz,14:23));
    wRES_C21toC30_L(iz) = sum(wRES_L(iz,24:33));
    wRES_C31plus_L(iz)  = sum(wRES_L(iz,34));
end  

figure(1)
plot(z,wRES_G(:,1),...
     z,wRES_G(:,2),...      
     z,wRES_G(:,3),...
     z,wRES_C1toC10_G,'g--',...
     z,wRES_C11toC20_G,'r-.',...
     z,wRES_C21toC30_G,'b-', ...
     z,wRES_C31plus_G,'k:' )
title('gas mass fractions')
legend('CO','H_2','H2O','C1-C10','C11-C20','C21-C30','C31+')
ylabel('gas phase mass fractions')
xlabel('reactor length')
grid on

figure(2)
plot(z,wRES_L(:,1),...
     z,wRES_L(:,2),...      
     z,wRES_L(:,3),...
     z,wRES_C1toC10_L,'g--',...
     z,wRES_C11toC20_L,'r-.',...
     z,wRES_C21toC30_L,'b-', ...
     z,wRES_C31plus_L,'k:' )
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


figure(6)
plot(z,vRES_L*liqDensityInit,'r')
title('Liquid mass flux')
ylabel('liquid mass flux (v_LS*rho_L [kg/(m^2s)])')
xlabel('reactor length')