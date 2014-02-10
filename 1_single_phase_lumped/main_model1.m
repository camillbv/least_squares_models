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
global rho_G

%% run file to get parameters
parameters

%% set feed gas phase concentrations
x0CO_G     = 0.3;  % kg/kg (unitless) gas feed mole fraction of CO     
x0H2_G     = 0.7;  % kg/kg (unitless) gas feed mole fraction of H2    
x0H2O_G    = 0;    % kg/kg (unitless) gas feed mole fraction of H2O    
x0C16H34_G = 0;    % kg/kg (unitless) gas feed mole fraction of C16H34 
x0 = [ x0CO_G  x0H2_G  x0H2O_G  x0C16H34_G]';

%% convert to weight fractions
w0   = x0.*Mw./(x0'*Mw); %  kg/kg (unitless)  gas feed weight fractions

%% calculate feed average molar mass
Mw_ave0   = x0'*Mw; % kg/kmol gas feed average molar mass

%% calculate feed gas density
rho_G    = ptot*Mw_ave0/R/T; % kg/m^3  gas feed total density

%% assemble initial conditions vector
y0 = w0';

%% set zspan 
zstart          = 0;                % m
zend            = L;                % m
zspan           = [zstart zend];    % m

%% call ode15s
[z,y] = ode15s('model_equations',zspan,y0);

%% unpack results
wRES = y;

%% convert results from weight to mole fractions
xRES = zeros(size(wRES));

for zstep = 1:length(z)
    xRES(zstep,:) = wRES(zstep,:)'./Mw/(sum(wRES(zstep,:)'./Mw));
end

%% calculate conversion of CO (%)
conversion     = 100*xRES(:,1);
cum_conversion = 100*(x0CO_G - xRES(:,1))./x0CO_G;

%% check the sum of mass fractions
massFractionsIn = sum(wRES(1,:))

dummy = size(wRES);
nRows = dummy(1);

massFractionsOut = sum(wRES(nRows,:))

%% plot the results
figure(1)

subplot(1,2,1)
plot(z,xRES)
xlabel('Reactor length [m]')
ylabel('Gas phase mole fractions (mol/mol)')
legend('CO','H_2','H_2O','C_{16}H_{34}')

subplot(1,2,2)
plot(z,conversion,'k--',z,cum_conversion,'r--')
xlabel('Reactor length [m]')
ylabel('Conversion (%)')
legend('Local Conversion of CO','Cumulative Conversion of CO')




