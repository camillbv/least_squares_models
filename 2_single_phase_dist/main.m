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
global rho_G Cmax

%% run file to get parameters
parameters

%% set feed gas phase concentrations
x0CO_G     = 0.3;  % mole/mole (unitless) gas feed mole fraction of CO     
x0H2_G     = 0.7;  % mole/mole (unitless) gas feed mole fraction of H2    
x0H2O_G    = 0;    % mole/mole (unitless) gas feed mole fraction of H2O    

x0 = [ x0CO_G  x0H2_G  x0H2O_G  zeros(Cmax+1,1)']';

%% convert to weight fractions
w0   = x0.*Mw./(x0'*Mw); %  kg/kg (unitless)  gas feed weight fractions

%% calculate feed average molar mass
Mw_ave0   = x0'*Mw; % kg/kmol gas feed average molar mass

%% calculate feed gas density
rho_G    = ptot*Mw_ave0/R/T*1000; % kg/m^3  gas feed total density

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

%% calculate conversion of CO (%) on mole basis
xconversion     = 100*xRES(:,1);
xcum_conversion = 100*(x0CO_G - xRES(:,1))./x0CO_G;

%% calculate conversion of CO (%) on mass basis
wconversion     = 100*wRES(:,1);
wcum_conversion = 100*(w0(1) - wRES(:,1))./w0(1);

%% calculate sum of mole fraction of alkanes
for zstep = 1:length(z)
    xalkanes(zstep) = sum(xRES(zstep,4:Cmax+4));
end

%% calculate sum of mass fraction of alkanes
for zstep = 1:length(z)
    walkanes(zstep) = sum(wRES(zstep,4:Cmax+4));
end

%% plot the results: mole fractions
figure(1)

subplot(1,3,1)
plot(z,xRES(:,1:3),z,xalkanes')
xlabel('Reactor length [m]')
ylabel('Gas phase mole fractions (mol/mol)')
legend('CO','H_2','H_2O','C1+')

subplot(1,3,2)
plot(z,xRES(:,4:Cmax+4))
xlabel('Reactor length [m]')
ylabel('Gas phase mole fractions (mol/mol)')
legend('C1','C2','C3',...
    'C4','C5','C6','C7','C8','C9','C10',...
    'C11','C12','C13','C14','C15','C16',...
    'C17','C18','C19','C20','C21','C22',...
    'C23','C24','C25','C26','C27','C28',...
    'C29','C30','C31+')

subplot(1,3,3)
plot(z,xconversion,'k--',z,xcum_conversion,'r--')
xlabel('Reactor length [m]')
ylabel('Conversion (%) on mole basis')
legend('Local Conversion of CO','Cumulative Conversion of CO')

%% plot the results: mass fractions
figure(2)


subplot(1,3,1)
plot(z,wRES(:,1:3),z,walkanes')
xlabel('Reactor length [m]')
ylabel('Gas phase mass fractions (kg/kg)')
legend('CO','H_2','H_2O','C1+')

subplot(1,3,2)
plot(z,wRES(:,4:Cmax+4))
xlabel('Reactor length [m]')
ylabel('Gas phase mass fractions (kg/kg)')
legend('C1','C2','C3',...
    'C4','C5','C6','C7','C8','C9','C10',...
    'C11','C12','C13','C14','C15','C16',...
    'C17','C18','C19','C20','C21','C22',...
    'C23','C24','C25','C26','C27','C28',...
    'C29','C30','C31+')

subplot(1,3,3)
plot(z,wconversion,'k--',z,wcum_conversion,'r--')
xlabel('Reactor length [m]')
ylabel('Conversion (%) on mass basis')
legend('Local Conversion of CO','Cumulative Conversion of CO')
