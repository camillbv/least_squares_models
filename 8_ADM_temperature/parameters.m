% camilla.berge.vik@ntnu.no, 07.02.2014

%% purpose: set parameters for the axial dispersion model
% Please see bottom of the file for references.


%% global variables

global a b tempGasInit tempSluInit nu Mw catDensity ...
       volFracSol volFracGas volFracLiq hL heatCapGas perimeter...
       gasConst supVelGasInit supVelLiqInit pTot factor_a factor_b ...
       kLa equiConst nLumps liqDensityInit dispCoefGas dispCoefLiq ...
       nComp deltaHr areaDensity heatCapSluInit sluDensityInit area ...
       volFracSlu tempSurr graConst ...
       condLiq condSol viscLiq einsteinK ...
       molecularDia mWoctacosane mDoctacosane Nparameter 

%           value     unit        description         source
%% kinetic parameters
a        = 53.11;  % mmol/(min gcat MPa^2) kin.par. Satterfield(1991), table V
b        = 22.26;  % 1/MPa     kin.par.             Satterfield(1991), table V
alpha    = 0.9;    % [-] parameter in the ASF distribution Krishna(1999), who cites [11] (see article)
deltaHr  = -167*10^6; % J/kmol  heat of reaction Visconti2013

%% unit conversion factor of kinetic expression
factor_a = 10^-15/60;
factor_b = 10^-6;

%% physical properties
molWtCO     =  28.0;  % kg/kmol  molar weight CO      SI Chemical Data, page 32
molWtH2     =   2.0;  % kg/kmol  molar weight H2      SI Chemical Data, page 40
molWtH2O    =  18.0;  % kg/kmol  molar weight H2O     SI Chemical Data, page 42
Z_CO        =  12.0;  % kg/kmol  atomic weight of carbon
Z_H         =   1.0;  % kg/kmol  atomic weight of hydrogen

%% stochiometric coefficients and molar weights
% non-lumped components
nuCO     = -1   ;  % mol/mol  stochio. coeff. CO   
nuH2     = -(3-alpha);  % mol/mol  stochio. coeff. H2
nuH2O    =  1;  % mol/mol  stochio. coeff. H2O

lumpKey     = [1 11 21 31]';
nLumps      = length(lumpKey);

% lumped components (alkanes)
nuAlkanes = zeros(nLumps,1);
avCNoAlkanes = zeros(nLumps,1);
molWtAlkanes = zeros(nLumps,1);
for lNo = 1:(nLumps-1)
    n = lumpKey(lNo);
    m = lumpKey(lNo+1);
    nuAlkanes(lNo) = (1 - alpha)*(alpha^(n-1)-alpha^m);
    avCNoAlkanes(lNo) = (n*alpha^(n-1) - (n-1)*alpha^n - (m+1)*alpha^m + m*(alpha^(m+1))) / ...
                            ((1-alpha)*(alpha^(n-1)-alpha^m));
    molWtAlkanes(lNo) = avCNoAlkanes(lNo)*Z_CO + (2*avCNoAlkanes(lNo)+2)*Z_H;
end
n = lumpKey(nLumps);
nuAlkanes(nLumps) = (1-alpha)*(alpha^(n-1));
avCNoAlkanes(nLumps) = n + alpha/(1-alpha);
molWtAlkanes(nLumps) = avCNoAlkanes(nLumps)*Z_CO + (2*avCNoAlkanes(nLumps)+2)*Z_H;

% assemble non-lumped and lumped components
nu  = [nuCO nuH2 nuH2O nuAlkanes']';
Mw  = [molWtCO   molWtH2   molWtH2O   molWtAlkanes']'; % kg/kmol  col. vector of molar masses

%% universal gas constant
gasConst        = 8.314*10^3;  % J/(K kmol) gas constant         SI Chemical Data, page 3
graConst        = 9.81;        % m/s^2      standard acceleration of gravity

%% k_i values; k_i = y_i/x_i
K = importdata('utregning_ki_3000_kPa.txt','\t',0); % import text file written from Excel sheet
equiConstAll = K.data;

disp('Manual selection of equiConst values!')
equiConst = equiConstAll([1:3 7 18 26 34]);

%% diffusion coefficients (all valid for 240 deg C)
diffCoefRef = 2*10^-9; % m^2/s reference diffusion coefficient in FT liquid Krishna (1999), page 284
diffCoefCO  = 45.5*10^-9; % m^2/s diffusion coefficient of CO in FT liquid  Krishna (1999), page 284
diffCoefH2  = 17.2*10^-9; % m^2/s diffusion coefficient of H2 in FT liquid  Krishna (1999), page 284

%% diffusion / mass transfer parameters
% components: [CO H2 H2O C1-C10 C11-C20 C21-C30 C31+]
% repr. comp: [CO H2 H2O C5     C14     C24     C28*] (This is the component
% on which the liquid properties are estimated)
molecularDia = [3.72 2.92 2.7 5.850 7.81 9.25 9.76]'; % Å
Nparameter   = 0.605; % backwards-calculation from eq 7, ErkeyRoddenAkgerman1990, CO in octacosane
% sources: CO, H2, C5: Bird, Stewart, Lightfoot (2nd ed), page 854 
%          C14, C24, C28: Erkey, Rodden, Akgerman (1990), table 2. 
%          Data for C24 extrapolated between C20 and C28.
%          H2O: Schatzberg (1967)
mWoctacosane = 394.8; % molecular weight of octacosane
mDoctacosane = 9.76;  % molecular diameter of octacosane

%% operating conditions
pTot            =  30*10^5;  % Pa       total pressure
tempCelsiusInit =  240;  % Celsius   Reactor temperature    Satterfield(1991), table V 
supVelGasInit   =  0.2;  % m/s 	  gas phase superficial velocity Krishna (1999)
supVelLiqInit   =  0.01;  % m/s       liquid phase superficial velocity
volFracSol      =  0.3; % kg/kg     volume fraction of solids in reactor Krishna(2001), page 245 
volFracGas      =  0.2;  % kg/kg     volume fraction of gas in reactor (see figure 4, Krishna 1999)
volFracLiq      =  1-volFracSol - volFracGas; % volume fraction of liquid in reactor
volFracSlu      = volFracLiq + volFracGas; % volume fraction of slurry in reactor

%% area density
sauterDiameter  = 10*10^-3; % m Sauter Mean Diameter ADJUSTED
areaDensity     = 6*volFracGas/sauterDiameter;

%% dispersion coefficients (1 m column diameter)
diaCol          = 8;   % m GUESSING
dispCoefLiq     = 0.68*diaCol^1.4*supVelGasInit^0.3;
dispCoefGas     = 21.7*diaCol^1.5*supVelGasInit^1.8;  

%% thermal conductivities
einsteinK = 4.5;        % Meisl (1967)
viscLiq   = 1.05*1e-3;  % Pa*s , Sehabiague 2013a, table 3, Heavy F-T Cut at 530 K
%viscLiq   = 2.9*1e-4;  % Pa*s , Sehabiague 2013a, table 3, Heavy F-T Cut at 530 K

condLiq   = 0.113;      % W/(mK), Maretto Krishna 1999, table 1
condSol   =   1.7;      % W/(mK), Maretto Krishna 1999, table 1

% %% mass transfer coefficients
% kLaCO =  volFracGas*0.5*sqrt(diffCoefCO/diffCoefRef); % 1/s volumetric mass transfer coefficient of CO Krishna (1999), page 284 
% kLaH2 =  volFracGas*0.5*sqrt(diffCoefH2/diffCoefRef); % 1/s volumetric mass transfer coefficient of H2 Krishna (1999), page 284 
% kLaH2O = 0.25; % set to an arbitrarily value (most of the water is in gas phase at reactor conditions)
% 
% % mass transfer coefficients for alkanes
% kLaAlkanes = zeros(nLumps,1);
% for i = 1:nLumps
%     kLaAlkanes(i) = 0.1;
% end
% kLa = [kLaCO kLaH2 kLaH2O kLaAlkanes']';

kLa =    [0.066   0.105   0.0852   0.0450   0.0328   0.027   0.025]'; % simualted

%% catalyst properties
catSkeletonDensity   =  2030;  % (kg solids)/(m^3 of particle excluding voids) 
catParticleDensity   =   647;  % (kg solids and liquids)/(m^3 of particle including voids) MarettoKrishna1999
                   
%% densities                   
catDensity  =  catSkeletonDensity*volFracSol; % kg catalyst / m^3   catalyst density in reactor

%% reactor properties
totalHeight 		 =     50; % m 		length of dispersion

%% temperature
tempGasInit     = tempCelsiusInit + 273; % K      Reactor temperature 
tempSluInit     = tempGasInit;

%% temperature equations parameters
heatCapLiq = 1500;  % J/(kg K)% heat capacity of liquid     MarettoKrishna1999
heatCapSol = 992;   % J/(kg K)% heat capacity of solid      MarettoKrishna1999
heatCapGas = 3000;  % J/(kg K)% heat capacity of methane (TO FIX!)

%% heat transfer parameters
hL =   6; % GUESSWORK

%% reactor cooling system
diaTub          = 0.114;    % m outer diameter of cooling tube
numTub          = 1200;     % - number of cooling tubes
perimeter       = diaCol   + numTub*diaTub;
area            = 4*(diaCol^2 - numTub*diaTub^2);
tempSurr        = 350;      % K cooling water temperature

%% inlet gas phase concentrations
molFracGasCOInit    = 0.3;   
molFracGasH2Init    = 0.6; 
molFracGasH2OInit   = 0.00001;     
molFracGasAlkInit   = zeros(nLumps,1)';
molFracGasInit      = [molFracGasCOInit molFracGasH2Init molFracGasH2OInit molFracGasAlkInit]';
molFracGasInit      = molFracGasInit./sum(molFracGasInit);  % scale it to sum to 1

%% inlet liquid phase concentrations
molFracLiqCOInit    = 0.01;
molFracLiqH2Init    = 0.01;
molFracLiqH2OInit   = 0;
molFracLiqAlkInit   = 1e-10*ones(nLumps,1); 
molFracLiqAlkInit(nLumps) = 1; 
molFracLiqInit      = [molFracLiqCOInit molFracLiqH2Init molFracLiqH2OInit molFracLiqAlkInit']';     % add water
molFracLiqInit      = molFracLiqInit./sum(molFracLiqInit);  % scale it to sum to 1

%% convert to weight fractions
wtFracGasInit   = molFracGasInit.*Mw./(molFracGasInit'*Mw); %  kg/kg (unitless)  gas feed weight fractions
wtFracLiqInit   = molFracLiqInit.*Mw./(molFracLiqInit'*Mw); %  kg/kg (unitless)  gas feed weight fractions

%% number of components
nCompGas = length(wtFracGasInit);
nCompLiq = length(wtFracLiqInit);
nComp    = nCompGas; 

%% feed average molar mass
avMolMassGasInit   = molFracGasInit'*Mw; % kg/kmol gas feed average molar mass

%% feed gas density
gasDensityInit    = pTot*avMolMassGasInit/gasConst/tempGasInit; % kg/m^3  gas feed total density
liqDensityInit    = 651; % ErkeyRoddenAkgerman (1988), table V, 534 K

sluDensityInit    = liqDensityInit*(1-liqDensityInit/catSkeletonDensity*volFracSol) + ...
                    catParticleDensity*volFracSol; %MarettoKrishna1999, equation 23

totDensityInit    = volFracGas*gasDensityInit + volFracSlu*sluDensityInit;


wtFracSolSluInit  = volFracSol*catParticleDensity/sluDensityInit;
heatCapSluInit    = wtFracSolSluInit*heatCapSol + (1-wtFracSolSluInit)*heatCapLiq;% heat capacity of slurry             
     


