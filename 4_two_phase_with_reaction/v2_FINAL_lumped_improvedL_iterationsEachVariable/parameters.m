% camilla.berge.vik@ntnu.no, 28.01.2013

%% purpose: set parameters for the axial dispersion model
% Please see bottom of the file for references.


%% global variables

global a b T nu Mw catDensity volFracSol volFracGas volFracLiq ...
       gasConst supVelGasInit supVelLiqInit pTot factor_a factor_b ...
       kLa equiConst nLumps liqDensityInit

%           value     unit        description         source

%% kinetic parameters
a        = 75.76;  % mmol/(min gcat MPa^2) kin.par. Satterfield(1991), table V
b        = 11.61;  % 1/MPa     kin.par.             Satterfield(1991), table V
alpha    = 0.9;    % [-] parameter in the ASF distribution Krishna(1999), who cites [11] (see article)

%% unit conversion factor of kinetic expression
factor_a = 10^-15/60;
factor_b = 10^-6;

%% physical properties
molWtCO     =  28.0;  % kg/kmol  molar weight CO      SI Chemical Data, page 32
molWtH2     =   2.0;  % kg/kmol  molar weight H2      SI Chemical Data, page 40
molWtH2O    =  18.0;  % kg/kmol  molar weight H2      SI Chemical Data, page 42
Z_CO        = 12.0;   % kg/kmol  atomic weight of carbon
Z_H         = 1.0;    % kg/kmol  atomic weight of hydrogen

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

%% k_i values; k_i = y_i/x_i
K = importdata('utregning_ki_3000_kPa.txt','\t',0); % import text file written from Excel sheet
equiConstAll = K.data;

disp('Manual selection of equiConst values!')
equiConst = equiConstAll([1:3 7 18 26 34]);

%% diffusion coefficients (all valid for 240 deg C)
diffCoefRef = 2*10^-9; % m^2/s reference diffusion coefficient in FT liquid Krishna (1999), page 284
diffCoefCO  = 45.5*10^-9; % m^2/s diffusion coefficient of CO in FT liquid  Krishna (1999), page 284
diffCoefH2  = 17.2*10^-9; % m^2/s diffusion coefficient of H2 in FT liquid  Krishna (1999), page 284

%% operating conditions
pTot            =  30*10^5;  % Pa       total pressure
TCelsius        =  240;  % Celsius   Reactor temperature    Satterfield(1991), table V 
supVelGasInit   =  0.5;  % m/s 	  gas phase superficial velocity Krishna (1999)
supVelLiqInit   =  0.01;  % m/s       liquid phase superficial velocity
volFracSol      =  0.35; % kg/kg     volume fraction of solids in reactor Krishna(2001), page 245 
volFracGas      =  0.3;  % kg/kg     volume fraction of gas in reactor (see figure 4, Krishna 1999)
volFracLiq      =  1-volFracSol - volFracGas; % volume fraction of liquid in reactor

%% mass transfer coefficients
kLaCO =  0.2*volFracGas*0.5*sqrt(diffCoefCO/diffCoefRef); % 1/s volumetric mass transfer coefficient of CO Krishna (1999), page 284 
kLaH2 =  0.2*volFracGas*0.5*sqrt(diffCoefH2/diffCoefRef); % 1/s volumetric mass transfer coefficient of H2 Krishna (1999), page 284 
kLaH2O = 0.05; % set to an arbitrarily value (most of the water is in gas phase at reactor conditions)

% mass transfer coefficients for alkanes
kLaAlkanes = zeros(nLumps,1);
for i = 1:nLumps
    kLaAlkanes(i) = 0.05;
end
kLa = [kLaCO kLaH2 kLaH2O kLaAlkanes']';

%% catalyst properties
catSkeletonDensity   =  3000;  % (kg solids)/(m^3 of particle excluding voids) 
                   
%% densities                   
catDensity  =  catSkeletonDensity*volFracSol; % kg catalyst / m^3   catalyst density in reactor

%% reactor properties
totalHeight 		 =     30; % m 		length of dispersion
T     = TCelsius + 273; % K      Reactor temperature 

%% inlet gas phase concentrations
molFracGasCOInit    = 0.3;   
molFracGasH2Init    = 0.7; 
molFracGasH2OInit   = 0.0;     
molFracGasAlkInit   = zeros(nLumps,1)';
molFracGasInit      = [molFracGasCOInit molFracGasH2Init molFracGasH2OInit molFracGasAlkInit]';
molFracGasInit      = molFracGasInit./sum(molFracGasInit);  % scale it to sum to 1

%% inlet liquid phase concentrations
molFracLiqCOInit    = 0;
molFracLiqH2Init    = 0;
molFracLiqH2OInit   = 0;
molFracLiqAlkInit   = 1e-10*ones(nLumps,1);    
molFracLiqAlkInit(nLumps-1) = 0.5;
molFracLiqAlkInit(nLumps) = 0.5; 
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


%% references
% C. Maretto, R. Krishna, Modelling of a bubble column slurry reactor for
% Fischer–Tropsch synthesis, Catalysis Today, Volume 52, Issues 2–3, 
% 14 September 1999, Pages 279-289

% C Maretto, R Krishna, Design and optimisation of a multi-stage 
% bubble column slurry reactor for Fischer–Tropsch synthesis, 
% Catalysis Today, Volume 66, Issues 2–4, 30 March 2001, Pages 241-248

% I. C. Yates, C. N. Satterfield, Intrinsic kinetics of the 
% Fischer-Tropsch synthesis on a cobalt catalyst, 
% Energy & Fuels 1991 5 (1), pages 168-173

% Aylward, G.H.; T.J.V. Findlay (2002). SI Chemical Data (5th edition ed.).
% Milton, Queensland: John Wiley & Sons Australia

% Hillestad, M.; Modeling the Fischer-Tropsch product distribution and 
% model impelmentation, note, December 16, 2011 (UNPUBLISHED)
