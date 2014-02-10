% camilla.berge.vik@ntnu.no, 19.09.2013

%% purpose: set parameters for the axial dispersion model
% Please see bottom of the file for references.


%% global variables

global a b T nu Mw catDensity volFracSol volFracGas volFracLiq ...
       gasConst supVelGasInit supVelLiqInit pTot factor_a factor_b Cmax ...
       kLa equiConst

%           value     unit        description         source

%% kinetic parameters
a        = 75.76;  % mmol/(min gcat MPa^2) kin.par. Satterfield(1991), table V
b        = 11.61;  % 1/MPa     kin.par.             Satterfield(1991), table V
alpha    = 0.9;    % [-] parameter in the ASF distribution Krishna(1999), who cites [11] (see article)

%% unit conversion factor of kinetic expression
factor_a = 10^-15/60;
factor_b = 10^-6;

%% stochiometric coefficients

nuCO     = -1   ;  % mol/mol  stochio. coeff. CO   
nuH2     = -(3-alpha);  % mol/mol  stochio. coeff. H2
nuH2O    =  1;  % mol/mol  stochio. coeff. H2O

Cmax = 30;  % highest carbon number not lumped

% stochiometric coefficients for non-lumped alkanes
nuAlkanes = zeros(Cmax,1);
for i = 1:Cmax
    nuAlkanes(i)     = (1-alpha)^2*alpha^(i-1);
end

% stochiometric coefficient for lumped alkanes
nuAlkanes(Cmax + 1) = (1-alpha)*alpha^(Cmax + 1 - 1);

nu = [nuCO nuH2 nuH2O nuAlkanes']';

%% universal gas constant
gasConst        = 8.314*10^3;  % J/(K kmol) gas constant         SI Chemical Data, page 3

%% physical properties
molWtCO     =  28.0;  % kg/kmol  molar weight CO      SI Chemical Data, page 32
molWtH2     =   2.0;  % kg/kmol  molar weight H2      SI Chemical Data, page 40
molWtH2O    =  18.0;  % kg/kmol  molar weight H2      SI Chemical Data, page 42
Z_CO     = 12.0;   % kg/kmol  atomic weight of carbon
Z_H      = 1.0;    % kg/kmol  atomic weight of hydrogen

%% molar weights for non-lumped alkanes
molWtAlkanes = zeros(Cmax+1,1);
molWtAlkanes(1) = 16; % kg/kmol molar weight CH4    SI Chemical Data, page 90
for i = 2:Cmax % excluding methane
    molWtAlkanes(i) = i*Z_CO + (2*i+2)*Z_H;
end

%% average molar weight for lumped alkanes
nAverage = (Cmax + 1) + alpha/(1-alpha);          % Hillestad (2011)
molWtAlkanes(Cmax + 1) = nAverage*Z_CO + (2*nAverage +2)*Z_H;

Mw  = [molWtCO   molWtH2   molWtH2O   molWtAlkanes']'; % kg/kmol  col. vector of molar masses

%%k_i values; k_i = y_i/x_i
K = importdata('utregning_ki_3000_kPa.txt','\t',0); % import text file written from Excel sheet
equiConst = K.data;
% for the case when the constant is equal to 0, set to an arbitrarily low
% value


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
kLaCO =  0.1*volFracGas*0.5*sqrt(diffCoefCO/diffCoefRef); % 1/s volumetric mass transfer coefficient of CO Krishna (1999), page 284 
kLaH2 =  0.1*volFracGas*0.5*sqrt(diffCoefH2/diffCoefRef); % 1/s volumetric mass transfer coefficient of H2 Krishna (1999), page 284 
kLaH2O = 0.005; % set to an arbitrarily value (most of the water is in gas phase at reactor conditions)

% mass transfer coefficients for alkanes
kLaAlkanes = zeros(Cmax+1,1);
for i = 1:Cmax+1
    kLaAlkanes(i) = 0.0005;
end
kLa = [kLaCO kLaH2 kLaH2O kLaAlkanes']';

%% catalyst properties
catSkeletonDensity   =  3000;  % (kg solids)/(m^3 of particle excluding voids) 
                   
%% densities                   
catDensity  =  catSkeletonDensity*volFracSol; % kg catalyst / m^3   catalyst density in reactor


%% reactor properties
totalHeight 		 =     30; % m 		length of dispersion
T     = TCelsius + 273; % K      Reactor temperature 

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
