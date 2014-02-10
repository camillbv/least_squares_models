% camilla.berge.vik@ntnu.no, 19.09.2013

%% purpose: set parameters
% Please see bottom of the file for references.

%% global variables

global a b T nu Mw rho_cat eps_S eps_G R v_GS ptot factor factor_a factor_b

%           value     unit        description         source

%% kinetic parameters
a        = 75.76;  % mmol/(min gcat Mpa^2) kin.par. Satterfield(1991), table V
b        = 11.61;  % 1/MPa     kin.par.             Satterfield(1991), table V

%% unit conversion factor of kinetic expression
factor   = 10^-3/60;
factor_a = 10^-12;
factor_b = 10^-6;

%% stochiometric coefficients using C16H34 as the only product
% CO + uH2  ->  1/16 C16H34 + H2O
% ---> u = (34/16 + 1) = 2.0625

nuCO     = -1   ;  % mol/mol  stochio. coeff. CO   
nuH2     = -2.0625;  % mol/mol  stochio. coeff. H2
nuH2O    =  1;  % mol/mol  stochio. coeff. H2O
nuC16H34 =  1/16   ;  % mol/mol  stochio. coeff. C16H34
nu  = [nuCO   nuH2   nuH2O   nuC16H34  ]'; %  -       col. vector of stoich. coeff.

%% universal constants
R        = 8.314*10^3;  % J/(K mol) gas constant         SI Chemical Data, page 3

%% physical properties
MwCO     =  28.0;  % kg/kmol  molar weight CO      SI Chemical Data, page 32
MwH2     =   2.0;  % kg/kmol  molar weight H2      SI Chemical Data, page 40
MwH2O    =  18.0;  % kg/kmol  molar weight H2      SI Chemical Data, page 42
MwC16H34 = 226.0;  % kg/kmol  molar weight C16H34  SI Chemical Data, page 90
Mw  = [MwCO   MwH2   MwH2O   MwC16H34  ]'; % kg/kmol  col. vector of molar masses

%% operating conditions
ptot     =  30*10^5;  % Pa       total pressure
T_C      =  240;  % Celsius   Reactor temperature    Satterfield(1991), table V 
v_GS     =  1;    % m/s 	  gas phase superficial velocity Krishna (1999)
eps_S    =  0.35; % kg/kg     mass fraction of solids in reactor Krishna(2001), page 245 
eps_G    =  0.3;  % kg/kg     mass fraction of gas in reactor (see figure 4, Krishna 1999)

%% catalyst
rho_SK   =  2030;  % (kg solids)/(m^3 of particle excluding voids) 
                   % catalyst skeleton density      Krishna(1999), table 1
rho_cat  =  rho_SK*eps_S; % kg catalyst / m^3   catalyst density in reactor

%% reactor properties
L 		 =     30; % m 		length of dispersion (Krishna1999)
T     = T_C + 273; % K      Reactor temperature 
