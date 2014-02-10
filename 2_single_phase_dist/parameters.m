% camilla.berge.vik@ntnu.no, 19.09.2013
% revised: 21.01.2014

%% purpose: set parameters for the axial dispersion model
% Please see bottom of the file for references.


%% global variables

global a b T nu Mw rho_cat eps_S eps_G R v_GS ptot factor Cmax factor_a factor_b

%           value     unit        description         source

%% kinetic parameters
a        = 75.76;  % mmol/(min gcat Mpa^2) kin.par. Satterfield(1991), table V
b        = 11.61;  % 1/MPa     kin.par.             Satterfield(1991), table V
alpha    = 0.9;    % [-] parameter in the ASF distribution Krishna(1999), who cites [11] (see article)

%% unit conversion factor of kinetic expression
factor = 10^-3/60;
factor_a = 10^-12;
factor_b = 10^-6;
%% stochiometric coefficients using C16H34 as the only product

nuCO     = -1   ;  % mol/mol  stochio. coeff. CO   
nuH2     = -(3-alpha);  % mol/mol  stochio. coeff. H2
nuH2O    =  1;  % mol/mol  stochio. coeff. H2O

Cmax = 30;  % highest carbon number not lumped

%% stochiometric coefficients for non-lumped alkanes
nu_alkanes = zeros(Cmax,1);
for i = 1:30
    nu_alkanes(i)     = (1-alpha)^2*alpha^(i-1);
end

%% stochiometric coefficient for lumped alkanes
nu_alkanes(Cmax + 1) = (1-alpha)*alpha^(Cmax + 1 - 1);

nu = [nuCO nuH2 nuH2O nu_alkanes']';

%% universal constants
R        = 8.314*10^3;  % J/(K mol) gas constant         SI Chemical Data, page 3

%% physical properties
MwCO     =  28.0;  % kg/kmol  molar weight CO      SI Chemical Data, page 32
MwH2     =   2.0;  % kg/kmol  molar weight H2      SI Chemical Data, page 40
MwH2O    =  18.0;  % kg/kmol  molar weight H2      SI Chemical Data, page 42
Z_CO     = 12.0;   % kg/kmol  atomic weight of carbon
Z_H      = 1.0;    % kg/kmol  atomic weight of hydrogen

%% molar weights for non-lumped alkanes
Mw_alkanes = zeros(Cmax,1);
Mw_alkanes(1) = 16; % kg/kmol molar weight CH4    SI Chemical Data, page 90
for i = 2:30
    Mw_alkanes(i) = i*Z_CO + (2*i+2)*Z_H;
end

%% average molar weight for lumped alkanes
n_average = (Cmax + 1) + alpha/(1-alpha);          % Hillestad (2011)
Mw_alkanes(Cmax + 1) = n_average*Z_CO + (2*n_average +2)*Z_H;

Mw  = [MwCO   MwH2   MwH2O   Mw_alkanes']'; % kg/kmol  col. vector of molar masses

%% operating conditions
ptot     =  30*10^5; % Pa       total pressure
T_C      =  240;  % Celsius   Reactor temperature    Satterfield(1991), table V 
v_GS     =  1;  % m/s 	  gas phase superficial velocity Krishna (1999)
eps_S    =  0.35; % kg/kg     mass fraction of solids in reactor Krishna(2001), page 245 
eps_G    =  0.3;  % kg/kg     mass fraction of gas in reactor (see figure 4, Krishna 1999)

%% catalyst
rho_SK   =  2030;  % (kg solids)/(m^3 of particle excluding voids) 
                   % catalyst skeleton density      Krishna(1999), table 1
rho_cat  =  rho_SK*eps_S; % kg catalyst / m^3   catalyst density in reactor

%% reactor properties
L 		 =     30; % m 		length of dispersion (Krishna1999)
T     = T_C + 273; % K      Reactor temperature 

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
