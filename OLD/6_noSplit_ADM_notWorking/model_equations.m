% camilla.berge.vik@ntnu.no, December 2013

%% purpose: model equations for gas plug flow model

function dydz = model_equations(z,y)

%% declare global variables
global  gasConst T nCompGas nCompLiq kLa Mw ...
        liqDensityInit  nu catDensity equiConst ...
        pTotInit supVelGasInit supVelLiqInit volFracLiqInit

%% unpack 
wtFracGas    = y(1:nCompGas);                 % gas weight fractions
wtFracLiq    = y(nCompGas+1:nCompGas+nCompLiq); % liquid weight fractions

%% wtFrac -> molFrac

molFracGas = wtFracGas./Mw/(sum(wtFracGas./Mw)); 
molFracLiq = wtFracLiq./Mw/(sum(wtFracLiq./Mw));

%% average molar mass
avMolMassGas   = molFracGas'*Mw;
avMolMassLiq   = molFracLiq'*Mw;

%% densities
gasDensity = pTotInit*avMolMassGas/gasConst/T; 
liqDensity = liqDensityInit; % constant 

%% reaction rate 
reactRate = 0;
              
%% d_wtFrac_dz    
TRANS = kLa.*(1./equiConst.*avMolMassGas./avMolMassLiq.*wtFracGas-wtFracLiq);

d_wtFracGas_dz = wtFracGas./supVelGasInit.*liqDensity./gasDensity.*sum(TRANS)...
                 - liqDensity./gasDensity./supVelGasInit.*TRANS;
d_wtFracLiq_dz = -wtFracLiq./supVelLiqInit.*sum(TRANS) ...
                 + 1./supVelLiqInit.*TRANS ...
                 + volFracLiqInit*nu.*reactRate./liqDensity./supVelLiqInit.*catDensity.*Mw;

%% assemble differentials
dydz = [d_wtFracGas_dz; d_wtFracLiq_dz]; 
    
end