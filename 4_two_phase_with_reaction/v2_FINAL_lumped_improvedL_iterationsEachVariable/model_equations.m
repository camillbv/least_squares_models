% camilla.berge.vik@ntnu.no, December 2013

%% purpose: model equations for gas plug flow model

function dydz = model_equations(z,y)

%% declare global variables
global  gasConst T ...
        nCompGas nCompLiq pTot kLa Mw ...
        liqDensityInit nu catDensity equiConst ...
        factor_a factor_b a b volFracLiq

%% unpack 
wtFracGas    = y(1:nCompGas);                 % gas weight fractions
wtFracLiq    = y(nCompGas+1:nCompGas+nCompLiq); % liquid weight fractions
supVelGas    = y(nCompGas+nCompLiq + 1);       % gas superficial velocity
supVelLiq    = y(nCompGas+nCompLiq + 2);       % liquid superficial velocity

%% wtFrac -> molFrac
molFracGas = wtFracGas./Mw/(sum(wtFracGas./Mw)); 
molFracLiq = wtFracLiq./Mw/(sum(wtFracLiq./Mw));

%% average molar mass
avMolMassGas   = molFracGas'*Mw;
avMolMassLiq   = molFracLiq'*Mw;

%% densities
gasDensity = pTot*avMolMassGas/gasConst/T; 
liqDensity = liqDensityInit; % constant 

%% reaction rate 
%reactRate = 0;   % this is a cold flow model 
reactRate = factor_a*(a*equiConst(1)*molFracLiq(1)*equiConst(2)*molFracLiq(2)*pTot^2 )/...
                  ((1+b*factor_b*equiConst(1)*molFracLiq(1)*pTot)^2);


              
%% d_wtFrac_dz    
TRANS = kLa.*(1./equiConst.*avMolMassGas./avMolMassLiq.*wtFracGas-wtFracLiq);

d_wtFracGas_dz = wtFracGas./supVelGas.*liqDensity./gasDensity.*sum(TRANS)...
                 - liqDensity./gasDensity./supVelGas.*TRANS;
d_wtFracLiq_dz = -wtFracLiq./supVelLiq.*sum(TRANS) ...
                 + 1./supVelLiq.*TRANS ...
                 + volFracLiq*nu.*reactRate./liqDensity./supVelLiq.*catDensity.*Mw;
             %pause

%% d_avMolMass_dz
SUMJ                = sum(wtFracGas./Mw);
SUMJd               = sum(d_wtFracGas_dz./Mw);
d_avMolMassGas_dz   = 1./(SUMJ.^2)*sum((d_wtFracGas_dz.*SUMJ - wtFracGas.*SUMJd));


%% d_gasDensity_dz
d_gasDensity_dz = pTot*d_avMolMassGas_dz/gasConst/T;

%% d_supVel_dz
d_supVelGas_dz = -supVelGas./gasDensity.*d_gasDensity_dz - liqDensity./gasDensity.*sum(TRANS);
d_supVelLiq_dz = sum(TRANS);

%% diagnostics: verify mass conservation
% disp('Mass conservation: The following two numbers should be one:')
% [sum(wtFracGas) sum(wtFracLiq)]

%% assemble differentials
dydz = [d_wtFracGas_dz; d_wtFracLiq_dz; d_supVelGas_dz; d_supVelLiq_dz];
%pause
end