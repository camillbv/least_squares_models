% camilla.berge.vik@ntnu.no, 23.09.2013

%% purpose: model equations for gas plug flow model

function dydz = model_equations(z,y)

%% declare global variables
global  gasConst ...
        nCompGas nCompLiq Mw ...
        liqDensityInit nu catDensity equiConst ...
        heatCapGas heatCapSluInit perimeter area ...
        hL areaDensity sluDensityInit tempSurr ...
        deltaHr pTot graConst condLiq ...
        condSol volFracSol volFracLiq viscLiq einsteinK ...
        molecularDia mWoctacosane mDoctacosane Nparameter 

%% unpack 
wtFracGas    = y(1:nCompGas);                   % gas weight fractions
wtFracLiq    = y(nCompGas+1:nCompGas+nCompLiq); % liquid weight fractions
supVelGas    = y(nCompGas+nCompLiq + 1);        % gas superficial velocity
supVelLiq    = y(nCompGas+nCompLiq + 2);        % liquid superficial velocity
tempGas      = y(nCompGas+nCompLiq + 3);        % gas temperature
tempSlu      = y(nCompGas+nCompLiq + 4);        % slurry temperature

supVelSlu=supVelLiq;
sluDensity = sluDensityInit;
heatCapSlu = heatCapSluInit;

%% wtFrac -> molFrac
molFracGas = wtFracGas./Mw/(sum(wtFracGas./Mw)); 
molFracLiq = wtFracLiq./Mw/(sum(wtFracLiq./Mw));

%% average molar mass
avMolMassGas   = molFracGas'*Mw;
avMolMassLiq   = molFracLiq'*Mw;

%% densities
gasDensity = pTot*avMolMassGas/gasConst/tempGas; 
liqDensity = liqDensityInit; % constant 

%% reaction rate in liquid (slurry) phase
aS = 2.592e-12*exp(-37.3*1000/(gasConst/1000)./tempSlu);
bS =  1.23e-12*exp( 68.5*1000/(gasConst/1000)./tempSlu);

reactRate = max((aS*equiConst(1)*molFracLiq(1)*equiConst(2)*molFracLiq(2)*pTot^2 )/...
                  ((1+bS*equiConst(1)*molFracLiq(1)*pTot)^2),0);
              
%% mass transfer between gas and slurry phase  
parStructkL = struct('liqDensity',liqDensity,'gasDensity',gasDensity, ...
    'viscLiq',viscLiq,'graConst',graConst, ...
    'molecularDia',molecularDia,'Mw',Mw,'mWoctacosane',mWoctacosane, ...
    'mDoctacosane',mDoctacosane,'Nparameter',Nparameter, ...
    'tempSlu',tempSlu);
kL = getMassTransCoeff(parStructkL)';

massTrans = kL.*areaDensity.*(1./equiConst.*avMolMassGas./avMolMassLiq.*wtFracGas-wtFracLiq);

%% d_wtFrac_dz
d_wtFracGas_dz = wtFracGas./supVelGas.*liqDensity./gasDensity.*sum(massTrans)...
                 - liqDensity./gasDensity./supVelGas.*massTrans;
d_wtFracLiq_dz = -wtFracLiq./supVelLiq.*sum(massTrans) ...
                 + 1./supVelLiq.*massTrans ...
                 + volFracLiq.*nu.*reactRate./liqDensity./supVelLiq.*catDensity.*Mw;

%% d_avMolMass_dz
SUMJ                = sum(wtFracGas./Mw);
SUMJd               = sum(d_wtFracGas_dz./Mw);
d_avMolMassGas_dz   = 1./(SUMJ.^2)*sum((d_wtFracGas_dz.*SUMJ - wtFracGas.*SUMJd));

%% d_temp_dz
parStructhW = struct('condLiq',condLiq,'condSol',condSol, ...
    'volFracSol',volFracSol,'viscLiq',viscLiq,'graConst',graConst, ...
    'supVelGas',supVelGas,'sluDensity',sluDensity, ...
    'einsteinK',einsteinK,'heatCapSlu',heatCapSlu);

hW = getWallHeatCoeff(parStructhW);

heatRemovalThroughCoolingTubes = hW*perimeter/area/(sluDensityInit*heatCapSluInit*supVelSlu)*(tempSurr-tempSlu);
heatTransLiqToGas = hL*areaDensity/(sluDensityInit*heatCapSluInit*supVelSlu)*(tempGas-tempSlu);

reactionHeat = volFracLiq/(sluDensityInit*heatCapSluInit*supVelSlu)*(-deltaHr)*reactRate*catDensity;

d_tempGas_dz = -hL*areaDensity/(gasDensity*heatCapGas*supVelGas)*(tempGas-tempSlu);
d_tempSlu_dz = reactionHeat + heatRemovalThroughCoolingTubes + heatTransLiqToGas;

%% d_gasDensity_dz
d_gasDensity_dz = pTot*d_avMolMassGas_dz/gasConst/tempGas - ...
                    avMolMassGas*pTot/gasConst.*d_tempGas_dz./tempGas.^2;

%% d_supVel_dz
d_supVelGas_dz = -supVelGas./gasDensity.*d_gasDensity_dz - liqDensity./gasDensity.*sum(massTrans);
d_supVelLiq_dz = sum(massTrans);
           
%% verify mass conservation
%disp('Mass conservation: The following two numbers should be one:')
%[sum(wtFracGas) sum(wtFracLiq)]

%% assemble differentials
dydz = [d_wtFracGas_dz; d_wtFracLiq_dz; ...
        d_supVelGas_dz; d_supVelLiq_dz; ...
        d_tempGas_dz; d_tempSlu_dz];
end