%% camilla.berge.vik@ntnu.no 23.01.2014
%% calculate the derivative of the gas density
%% in each spatial point

function d_gasDensity_dz = get_d_gasDensity_dz(par)

    %% unpack parameters
    d_avMolMassGas_dz   = par.d_avMolMassGas_dz;
    pTot                = par.pTot;
    gasConst            = par.gasConst;
    T                   = par.T;
    
    %% perform calculations
    d_gasDensity_dz = pTot*d_avMolMassGas_dz./gasConst/T;
    
end