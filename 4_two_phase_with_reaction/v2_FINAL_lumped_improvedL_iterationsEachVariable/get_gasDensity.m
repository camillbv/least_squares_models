%% camilla.berge.vik@ntnu.no 22.01.2014
%% calculate gas density in each spatial point

function gasDensity = get_gasDensity(par)

    %% unpack parameters
    avMolMassGas = par.avMolMassGas;
    pTot         = par.pTot;
    gasConst     = par.gasConst;
    T            = par.T;
    
    %% calculate gas density
    gasDensity = pTot*avMolMassGas/gasConst/T;

end