%% camilla.berge.vik@ntnu.no 23.01.2014
%% calculate the derivative of the superficial gas velocity
%% in each spatial point

function d_supVelGas_dz = get_d_supVelGas_dz(par)

    %% unpack parameters
    supVelGas       = par.supVelGas;
    gasDensity      = par.gasDensity;
    liqDensity      = par.liqDensity;
    d_gasDensity_dz = par.d_gasDensity_dz;
    kLa             = par.kLa;
    equiConst       = par.equiConst;
    avMolMassGas    = par.avMolMassGas;
    avMolMassLiq    = par.avMolMassLiq;
    wtFracGas       = par.wtFracGas;
    wtFracLiq       = par.wtFracLiq;
    
    %% derive additional variables
    N = size(wtFracGas,1);
    nCompGas = size(wtFracGas,2);
    
    %% perform calculations
    massTransferTerms = zeros(N,nCompGas);
    for zPoint = 1:N
        massTransferTerms(zPoint,:) = kLa.*( 1./equiConst.*...
            avMolMassGas(zPoint)./avMolMassLiq(zPoint).*...
            wtFracGas(zPoint,:)'-wtFracLiq(zPoint,:)' );
    end
    
    d_supVelGas_dz = zeros(N,1);
    for zPoint = 1:N
        d_supVelGas_dz(zPoint) = -supVelGas(zPoint)./gasDensity(zPoint).*d_gasDensity_dz(zPoint) - ...
            liqDensity./gasDensity(zPoint).*sum(massTransferTerms(zPoint,:));
    end
    
    
end