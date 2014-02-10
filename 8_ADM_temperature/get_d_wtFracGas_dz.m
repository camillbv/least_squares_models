%% camilla.berge.vik@ntnu.no 22.01.2014
%% calculate the derivative of gas weight fraction 
%% for all components, in all spatial points

function d_wtFracGas_dz = get_d_wtFracGas_dz(par)

   
    %% unpack parameters
    wtFracGas       = par.wtFracGas;
    wtFracLiq       = par.wtFracLiq;
    supVelGas       = par.supVelGas;
    gasDensity      = par.gasDensity;
    liqDensity      = par.liqDensity;
    kLa             = par.kLa;
    equiConst       = par.equiConst;
    avMolMassGas    = par.avMolMassGas;
    avMolMassLiq    = par.avMolMassLiq;
    
    %% derive additional variables
    N = size(wtFracGas,1);
    nCompGas = size(wtFracGas,2);
    
    %% perform calculations
    d_wtFracGas_dz = zeros(N,nCompGas);
    
    massTransferTerms = zeros(N,nCompGas);
    for zPoint = 1:N
        massTransferTerms(zPoint,:) = kLa.*( 1./equiConst.*...
            avMolMassGas(zPoint)./avMolMassLiq(zPoint).*...
            wtFracGas(zPoint,:)'-wtFracLiq(zPoint,:)' );
    end
    
    for zPoint = 1:N
        d_wtFracGas_dz(zPoint,:) = wtFracGas(zPoint,:)'./supVelGas(zPoint).* ...
            liqDensity./gasDensity(zPoint).*sum(massTransferTerms(zPoint,:)) - ...
            liqDensity./gasDensity(zPoint)./supVelGas(zPoint).*massTransferTerms(zPoint,:)';
    end
    
end