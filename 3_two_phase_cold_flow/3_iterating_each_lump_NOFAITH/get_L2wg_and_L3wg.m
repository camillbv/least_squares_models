%% camilla.berge.vik@ntnu.no 27.01.2014
%% calculate the L2* operator for gas weight fraction equation

function [L2wg,L3wg] = get_L2wg_and_L3wg(par)

   
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
    L2wg = zeros(N,nCompGas);
    L3wg = zeros(N,nCompGas);
    
    massTransferTerms = zeros(N,nCompGas);
    for zPoint = 1:N
        massTransferTerms(zPoint,:) = kLa.*( 1./equiConst.*...
            avMolMassGas(zPoint)./avMolMassLiq(zPoint).*...
            wtFracGas(zPoint,:)'-wtFracLiq(zPoint,:)' );
    end
    
    for zPoint = 1:N
        L2wg(zPoint,:) = -1./supVelGas(zPoint).*liqDensity./gasDensity(zPoint).*sum(massTransferTerms(zPoint,:));
        L3wg(zPoint,:) =  1./supVelGas(zPoint).*liqDensity./gasDensity(zPoint).*kLa./equiConst.*...
            avMolMassGas(zPoint)./avMolMassLiq(zPoint);
    end

    L2wg = diag(reshape(L2wg,N*nCompGas,1));
    L3wg = diag(reshape(L3wg,N*nCompGas,1));
    
end