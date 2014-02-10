%% camilla.berge.vik@ntnu.no 27.01.2014
%% calculate the L2* operator for gas weight fraction equation
%% allowing for calculating selected components only

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
    compIndices     = par.compIndices;

    %% derive additional variables
    N           = size(wtFracGas,1);
    nCompGas    = size(wtFracGas,2);
    nCG         = length(compIndices);
    
    %% perform calculations for all components
    massTransferTerms = zeros(N,nCompGas);
    for zPoint = 1:N
        massTransferTerms(zPoint,:) = kLa.*( 1./equiConst.*...
            avMolMassGas(zPoint)./avMolMassLiq(zPoint).*...
            wtFracGas(zPoint,:)'-wtFracLiq(zPoint,:)' );
    end

    %% perform calculations for selected components
    kLa         = kLa(compIndices);
    equiConst   = equiConst(compIndices);

    L2wg = zeros(N,nCG);
    L3wg = zeros(N,nCG);
    for zPoint = 1:N
        L2wg(zPoint,:) =  -1./supVelGas(zPoint).*liqDensity./gasDensity(zPoint).*sum(massTransferTerms(zPoint,:));
        L3wg(zPoint,:) =  1./supVelGas(zPoint).*liqDensity./gasDensity(zPoint).*kLa./equiConst.*...
            avMolMassGas(zPoint)./avMolMassLiq(zPoint);
    end

    L2wg = diag(reshape(L2wg,N*nCG,1));
    L3wg = diag(reshape(L3wg,N*nCG,1));

end