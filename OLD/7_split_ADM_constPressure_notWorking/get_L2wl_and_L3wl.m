%% camilla.berge.vik@ntnu.no 27.01.2014
%% calculate the L2* and L3* operator for liquid weight fraction equation
%% allowing for calculating selected components only

function [L2wl,L3wl] = get_L2wl_and_L3wl(par)

   
    %% unpack parameters
    wtFracGas       = par.wtFracGas;
    wtFracLiq       = par.wtFracLiq;
    supVelLiq       = par.supVelLiq;
    kLa             = par.kLa;
    equiConst       = par.equiConst;
    avMolMassGas    = par.avMolMassGas;
    avMolMassLiq    = par.avMolMassLiq;
    compIndices     = par.compIndices;
    
    %% derive additional variables
    N           = size(wtFracLiq,1);
    nCompLiq    = size(wtFracLiq,2);

    %% perform calculations for all components
    massTransferTerms = zeros(N,nCompLiq);
    for zPoint = 1:N
        massTransferTerms(zPoint,:) = kLa.*( 1./equiConst.*...
            avMolMassGas(zPoint)./avMolMassLiq(zPoint).*...
            wtFracGas(zPoint,:)'-wtFracLiq(zPoint,:)' );
    end
    
    %% perform calculations for selected components    
    nCL         = length(compIndices);
    kLa         = kLa(compIndices);
    L2wl        = zeros(N,nCL);
    L3wl        = zeros(N,nCL);
    for zPoint = 1:N
        L2wl(zPoint,:) = 1./supVelLiq(zPoint).*sum(massTransferTerms(zPoint,:));
        L3wl(zPoint,:) = 1./supVelLiq(zPoint).*kLa;
    end

    L2wl = diag(reshape(L2wl,N*nCL,1));
    L3wl = diag(reshape(L3wl,N*nCL,1));
    
end