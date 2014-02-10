%% camilla.berge.vik@ntnu.no 27.01.2014
%% calculate the L2* and L3* operator for liquid weight fraction equation

function [L2wl,L3wl] = get_L2wl_and_L3wl(par)

   
    %% unpack parameters
    wtFracGas       = par.wtFracGas;
    wtFracLiq       = par.wtFracLiq;
    supVelLiq       = par.supVelLiq;
    kLa             = par.kLa;
    equiConst       = par.equiConst;
    avMolMassGas    = par.avMolMassGas;
    avMolMassLiq    = par.avMolMassLiq;
    
    %% derive additional variables
    N = size(wtFracLiq,1);
    nCompLiq = size(wtFracLiq,2);
    
    %% perform calculations
    L2wl = zeros(N,nCompLiq);
    L3wl = zeros(N,nCompLiq);
    
    massTransferTerms = zeros(N,nCompLiq);
    for zPoint = 1:N
        massTransferTerms(zPoint,:) = kLa.*( 1./equiConst.*...
            avMolMassGas(zPoint)./avMolMassLiq(zPoint).*...
            wtFracGas(zPoint,:)'-wtFracLiq(zPoint,:)' );
    end
    
    for zPoint = 1:N
        L2wl(zPoint,:) = -1./supVelLiq(zPoint).*sum(massTransferTerms(zPoint,:));
        L3wl(zPoint,:) = 1./supVelLiq(zPoint).*kLa;
    end

    L2wl = diag(reshape(L2wl,N*nCompLiq,1));
    L3wl = diag(reshape(L3wl,N*nCompLiq,1));
    
end