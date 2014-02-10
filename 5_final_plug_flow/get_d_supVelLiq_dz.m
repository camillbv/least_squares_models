%% camilla.berge.vik@ntnu.no 23.01.2014
%% calculate the derivative of the superficial liquid velocity
%% in each spatial point

function d_supVelLiq_dz = get_d_supVelLiq_dz(par)

    %% unpack parameters
    kLa             = par.kLa;
    equiConst       = par.equiConst;
    avMolMassGas    = par.avMolMassGas;
    avMolMassLiq    = par.avMolMassLiq;
    wtFracGas       = par.wtFracGas;
    wtFracLiq       = par.wtFracLiq;

    
    %% derive additional variables
    N = size(wtFracLiq,1);
    nCompLiq = size(wtFracLiq,2);
    
    %% perform calculations
    massTransferTerms = zeros(N,nCompLiq);
    for zPoint = 1:N
        massTransferTerms(zPoint,:) = kLa.*( 1./equiConst.*...
            avMolMassGas(zPoint)./avMolMassLiq(zPoint).*...
            wtFracGas(zPoint,:)'-wtFracLiq(zPoint,:)' );
    end
    
    d_supVelLiq_dz = zeros(N,1); 
    for zPoint = 1:N
        d_supVelLiq_dz(zPoint) = sum(massTransferTerms(zPoint,:));
    end
    
end