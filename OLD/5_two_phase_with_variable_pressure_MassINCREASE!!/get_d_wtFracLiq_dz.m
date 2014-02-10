%% camilla.berge.vik@ntnu.no 22.01.2014
%% calculate the derivative of gas weight fraction 
%% for all components, in all spatial points

function d_wtFracLiq_dz = get_d_wtFracLiq_dz(par)

    %% unpack parameters
    wtFracGas       = par.wtFracGas;
    wtFracLiq       = par.wtFracLiq;
    supVelLiq       = par.supVelLiq;
    liqDensity      = par.liqDensity;
    kLa             = par.kLa;
    equiConst       = par.equiConst;
    Mw              = par.Mw;
    avMolMassGas    = par.avMolMassGas;
    avMolMassLiq    = par.avMolMassLiq;
    nu              = par.nu;
    reactRate       = par.reactRate;
    catDensity      = par.catDensity;
    
    %% derive additional variables
    N = size(wtFracLiq,1);
    nCompLiq = size(wtFracLiq,2);
    
    %% perform calculations
    d_wtFracLiq_dz = zeros(N,nCompLiq);
    
    massTransferTerms = zeros(N,nCompLiq);
    for zPoint = 1:N
        massTransferTerms(zPoint,:) = kLa.*( 1./equiConst.*...
            avMolMassGas(zPoint)./avMolMassLiq(zPoint).*...
            wtFracGas(zPoint,:)'-wtFracLiq(zPoint,:)' );
    end
  
    for zPoint = 1:N
        d_wtFracLiq_dz(zPoint,:) = -wtFracLiq(zPoint,:)'./supVelLiq(zPoint).*sum(massTransferTerms(zPoint,:)) ...
            + 1./supVelLiq(zPoint).*massTransferTerms(zPoint,:)' ...
            + nu.*reactRate(zPoint)./liqDensity./supVelLiq(zPoint).*catDensity.*Mw;
    end
    
end