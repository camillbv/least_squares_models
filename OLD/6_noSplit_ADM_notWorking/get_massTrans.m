function massTrans = get_massTrans(par)
 
    %% unpack parameters
    wtFracGas       = par.wtFracGas;
    wtFracLiq       = par.wtFracLiq;
    kLa             = par.kLa;
    equiConst       = par.equiConst;
    avMolMassGas    = par.avMolMassGas;
    avMolMassLiq    = par.avMolMassLiq;
    
    %% derive additional variables
    N = size(wtFracGas,1);
    nCompGas = size(wtFracGas,2);
    
    %% perform calculations   
    massTrans = zeros(N,nCompGas);
    for zPoint = 1:N
        massTrans(zPoint,:) = kLa.*( 1./equiConst.*...
            avMolMassGas(zPoint)./avMolMassLiq(zPoint).*...
            wtFracGas(zPoint,:)'-wtFracLiq(zPoint,:)' );
    end


end