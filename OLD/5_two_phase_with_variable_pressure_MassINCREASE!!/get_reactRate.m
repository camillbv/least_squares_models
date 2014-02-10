%% camilla.berge.vik@ntnu.no 24.01.2014
%% calculate the reaction rate in each spatial point
function reactRate = get_reactRate(par)
    %% unpack variables
    volFracLiq  = par.volFracLiq;
    factor_a    = par.factor_a;
    factor_b    = par.factor_b;
    a           = par.a;
    b           = par.b;
    equiConst   = par.equiConst;
    molFracLiq  = par.molFracLiq;
    pTot        = par.pTot;
    
    %% derive additional parameters
    molFracLiqCO = molFracLiq(:,1);
    molFracLiqH2 = molFracLiq(:,2);
    N = length(molFracLiqCO);
    
    %% perform calculations
    reactRate = zeros(N,1);
    for zPoint = 1:N
            reactRate(zPoint) = max(volFracLiq(zPoint)*factor_a*(a*equiConst(1)*molFracLiqCO(zPoint)*equiConst(2)*molFracLiqH2(zPoint)*(pTot(zPoint))^2 )/...
                ((1+b*factor_b*equiConst(1)*molFracLiqCO(zPoint)*pTot(zPoint))^2),0);
    end
 
end