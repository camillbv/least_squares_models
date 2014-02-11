%% camilla.berge.vik@ntnu.no 24.01.2014
%% calculate the reaction rate in each spatial point
function reactRate = get_reactRate(par)
    %% unpack variables
    aS          = par.aS;
    bS          = par.bS;
    equiConst   = par.equiConst;
    molFracLiq  = par.molFracLiq;
    pTot        = par.pTot;
    
    %% derive additional parameters
    molFracLiqCO = molFracLiq(:,1);
    molFracLiqH2 = molFracLiq(:,2);
    N = length(molFracLiqCO);
    
    %% perform calculations
    reactRate = zeros(N,1);
    for zP = 1:N
            reactRate(zP) = max((aS(zP)*equiConst(1)*molFracLiqCO(zP)*equiConst(2)*molFracLiqH2(zP)*pTot^2 )/...
                ((1+bS(zP)*equiConst(1)*molFracLiqCO(zP)*pTot)^2),0);
    end
 
end