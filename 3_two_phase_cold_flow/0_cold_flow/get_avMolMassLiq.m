%% camilla.berge.vik@ntnu.no 22.01.2014
%% calculate average liquid molar mass in each spatial point

function avMolMassLiq = get_avMolMassLiq(par)

    %% unpack parameters
    molFracLiq = par.molFracLiq;
    Mw = par.Mw;
    
    %% derive additional parameters
    N = size(molFracLiq,1);
    
    %% perform calculations for each spatial point
    avMolMassLiq = zeros(N,1);
    for zPoint = 1:N
        avMolMassLiq(zPoint) = molFracLiq(zPoint,:)*Mw;
    end
    
end