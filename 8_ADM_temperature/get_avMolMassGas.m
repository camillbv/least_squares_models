%% camilla.berge.vik@ntnu.no 22.01.2014
%% calculate average gas molar mass in each spatial point

function avMolMassGas = get_avMolMassGas(par)

    %% unpack parameters
    molFracGas = par.molFracGas;
    Mw = par.Mw;
    
    %% derive additional parameters
    N = size(molFracGas,1);
    
    %% perform calculations for each spatial point
    avMolMassGas = zeros(N,1);
    for zPoint = 1:N
        avMolMassGas(zPoint) = molFracGas(zPoint,:)*Mw;
    end
    
end