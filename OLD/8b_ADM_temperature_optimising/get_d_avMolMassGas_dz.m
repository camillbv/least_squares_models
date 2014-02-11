%% camilla.berge.vik@ntnu.no 23.01.2014
%% calculate the derivative of average gas molar mass
%% for cold flow two phase model

function d_avMolMassGas_dz = get_d_avMolMassGas_dz(par)

    %% pack out parameters
    wtFracGas       = par.wtFracGas;
    d_wtFracGas_dz  = par.d_wtFracGas_dz;
    Mw              = par.Mw;

    %% defive additional parameters
    N               = size(wtFracGas,1);
    
    %% perform calculations
    d_avMolMassGas_dz = zeros(N,1);
    
    for zPoint = 1:N
        SUMJ                        = sum(wtFracGas(zPoint,:)'./Mw);
        SUMJd                       = sum(d_wtFracGas_dz(zPoint,:)'./Mw);
        d_avMolMassGas_dz(zPoint)   = 1./(SUMJ.^2)*sum((d_wtFracGas_dz(zPoint,:)'.*SUMJ - wtFracGas(zPoint,:)'.*SUMJd));        
    end

end