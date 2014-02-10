%% camilla.berge.vik@ntnu.no 28.01.2014
%% calculate the derivative of the total pressure
%% in each spatial point
function d_pTot_dz = get_d_pTot_dz(par)

    %% unpack parameters
    totDensityInit = par.totDensityInit;
    graConst = par.graConst;

    %% perform calculations
    d_pTot_dz = -totDensityInit*graConst;
end