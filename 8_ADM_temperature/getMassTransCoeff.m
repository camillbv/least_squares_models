%% header
%% purpose: calculate kL for CO using the relation by Calderbank and 
%% Moo-Young (1960), equation (2).

function kL = getMassTransCoeff(par)

    %% unpack parameters
    liqDensity  = par.liqDensity;
    gasDensity  = par.gasDensity;
    diffCoefCO  = par.diffCoefCO;
    viscLiq     = par.viscLiq;    
    graConst    = par.graConst;
        
    %% derive additional parameters
    N = length(gasDensity);
    
    %% in case of constant liquid density
    if size(liqDensity,1) == 1
        liqDensity = liqDensity*ones(N,1);
    end
    % verify gas and liquid density vectors are of the same length
    if size(gasDensity,1) ~= size(liqDensity,1)
        disp('error inside getMassTransCoeff: gas and slurry densities not of same length!')
    end
        
    %% perform calculations
    densityDiff = liqDensity - gasDensity;
    Sc          = viscLiq./(liqDensity.*diffCoefCO);
    kL          = 0.42./(Sc.^(1/2)).*((densityDiff.*viscLiq.*graConst./(liqDensity.^2)).^(1/3));


end