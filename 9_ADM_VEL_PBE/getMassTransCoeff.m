%% header
%% purpose: calculate kL for CO using the relation by Calderbank and 
%% Moo-Young (1960), equation (2).

function kL = getMassTransCoeff(par)

    %% unpack parameters
    liqDensity   = par.liqDensity;
    gasDensity   = par.gasDensity;
    viscLiq      = par.viscLiq;    
    graConst     = par.graConst;
    molecularDia = par.molecularDia; 
    Mw           = par.Mw;
    mWoctacosane = par.mWoctacosane;
    mDoctacosane = par.mDoctacosane;
    Nparameter   = par.Nparameter;
    tempSlu      = par.tempSlu;
    
    %% derive additional parameters
    N = length(gasDensity);
    nComp = length(Mw);
    
    %% in case of constant liquid density
    if size(liqDensity,1) == 1
        liqDensity = liqDensity*ones(N,1);
    end
    % verify gas and liquid density vectors are of the same length
    if size(gasDensity,1) ~= size(liqDensity,1)
        disp('error inside getMassTransCoeff: gas and slurry densities not of same length!')
    end
        
    %% perform calculations
    % calculate diffusion coefficients of each species in liquid at current
    % temperature, according to Erkey, Rodden and Akgerman (1990)
    
    beta        = 94.5./(Mw.^0.239)./(mWoctacosane^0.781)./((molecularDia.*mDoctacosane).^1.134);
    b           = 1.206 + 0.0632*(molecularDia./mDoctacosane);
    Vnull       = (Nparameter.*mDoctacosane^3)/sqrt(2);
    vD          = b.*Vnull;
    diffCoef    = zeros(N,nComp);
    molarVolume = mWoctacosane./liqDensity*1000;
    
    for zP = 1:N
        diffCoef(zP,:)    = 10^-9*sqrt(tempSlu(zP)).*beta.*(molarVolume(zP)-vD);
    end
        
    densityDiff = liqDensity - gasDensity;
    
    kL = zeros(N,nComp);
    for zP = 1:N
        Sc          = viscLiq./(liqDensity(zP).*diffCoef(zP,:));
        kL(zP,:)    = 0.42./(Sc.^(1/2)).* ...
            (densityDiff(zP).*viscLiq.*graConst./liqDensity(zP).^2).^(1/3);
    end

end