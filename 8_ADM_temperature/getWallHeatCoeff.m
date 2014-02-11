function hW = getWallHeatCoeff(par)

%% unpack parameters
    condLiq     = par.condLiq;
    condSol     = par.condSol;
    volFracSol  = par.volFracSol;
    viscLiq     = par.viscLiq;
    graConst    = par.graConst;
    supVelGas   = par.supVelGas;
    sluDensity  = par.sluDensity;
    einsteinK   = par.einsteinK; 
    heatCapSlu  = par.heatCapSlu;
    
    %% calculate slurry viscosity
    viscSlu = viscLiq*(1+einsteinK*volFracSol)
    
    %% calculate thermal conductivity
    condSlu = condLiq* ...
             (2*condLiq+condSol-2*volFracSol*(condLiq-condSol))/ ...
             (  condLiq+condSol-  volFracSol*(condLiq-condSol))
         
    %% calculate hW
    hW = 0.1*(condSlu^0.5*sluDensity^0.75*heatCapSlu^0.5* ...
              viscSlu^(-0.25)*graConst^0.25*supVelGas^0.25)
         
end