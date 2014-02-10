%% header
% purpose: solve two-phase plug flow model
% for Fischer-Tropsch reactor
% using least-squares method
% camilla.berge.vik@ntnu.no
% 22.01.2013

%% clear memory
clc
clear all
close all

%% get parameters
parameters

global  gasConst T ...
    nCompGas nCompLiq pTot kLa Mw ...
    liqDensityInit nu equiConst

%% set numerical parameters
P = 30;
N = P+1;
nVariables = 70;

%% inlet gas phase concentrations
molFracGasCOInit    = 0.25;
molFracGasH2Init    = 0.75;
molFracGasH2OInit   = 0.0;
molFracGasAlkInit   = zeros(Cmax+1,1)';
molFracGasInit      = [molFracGasCOInit molFracGasH2Init molFracGasH2OInit molFracGasAlkInit]';
molFracGasInit      = molFracGasInit./sum(molFracGasInit);  % scale it to sum to 1

%% inlet liquid phase concentrations
molFracLiqCOInit    = 0;
molFracLiqH2Init    = 0;
molFracLiqH2OInit   = 0;
molFracLiqAlkInit   = nu(4:length(nu)).*ones(Cmax+1,1);    % prod. gen. by one kmole reacted CO
molFracLiqInit      = [molFracLiqCOInit molFracLiqH2Init molFracLiqH2OInit molFracLiqAlkInit']';     % add water
molFracLiqInit      = molFracLiqInit./sum(molFracLiqInit);  % scale it to sum to 1

%% convert to weight fractions
wtFracGasInit   = molFracGasInit.*Mw./(molFracGasInit'*Mw); %  kg/kg (unitless)  gas feed weight fractions
wtFracLiqInit   = molFracLiqInit.*Mw./(molFracLiqInit'*Mw); %  kg/kg (unitless)  gas feed weight fractions

%% number of components
nCompGas = length(wtFracGasInit);
nCompLiq = length(wtFracLiqInit);

%% feed average molar mass
avMolMassGasInit   = molFracGasInit'*Mw; % kg/kmol gas feed average molar mass
avMolMassLiqInit   = molFracLiqInit'*Mw; % kg/kmol liq feed average molar mass

%% feed gas density
gasDensityInit    = pTot*avMolMassGasInit/gasConst/T; % kg/m^3  gas feed total density
liqDensityInit    = 600; % estimate based on HYSYS flash simulations

%% calculate points and weights in the reference domain
[z_GLL, wz_GLL]       = GaussLobattoLegendre(N);    % GLL points

%% map points and weights into physical domain
[COLUMN,wCOLUMN] = map(z_GLL,  wz_GLL,  0, totalHeight); % COLUMN length

%% derivative matrix
D_ref  = LagrangeDerivativeMatrix_GLL(N);   % reference domain -1,1
D      = 2/totalHeight*D_ref;                         % physical domain   0,totalHeight

%% calculate constant parameters throughout iteration:
L = D;
W = diag(wCOLUMN);
A = L'*W*L;
B = zeros(N,N);
B(1,1)=1;

LLw = zeros(N*nCompGas,N*nCompGas);
WWw = zeros(N*nCompGas,N*nCompGas);
BBw = zeros(N*nCompGas,N*nCompGas);
for varNum = 1:nCompGas
    indices = (varNum-1)*N+1:varNum*N;
    LLw(indices,indices)=L;
    WWw(indices,indices)=W;
    BBw(indices,indices)=B;
    
end
AAw = LLw'*WWw*LLw;

%% cold flow
reactRate = zeros(N,1);

%% run ode15s to get initial estimates for weight fractions

%% assemble initial conditions vector
y0 = [wtFracGasInit' wtFracLiqInit' supVelGasInit supVelLiqInit];

%% set height
heightValues    = [0:totalHeight/(N-1):totalHeight];

%% call ode15s
[zODE,yODE] = ode15s('model_equations',heightValues,y0);

%% assemble initial conditions vector
globalfInit = yODE;

bigFg = zeros(N,nVariables);
for varNum = 1:nVariables
    bigFg(1,varNum) = globalfInit(1,varNum);
end
FFgwG = bigFg(:,1:nCompGas);
FFgwL = bigFg(:,nCompGas+1:nCompGas+nCompLiq);
FFgvG = bigFg(:,nCompGas+nCompLiq+1);
FFgvL = bigFg(:,nCompGas+nCompLiq+2);

%% reshape F's into long vectors instead of matrices
FFgwG = reshape(FFgwG,N*nCompGas,1);
FFgwL = reshape(FFgwL,N*nCompLiq,1);

%% ---- START OF ITERATION LOOP

%% set overall iteration parameters
max_iter            = 10;
tol                 = 1e-13;
iteration_error     = 1;
L2_norm_tot         = 1;
iter                = 0;

solVec              = globalfInit; % initiate solution vector

%% print iteration message to screen
disp(['Starting ',num2str(max_iter), ' iterations in overall loop ... '])
disp(['Specified tolerance is ',num2str(tol),'.'])
disp('...')
tic

while  sum(L2_norm_tot) > tol && iter < max_iter % criteria to continue iteration
    
    %% increase iteration number
    iter                = iter+1;
    
    %% display iteration message
    disp('...')
    disp(['Starting iteration # ',num2str(iter),'...'])
    
    fB = solVec; % pick up old value to calculate errors at end of each overall iteration
    % fB contains all variables in all points; it is a matrix with dimension N*nVariables !!
    
    
    disp('***')
    disp('***')
    disp('GAS MASS FRACTIONS')
    disp('***')
    disp('***')
    
    
    %% pack out from solution vector
    wtFracGas   = solVec(:,1:nCompGas);
    wtFracLiq   = solVec(:,nCompGas+1:nCompGas+nCompLiq);
    supVelGas   = solVec(:,nCompGas+nCompLiq+1);
    supVelLiq   = solVec(:,nCompGas+nCompLiq+2);
    
    %% fetch appropriate L, W, F_gamma, B
    L = LLw;
    W = WWw;
    B = BBw;
    A = L'*W*L;
    
    %% iteration properties for variable iteration
    iter_variable           = 0;
    L2_norm_tot_variable    = 1;
    max_iter_variable       = 30;
    tol_variable            = 1e-12;
    
    %% iteration loop for gas mass fractions    
    while L2_norm_tot_variable > tol_variable && iter_variable < max_iter_variable
        iter_variable = iter_variable + 1;
        disp('   ...')
        disp(['   Iteration # ',num2str(iter_variable),'...'])
        
        %% convert to mole fractions
        molFracGas = zeros(N,nCompGas);
        molFracLiq = zeros(N,nCompLiq);
        for zPoint = 1:N
            molFracGas(zPoint,:) = wtFracGas(zPoint,:)'./Mw/(sum(wtFracGas(zPoint,:)'./Mw));
            molFracLiq(zPoint,:) = wtFracLiq(zPoint,:)'./Mw/(sum(wtFracLiq(zPoint,:)'./Mw));
        end
        
        %% update parameters dependent on z or the solution vector
        parStructAvMolMassGas = struct('molFracGas',molFracGas,'Mw',Mw);
        parStructAvMolMassLiq = struct('molFracLiq',molFracLiq,'Mw',Mw);
        avMolMassGas = get_avMolMassGas(parStructAvMolMassGas);
        avMolMassLiq = get_avMolMassLiq(parStructAvMolMassLiq);
        
        parStructGasDensity = struct('avMolMassGas',avMolMassGas,'pTot',pTot,'gasConst',gasConst,'T',T);
        gasDensity = get_gasDensity(parStructGasDensity);
        liqDensity = liqDensityInit; % constant density
        
        parStructWtFracGasDeriv = struct('wtFracGas',wtFracGas,'wtFracLiq',wtFracLiq, ...
            'supVelGas',supVelGas, ...
            'gasDensity',gasDensity,'liqDensity',liqDensity, ...
            'kLa',kLa,'equiConst',equiConst, ...
            'avMolMassGas',avMolMassGas,'avMolMassLiq',avMolMassLiq);
        d_wtFracGas_dz=get_d_wtFracGas_dz(parStructWtFracGasDeriv);
        
        %% calculate g 
        g = reshape(d_wtFracGas_dz,N*nCompGas,1);
        
        %% calculate F
        F = L'*W*g;
        
        %% pick up F_gamma
        F_gamma = FFgwG;
        
        %% solve for f
        f_new = (A+B)\(F+F_gamma);
        
        %% underrelaxation
        underrelaxation = 0.1;
        f = f_new*underrelaxation + reshape(wtFracGas,N*nCompGas,1)*(1-underrelaxation);
        
        %% plot mass fraction for CO for diagnostics
        figure(1)
        subplot(2,2,1)
        plot(f(1:N))
        title('CO gas mass fractions for each iteration')
        hold on
        
        %% re-shape wtfracgas to long vector for residual calculation
        wtFracGas = reshape(wtFracGas,N*nCompGas,1);
        
        %% calculate and evaluate errors and residuals
        res_variable                = L*f - g;
        L2_norm_res_variable        = sqrt((L*f-g)'*W*(L*f-g));
        L2_norm_resBC_variable      = sqrt((B*f-F_gamma)'*W*(B*f-F_gamma));
        L2_norm_tot_variable        = L2_norm_res_variable+L2_norm_resBC_variable;
        error_variable              = f-wtFracGas;
        iteration_error_variable    = sqrt((f-wtFracGas)'*W*(f-wtFracGas));
                
        %% plot some diagnostics
        subplot(2,2,2)
        semilogy(iter_variable, L2_norm_tot_variable,'ko',iter_variable,iteration_error_variable,'k*')
        title('L2 residual and iteration error')
        xlabel('iteration_variable number (gass mass variable iteration only)')
        legend('L2 norm','iteration error')
        hold on
        
        subplot(2,2,3)
        plot(res_variable)
        title('local residual = L*f - g for gas mass fractions')
        xlabel('z index')
        hold on
        
        subplot(2,2,4)
        plot(error_variable)
        title('local iteration error = f-wtFracGas for gas mass fractions')
        xlabel('z index')
        hold on
        
        %% update and reshape wtFracGas into matrix again for next iteration
        wtFracGas = reshape(f,N,nCompGas); 
        
        if iter_variable == max_iter_variable
            disp(['Maximal number of iterations (',num2str(max_iter_variable),') reached for variable # ',num2str(varNum),'.'])
        end
        
    end  %% end of iteration on gas mass fractions
    
    disp('***')
    disp('***')
    disp('LIQUID MASS FRACTIONS')
    disp('***')
    disp('***')
    
    %% update solution vector with new estimates for gas weight fractions
    solVec(:,1:nCompGas)=wtFracGas;
    
    %% pack out for liquid weight fractions
    wtFracGas   = solVec(:,1:nCompGas);
    wtFracLiq   = solVec(:,nCompGas+1:nCompGas+nCompLiq);
    supVelGas   = solVec(:,nCompGas+nCompLiq+1);
    supVelLiq   = solVec(:,nCompGas+nCompLiq+2);
    
    %% iteration properties for variable iteration
    iter_variable           = 0;
    L2_norm_tot_variable    = 1;
    max_iter_variable       = 30;
    tol_variable            = 1e-12;
    
    while L2_norm_tot_variable > tol_variable && iter_variable < max_iter_variable
        
        iter_variable = iter_variable + 1;
        disp('   ...')
        disp(['   Iteration # ',num2str(iter_variable),' for liquid mass fractions'])
        

        
        %% convert to mole fractions
        molFracGas = zeros(N,nCompGas);
        molFracLiq = zeros(N,nCompLiq);
        
        for zPoint = 1:N
            molFracGas(zPoint,:) = wtFracGas(zPoint,:)'./Mw/(sum(wtFracGas(zPoint,:)'./Mw));
            molFracLiq(zPoint,:) = wtFracLiq(zPoint,:)'./Mw/(sum(wtFracLiq(zPoint,:)'./Mw));
        end
        
        %% update parameters dependent on z or the solution vector
        
        parStructAvMolMassGas = struct('molFracGas',molFracGas,'Mw',Mw);
        parStructAvMolMassLiq = struct('molFracLiq',molFracLiq,'Mw',Mw);
        
        
        avMolMassGas = get_avMolMassGas(parStructAvMolMassGas);
        avMolMassLiq = get_avMolMassLiq(parStructAvMolMassLiq);
        
        
        parStructWtFracLiqDeriv = struct('wtFracGas',wtFracGas,'wtFracLiq',wtFracLiq, ...
            'supVelLiq',supVelLiq,'liqDensity',liqDensity, ...
            'kLa',kLa,'equiConst',equiConst,'Mw',Mw, ...
            'avMolMassGas',avMolMassGas,'avMolMassLiq',avMolMassLiq, ...
            'nu',nu,'reactRate',reactRate,'catDensity',catDensity);
        
        d_wtFracLiq_dz = get_d_wtFracLiq_dz(parStructWtFracLiqDeriv);
        
        %% calculate g
        g = reshape(d_wtFracLiq_dz,N*nCompLiq,1);
        
        %% calculate F
        F = L'*W*g;
        
        %% pick up F_gamma
        F_gamma = FFgwL;
        
        %% solve for f
        f_new = (A+B)\(F+F_gamma);
        
        %% set underrelaxation here if desired
        underrelaxation = 0.02;
        f = f_new*underrelaxation + reshape(wtFracLiq,N*nCompLiq,1)*(1-underrelaxation);
        
        %% plot current f
        figure(2)
        subplot(2,2,1)
        plot(f(1:N))
        title('CO gas mass fractions for each iteration')
        hold on
        
        %% re-shape wtfracgas for residual calculation
        wtFracLiq = reshape(wtFracLiq,N*nCompLiq,1);
        
        %% calculate and evaluate errors and residuals
        res_variable                = L*f - g;
        L2_norm_res_variable        = sqrt((L*f-g)'*W*(L*f-g));
        L2_norm_resBC_variable      = sqrt((B*f-F_gamma)'*W*(B*f-F_gamma));
        L2_norm_tot_variable        = L2_norm_res_variable+L2_norm_resBC_variable;
        error_variable              = f-wtFracLiq;
        iteration_error_variable    = sqrt((f-wtFracLiq)'*W*(f-wtFracLiq));
        
        
        %% plot some intermediate diagnostics
        
        subplot(2,2,2)
        semilogy(iter_variable, L2_norm_tot_variable,'ko',iter_variable,iteration_error_variable,'k*')
        title('L2 residual and iteration error')
        xlabel('iteration_variable number (liquid mass variable iteration only)')
        legend('L2 norm','iteration error')
        hold on
        
        subplot(2,2,3)
        plot(res_variable)
        title('local residual = L*f - g for liquid mass fractions')
        xlabel('z index')
        hold on
        
        subplot(2,2,4)
        plot(error_variable)
        title('local iteration error = f-wtFracLiq for liquid mass fractions')
        xlabel('z index')
        hold on
        
        wtFracLiq = reshape(f,N,nCompLiq); % set new value as old before entering next iteration step
        
        if iter_variable == max_iter_variable
            disp(['Maximal number of iterations (',num2str(max_iter_variable),') reached for gas mass fractions'])
        end
        
    end %% end of liquid mass fractions iterations
    
    %% update solution vector with new estimates for gas weight fractions
    solVec(:,nCompGas+1:nCompGas+nCompLiq)=wtFracLiq;
    
%     disp('One round of gas mass iterations and one round of liquid mass iterations finished.')
%     disp('Press a button to continue to superficial velocity')
%     pause
    
    disp('***')
    disp('***')
    disp('SUPERFICIAL VELOCITIES')
    disp('***')
    disp('***')
    
    %% pack out for superficial velocity
    wtFracGas   = solVec(:,1:nCompGas);
    wtFracLiq   = solVec(:,nCompGas+1:nCompGas+nCompLiq);
    supVelGas   = solVec(:,nCompGas+nCompLiq+1);
    supVelLiq   = solVec(:,nCompGas+nCompLiq+2);
    
    %% convert to mole fractions
    molFracGas = zeros(N,nCompGas);
    molFracLiq = zeros(N,nCompLiq);
    
    for zPoint = 1:N
        molFracGas(zPoint,:) = wtFracGas(zPoint,:)'./Mw/(sum(wtFracGas(zPoint,:)'./Mw));
        molFracLiq(zPoint,:) = wtFracLiq(zPoint,:)'./Mw/(sum(wtFracLiq(zPoint,:)'./Mw));
    end
    
    %% update parameters dependent on z or the solution vector
    parStructAvMolMassGas = struct('molFracGas',molFracGas,'Mw',Mw);
    parStructAvMolMassLiq = struct('molFracLiq',molFracLiq,'Mw',Mw);
    
    avMolMassGas = get_avMolMassGas(parStructAvMolMassGas);
    avMolMassLiq = get_avMolMassLiq(parStructAvMolMassLiq);
    parStructGasDensity = struct('avMolMassGas',avMolMassGas,'pTot',pTot,'gasConst',gasConst,'T',T);
    gasDensity = get_gasDensity(parStructGasDensity);
    liqDensity = liqDensityInit; % constant density
    
    parStructWtFracGasDeriv = struct('wtFracGas',wtFracGas,'wtFracLiq',wtFracLiq, ...
        'supVelGas',supVelGas, ...
        'gasDensity',gasDensity,'liqDensity',liqDensity, ...
        'kLa',kLa,'equiConst',equiConst, ...
        'avMolMassGas',avMolMassGas,'avMolMassLiq',avMolMassLiq);
    
    d_wtFracGas_dz=get_d_wtFracGas_dz(parStructWtFracGasDeriv);
    
    parStructAvMolMassGasDeriv = struct('wtFracGas',wtFracGas, ...
        'd_wtFracGas_dz',d_wtFracGas_dz,'Mw',Mw);
    d_avMolMassGas_dz = get_d_avMolMassGas_dz(parStructAvMolMassGasDeriv);
    
    parStructGasDensityDeriv = struct('d_avMolMassGas_dz',d_avMolMassGas_dz, ...
        'pTot',pTot,'gasConst',gasConst,'T',T);
    
    d_gasDensity_dz = get_d_gasDensity_dz(parStructGasDensityDeriv);
    
    
    %% calculate L
    L = zeros(2*N,2*N);
    
    %% calculate L for gas superficial velocity
    L(1:N,1:N) = D + diag(1./gasDensity.*d_gasDensity_dz);
        
    %% calculate L for liquid superficial velocity
    L(N+1:2*N,N+1:2*N)=D;
    
    %% calculate W
    W = zeros(2*N,2*N);
    W(1:N,1:N) = diag(wCOLUMN);
    W(N+1:2*N,N+1:2*N)=diag(wCOLUMN);
    
    %% calculate A
    A = L'*W*L;
    
    %% calculate B
    B = zeros(2*N,2*N);
    B(1,1) = 1;
    B(N+1,N+1)=1;
    
    %% pick up F_gamma
    F_gamma = reshape(bigFg(:,69:70),2*N,1);
    
    %% calculate g
    
    parStructSupVelGasRHSDeriv = struct('gasDensity',gasDensity,...
        'liqDensity',liqDensity,'kLa',kLa,'equiConst',equiConst, ...
        'avMolMassGas',avMolMassGas,'avMolMassLiq',avMolMassLiq,...
        'wtFracGas',wtFracGas,'wtFracLiq',wtFracLiq);
    
    d_supVelGasRHS_dz = get_d_supVelGasRHS_dz(parStructSupVelGasRHSDeriv);
    
    parStructSupVelLiqDeriv = struct('kLa',kLa,'equiConst',equiConst, ...
        'avMolMassGas',avMolMassGas,'avMolMassLiq',avMolMassLiq,...
        'wtFracGas',wtFracGas,'wtFracLiq',wtFracLiq);
    
    d_supVelLiq_dz = get_d_supVelLiq_dz(parStructSupVelLiqDeriv);
    
    gG = d_supVelGasRHS_dz;
    gL = d_supVelLiq_dz;
    
    g = [gG; gL];
    
    %% calculate F
    F = L'*W*g;
    
    %% solve for f
    f_new = (A+B)\(F+F_gamma);
    
    %% pick up f_old for plotting
    f_old = [supVelGas; supVelLiq];
    
    %% set underrelaxation here
    underrelaxation = 0.05;
    f = f_new*underrelaxation + f_old*(1-underrelaxation);
    
    
    %% calculate and evaluate errors and residuals
    res_variable                = L*f - g;
    L2_norm_res_variable        = sqrt((L*f-g)'*W*(L*f-g));
    L2_norm_resBC_variable      = sqrt((B*f-F_gamma)'*W*(B*f-F_gamma));
    L2_norm_tot_variable        = L2_norm_res_variable+L2_norm_resBC_variable;
    error_variable              = f-f_old;
    iteration_error_variable    = sqrt((f-f_old)'*W*(f-f_old));
    
    %% plot    
    figure(3)
    subplot(2,2,1)
    plot(f_new(1:N))
    title('Superficial gas velocity')
    hold on
    
    subplot(2,2,2)
    plot(f_new(N+1:2*N))
    title('Superficial liquid velocity')
    hold on
    
    subplot(2,2,3)
    plot(res_variable)
    title('local residual = L*f - g for superficial velocities')
    xlabel('z index')
    hold on
    
    subplot(2,2,4)
    plot(error_variable)
    title('local iteration error = f_new-f_old for superficial velocities')
    xlabel('z index')
    hold on
    
    %% update solution vector
    solVec(:,nCompGas+nCompLiq+1) = f_new(1:N);
    solVec(:,nCompGas+nCompLiq+2) = f_new(N+1:2*N);
    
    
end

%% unpack and plot results
%% unpack results
wRES_G_LSQ = solVec(:,1:nCompGas);               
wRES_L_LSQ = solVec(:,nCompGas+1:nCompGas+nCompLiq); 
vRES_G_LSQ = solVec(:,nCompGas+nCompLiq + 1);    
vRES_L_LSQ = solVec(:,nCompGas+nCompLiq + 2); 

%% lump products for better readability in plots
for iz = 1:size(wRES_G_LSQ,1)
    % gas
    wRES_C1toC10_G_LSQ(iz) =  sum(wRES_G_LSQ(iz,4:13));
    wRES_C11toC20_G_LSQ(iz) = sum(wRES_G_LSQ(iz,14:23));
    wRES_C21toC30_G_LSQ(iz) = sum(wRES_G_LSQ(iz,24:33));
    wRES_C31plus_G_LSQ(iz)  = sum(wRES_G_LSQ(iz,34));
    
    % liquid
    wRES_C1toC10_L_LSQ(iz)  = sum(wRES_L_LSQ(iz,4:13));
    wRES_C11toC20_L_LSQ(iz) = sum(wRES_L_LSQ(iz,14:23));
    wRES_C21toC30_L_LSQ(iz) = sum(wRES_L_LSQ(iz,24:33));
    wRES_C31plus_L_LSQ(iz)  = sum(wRES_L_LSQ(iz,34));
end  

%% mass fractions, gas
figure(4)
plot(COLUMN,wRES_G_LSQ(:,1),...
     COLUMN,wRES_G_LSQ(:,2),...      
     COLUMN,wRES_G_LSQ(:,3),...
     COLUMN,wRES_C1toC10_G_LSQ,'g--',...
     COLUMN,wRES_C11toC20_G_LSQ,'r-.',...
     COLUMN,wRES_C21toC30_G_LSQ,'b-', ...
     COLUMN,wRES_C31plus_G_LSQ,'k:' )
title('gas mass fractions')
legend('CO','H_2','H2O','C1-C10','C11-C20','C21-C30','C31+')
ylabel('gas phase mass fractions')
xlabel('reactor length')
grid on

%% mass fractions, liquid
figure(5)
plot(COLUMN,wRES_L_LSQ(:,1),...
     COLUMN,wRES_L_LSQ(:,2),...      
     COLUMN,wRES_L_LSQ(:,3),...
     COLUMN,wRES_C1toC10_L_LSQ,'g--',...
     COLUMN,wRES_C11toC20_L_LSQ,'r-.',...
     COLUMN,wRES_C21toC30_L_LSQ,'b-', ...
     COLUMN,wRES_C31plus_L_LSQ,'k:' )
title('liquid mass fractions')
legend('CO','H_2','H2O','C1-C10','C11-C20','C21-C30','C31+')
ylabel('liquid phase mass fractions')
xlabel('reactor length')
grid on

%% gas velocity
figure(6)
plot(COLUMN,vRES_G_LSQ)
title('Superficial gas velocity')
xlabel('reactor length [m]')
ylabel('velocity [m/s]')

%% liquid velocity
figure(7)
plot(COLUMN,vRES_L_LSQ)
title('Superficial liquid velocity')
xlabel('reactor length [m]')
ylabel('velocity [m/s]')

%% liquid mass flux
figure(8)
plot(COLUMN,vRES_L_LSQ*liqDensityInit,'r')
title('Liquid mass flux [kg/s]')
ylabel('liquid mass flux (v_LS*rho_L [kg/(m^2s)])')
xlabel('reactor length [m]')

