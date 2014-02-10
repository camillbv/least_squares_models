%% header
% purpose: solve two-phase plug flow model
% for Fischer-Tropsch reactor
% using least-squares method
% camilla.berge.vik@ntnu.no
% 27.01.2013

%% clear memory
clc
clear all
close all

global  gasConst T ...
    nCompGas nCompLiq pTot kLa Mw ...
    liqDensityInit nu equiConst

format long

%% set parameters and initial conditions
parameters
 
%% set numerical parameters
P = 50;
N = P+1;
nVar = (3+nLumps)*2 + 2 + 2;

%% calculate points and weights in the reference domain
[z_GLL, wz_GLL]       = GaussLobattoLegendre(N);    % GLL points

%% map points and weights into physical domain
[COLUMN,wCOLUMN] = map(z_GLL,  wz_GLL,  0, totalHeight); % COLUMN length

%% derivative matrix
D_ref  = LagrangeDerivativeMatrix_GLL(N);   % reference domain -1,1
D      = 2/totalHeight*D_ref;               % physical domain   0,totalHeight

%% run ode15s to get initial estimates for weight fractions
y0 = [wtFracGasInit' wtFracLiqInit' supVelGasInit supVelLiqInit volFracLiqInit pTotInit];
[zODE,yODE] = ode15s('model_equations',COLUMN,y0);

%% use ODE estimate as initial guess
solVec = reshape(yODE,N,nVar);

%% boundary conditions
BC = y0';

%% ---- START OF ITERATION LOOP

%% set overall iteration parameters
max_iter            = 3;
tol                 = 1e-13;
iteration_error     = 1;
L2_norm_tot         = 1;
iter                = 0;

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
    
    %% iteration properties for variable iteration
    iter_variable           = 0;
    L2_norm_tot_variable    = 1;
    max_iter_variable       = 20;
    tol_variable            = 1e-12;
    
    disp('***')
    disp('***')
    disp('GAS MASS FRACTIONS')
    disp('***')
    disp('***')
    
    %% iteration loop for gas weight fractions
    while L2_norm_tot_variable > tol_variable && iter_variable < max_iter_variable
        iter_variable = iter_variable + 1;
        
        disp(['   Iteration # ',num2str(iter_variable),'...'])        

        %% pack out from solution vector
        wtFracGas   = solVec(:,1:nCompGas);
        wtFracLiq   = solVec(:,nCompGas+1:nCompGas+nCompLiq);
        supVelGas   = solVec(:,nCompGas+nCompLiq+1);
        supVelLiq   = solVec(:,nCompGas+nCompLiq+2);
        volFracLiq  = solVec(:,nCompGas+nCompLiq+3);
        pTot        = solVec(:,nCompGas+nCompLiq+4);
        
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
        
        parStructL2wgL3wg = struct('wtFracGas',wtFracGas,'wtFracLiq',wtFracLiq, ...
            'supVelGas',supVelGas, ...
            'gasDensity',gasDensity,'liqDensity',liqDensity, ...
            'kLa',kLa,'equiConst',equiConst, ...
            'avMolMassGas',avMolMassGas,'avMolMassLiq',avMolMassLiq);
        
        %% calculate L, W, B, A
        
        L1wg = zeros(N*nCompGas);
        W    = zeros(N*nCompGas);
        B    = zeros(N*nCompGas);
        for gasComp = 1:nCompGas
            ind = [(gasComp-1)*N+1:gasComp*N];
            L1wg(ind,ind) = D;
            W(ind,ind) = diag(wCOLUMN);
            B(ind(1),ind(1)) = 1;
        end
        
        [L2wg,L3wg] = get_L2wg_and_L3wg(parStructL2wgL3wg);
                    
        L = L1wg + L2wg + L3wg;

        A = L'*W*L;
        
        %% calculate g
        gMatrix = zeros(N,nCompGas);
        for zPoint = 1:N
            gMatrix(zPoint,:) = 1./supVelGas(zPoint).*liqDensity./gasDensity(zPoint).*kLa.*wtFracLiq(zPoint,:)';
        end
        g = reshape(gMatrix,N*nCompGas,1);
        
        %% calculate F
        F = L'*W*g;
        
        %% pick up F_gamma
        F_gamma = zeros(N,nCompGas);
        F_gamma(1,:) = BC(1:nCompGas);
        F_gamma = reshape(F_gamma,N*nCompGas,1);
        
        %% solve for f
        f_new = (A+B)\(F+F_gamma);
        %pause

        f_old = reshape(wtFracGas,N*nCompGas,1);
     %   pause
        
        %% underrelaxation
        underrelaxation = 0.1;
        f = f_new*underrelaxation + f_old*(1-underrelaxation);
        
        %% plot weight fraction for CO for diagnostics
        figure(1)
        subplot(2,2,1)
        plot(f(1:N))
        title('CO gas weight fractions for each iteration')
        hold on
        
        %% calculate and evaluate errors and residuals
        res_variable                = L*f - g;
        L2_norm_res_variable        = sqrt((L*f-g)'*W*(L*f-g));
        L2_norm_resBC_variable      = sqrt((B*f-F_gamma)'*W*(B*f-F_gamma));
        L2_norm_tot_variable        = L2_norm_res_variable+L2_norm_resBC_variable;
        error_variable              = f-f_old;
        iteration_error_variable    = sqrt((f-f_old)'*W*(f-f_old));
        
        %% plot some diagnostics
        subplot(2,2,2)
        semilogy(iter_variable, L2_norm_tot_variable,'ko',iter_variable,iteration_error_variable,'k*')
        title('L2 residual and iteration error')
        xlabel('iteration_variable number (gas weight fraction variable iteration only)')
        legend('L2 norm','iteration error')
        hold on
        
        subplot(2,2,3)
        plot(res_variable)
        title('local residual = L*f - g for gas weight fractions')
        xlabel('z index')
        hold on
        
        subplot(2,2,4)
        plot(error_variable)
        title('local iteration error = f-wtFracGas for gas weight fractions')
        xlabel('z index')
        hold on
        
        %% update and reshape wtFracGas into matrix again for next iteration
        wtFracGas = reshape(f,N,nCompGas);
        solVec(:,1:nCompGas)=wtFracGas;
        
    end 
    
    
    disp('***')
    disp('***')
    disp('LIQUID MASS FRACTIONS')
    disp('***')
    disp('***')
  %  pause
    
    
    %% iteration properties for variable iteration
    iter_variable           = 0;
    L2_norm_tot_variable    = 1;
    max_iter_variable       = 20;
    tol_variable            = 1e-12;
    
    %% iteration loop for liquid weight fractions
    while L2_norm_tot_variable > tol_variable && iter_variable < max_iter_variable
        
        iter_variable = iter_variable + 1;

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
        
        %% calculate reaction rate
        parStructReactRate = struct('volFracLiq',volFracLiq,'factor_a',factor_a,'factor_b',factor_b, ...
            'a',a,'b',b,'equiConst',equiConst,'molFracLiq',molFracLiq, ...
            'pTot',pTot);
        
        reactRate = get_reactRate(parStructReactRate);       
        
       parStructL2wlL3wl = struct('wtFracGas',wtFracGas,'wtFracLiq',wtFracLiq, ...
            'supVelLiq',supVelLiq,'liqDensity',liqDensity, ...
            'kLa',kLa,'equiConst',equiConst,'Mw',Mw, ...
            'avMolMassGas',avMolMassGas,'avMolMassLiq',avMolMassLiq);        
        
        %% calculate L, A
        L1wl = L1wg;
        [L2wl,L3wl] = get_L2wl_and_L3wl(parStructL2wlL3wl);
        
        L = L1wl + L2wl + L3wl;
        
        A = L'*W*L;
        
        %% calculate g
        gMatrix = zeros(N,nCompLiq);
        for zPoint = 1:N
            gMatrix(zPoint,:) = 1./supVelLiq(zPoint).*kLa*1./equiConst.* ...
                avMolMassGas(zPoint)./avMolMassLiq(zPoint).*wtFracGas(zPoint,:)' + ...
                 volFracLiq(zPoint)*nu.*reactRate(zPoint)./liqDensity./supVelLiq(zPoint).*catDensity.*Mw;
        end
        
        g = reshape(gMatrix,N*nCompLiq,1);
        
        %% calculate F
        F = L'*W*g;
        
        %% pick up F_gamma

        F_gamma = zeros(N,nCompLiq);
        F_gamma(1,:) = BC(nCompGas + 1:nCompGas + nCompLiq);
        F_gamma = reshape(F_gamma,N*nCompLiq,1);
        
        %% solve for f
        f_new = (A+B)\(F+F_gamma);
  %      f_new = max(f_new,0);
        
        %pause
        
        %% set underrelaxation here if desired
        underrelaxation = 0.1;
        f_old = reshape(wtFracLiq,N*nCompLiq,1);
        f = f_new*underrelaxation + f_old*(1-underrelaxation);
        
        %% plot current f
        figure(2)
        subplot(2,2,1)
        plot(f(1:N))
        title('CO liquid weight fractions for each iteration')
        hold on
        
        %% calculate and evaluate errors and residuals
        res_variable                = L*f - g;
        L2_norm_res_variable        = sqrt((L*f-g)'*W*(L*f-g));
        L2_norm_resBC_variable      = sqrt((B*f-F_gamma)'*W*(B*f-F_gamma));
        L2_norm_tot_variable        = L2_norm_res_variable+L2_norm_resBC_variable;
        error_variable              = f-f_old;
        iteration_error_variable    = sqrt((f-f_old)'*W*(f-f_old));
        
        %% plot some intermediate diagnostics
        subplot(2,2,2)
        semilogy(iter_variable, L2_norm_tot_variable,'ko',iter_variable,iteration_error_variable,'k*')
        title('L2 residual and iteration error')
        xlabel('iteration_variable number (liquid weight fraction variable iteration only)')
        legend('L2 norm','iteration error')
        hold on
        
        subplot(2,2,3)
        plot(res_variable)
        title('local residual = L*f - g for liquid weight fractions')
        xlabel('z index')
        hold on
        
        subplot(2,2,4)
        plot(error_variable)
        title('local iteration error = f-wtFracLiq for liquid weight fractions')
        xlabel('z index')
        hold on
        
        wtFracLiq = reshape(f,N,nCompLiq); % set new value as old before entering next iteration step
        solVec(:,nCompGas+1:nCompGas+nCompLiq)=wtFracLiq;
        
        if iter_variable == max_iter_variable
            disp(['Maximal number of iterations (',num2str(max_iter_variable),') reached for variable weight fractions'])
        end
        
        
        
    end  %% end of iteration on weight fractions
    
    
    %% update solution vector with new estimates for weight fractions
    solVec(:,1:nCompGas)=wtFracGas;
    solVec(:,nCompGas+1:nCompGas+nCompLiq)=wtFracLiq;
   
    disp('***')
    disp('***')
    disp('SUPERFICIAL VELOCITIES')
    disp('***')
    disp('***')
    
 %   pause
    
    %% pack out for superficial velocity
    wtFracGas   = solVec(:,1:nCompGas);
    wtFracLiq   = solVec(:,nCompGas+1:nCompGas+nCompLiq);
    supVelGas   = solVec(:,nCompGas+nCompLiq+1);
    supVelLiq   = solVec(:,nCompGas+nCompLiq+2);
 %   volFracLiq   = solVec(:,nCompGas+nCompLiq+3);
    pTot        = solVec(:,nCompGas+nCompLiq+4);
    
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
    
    parStructPtotDeriv = struct('totDensityInit',totDensityInit,'graConst',graConst);
    
    d_pTot_dz = get_d_pTot_dz(parStructPtotDeriv);
    
    parStructAvMolMassGasDeriv = struct('wtFracGas',wtFracGas, ...
        'd_wtFracGas_dz',d_wtFracGas_dz,'Mw',Mw);
    d_avMolMassGas_dz = get_d_avMolMassGas_dz(parStructAvMolMassGasDeriv);
    
    parStructGasDensityDeriv = struct('d_avMolMassGas_dz',d_avMolMassGas_dz, ...
        'avMolMassGas',avMolMassGas,'pTot',pTot,'gasConst',gasConst,'T',T, ...
        'd_pTot_dz',d_pTot_dz);
    
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
    F_gamma = zeros(N,2);
    F_gamma(1,:) = BC(nVar-3:nVar-2);
    F_gamma = reshape(F_gamma,N*2,1);
    
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
    underrelaxation = 0.5;
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
    solVec(:,nCompGas+nCompLiq+1) = f(1:N);
    solVec(:,nCompGas+nCompLiq+2) = f(N+1:2*N);
    
    disp('***')
    disp('***')
    disp('VOLUME FRACTIONS')
    disp('***')
    disp('***')
    
        %% pack out for superficial velocity
    wtFracGas   = solVec(:,1:nCompGas);
    wtFracLiq   = solVec(:,nCompGas+1:nCompGas+nCompLiq);
    supVelGas   = solVec(:,nCompGas+nCompLiq+1);
    supVelLiq   = solVec(:,nCompGas+nCompLiq+2);
    volFracLiq   = solVec(:,nCompGas+nCompLiq+3);
    
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
     
    %% calculate L, A, B
    L1 = D;
    sumTRANS = zeros(N,1);
    for zPoint = 1:N
        sumTRANS(zPoint) = sum(  kLa.*( 1./equiConst.*...
            avMolMassGas(zPoint)./avMolMassLiq(zPoint).*...
            wtFracGas(zPoint,:)'-wtFracLiq(zPoint,:)' ) );
    end
    
    L2 = diag(1./supVelLiq.*sumTRANS);
    
    L = L1 + L2;
    
    W = diag(wCOLUMN);
    
    A = L'*W*L;
    
    B = zeros(N,N);
    B(1,1) = 1;
    
    %% calculate g
    
    g = zeros(N,1);
    
    F = L'*W*g;
    
    F_gamma = zeros(N,1);
    F_gamma(1) = BC(nVar-1);
    
        %% solve for f
    f_new = (A+B)\(F+F_gamma);
    
    %% pick up f_old for plotting
    f_old = volFracLiq;
    
    %% set underrelaxation here
    underrelaxation = 0.5;
    f = f_new*underrelaxation + f_old*(1-underrelaxation);
        
    %% calculate and evaluate errors and residuals
    res_variable                = L*f - g;
    L2_norm_res_variable        = sqrt((L*f-g)'*W*(L*f-g));
    L2_norm_resBC_variable      = sqrt((B*f-F_gamma)'*W*(B*f-F_gamma));
    L2_norm_tot_variable        = L2_norm_res_variable+L2_norm_resBC_variable;
    error_variable              = f-f_old;
    iteration_error_variable    = sqrt((f-f_old)'*W*(f-f_old));
    
    %% plot
    figure(4)
    subplot(3,1,1)
    plot(f_new(1:N))
    title('liquid volume fraction')
    hold on
    
    subplot(3,1,2)
    plot(res_variable)
    title('local residual = L*f - g for liquid volume fraction')
    xlabel('z index')
    hold on
    
    subplot(3,1,3)
    plot(error_variable)
    title('local iteration error = f_new-f_old for liquid volume fraction')
    xlabel('z index')
    hold on
    
    %% update solution vector
    solVec(:,nCompGas+nCompLiq+3) = f(1:N);

 %   pause
    
    disp('***')
    disp('***')
    disp('PRESSURE')
    disp('***')
    disp('***')
    
    % no need to solve 
    pTot = pTotInit + d_pTot_dz*COLUMN;
    
end

%% unpack and plot results
%% unpack results
wRES_G = solVec(:,1:nCompGas);
wRES_L = solVec(:,nCompGas+1:nCompGas+nCompLiq);
vRES_G = solVec(:,nCompGas+nCompLiq + 1);
vRES_L = solVec(:,nCompGas+nCompLiq + 2);
epsRES_L    = solVec(:,nCompGas+nCompLiq + 3);
pRES        = solVec(:,nCompGas+nCompLiq + 4);

%% make plots
figure(5)
plot(COLUMN,wRES_G(:,1),...
     COLUMN,wRES_G(:,2),...      
     COLUMN,wRES_G(:,3),...
     COLUMN,wRES_G(:,4),'g--',...
     COLUMN,wRES_G(:,5),'r-.',...
     COLUMN,wRES_G(:,6),'b-', ...
     COLUMN,wRES_G(:,7),'k:' )
title('gas mass fractions')
legend('CO','H_2','H2O','C1-C10','C11-C20','C21-C30','C31+')
ylabel('gas phase mass fractions')
xlabel('reactor length')
grid on

figure(6)
plot(COLUMN,wRES_L(:,1),...
     COLUMN,wRES_L(:,2),...      
     COLUMN,wRES_L(:,3),...
     COLUMN,wRES_L(:,4),'g--',...
     COLUMN,wRES_L(:,5),'r-.',...
     COLUMN,wRES_L(:,6),'b-', ...
     COLUMN,wRES_L(:,7),'k:' )
title('liquid mass fractions')
legend('CO','H_2','H2O','C1-C10','C11-C20','C21-C30','C31+')
ylabel('liquid phase mass fractions')
xlabel('reactor length')
grid on

%% calculate conversion of CO (%) on mass basis
wconversion     = 100*wRES_G(:,1);
wcum_conversion = 100*(wtFracGasInit(1) - wRES_G(:,1))./wtFracGasInit(1);

figure(7)
plot(COLUMN,wconversion,COLUMN,wcum_conversion,'r')
title('Conversion on mass basis (%)')
ylabel('conversion (%)')
xlabel('reactor length')
legend('Conversion','Cumulative Conversion')

figure(8)
plot(COLUMN,vRES_G)
title('Superficial gas velocity')

figure(9)
plot(COLUMN,vRES_L)
title('Superficial liquid velocity')
xlabel('reactor length')
ylabel('liquid velocity')

figure(10)
plot(COLUMN,epsRES_L)
title('Liquid volume fraction')
xlabel('reactor length')
ylabel('volume fraction [-]')
figure(11)
plot(COLUMN,pRES)
title('Pressure profile')
xlabel('reactor length')
ylabel('pressure [Pa]')
