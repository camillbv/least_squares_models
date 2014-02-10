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
P = 4;
N = P+1;
nVar = (3+nLumps)*2 + 2;

%% calculate points and weights in the reference domain
[z_GLL, wz_GLL]       = GaussLobattoLegendre(N);    % GLL points

%% map points and weights into physical domain
[COLUMN,wCOLUMN] = map(z_GLL,  wz_GLL,  0, totalHeight); % COLUMN length

%% derivative matrix
D_ref  = LagrangeDerivativeMatrix_GLL(N);   % reference domain -1,1
D      = 2/totalHeight*D_ref;               % physical domain   0,totalHeight

%% run ode15s to get initial estimates for weight fractions
y0 = [wtFracGasInit' wtFracLiqInit' supVelGasInit supVelLiqInit];
[zODE,yODE] = ode15s('model_equations',COLUMN,y0);

%% use ODE estimate as initial guess
solVec = reshape(yODE,N,nVar);

%% boundary conditions
BC = y0';

%% set overall iteration parameters
max_iter            = 3;
tol                 = 1e-6;
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
    
    %% iteration properties for weight fractions
    iter_var           = 0;
    L2_norm_tot_var    = 1;
    max_iter_var       = 20;
    tol_var            = 1e-6;
    
    %% -- weight fractions
    while L2_norm_tot_var > tol_var && iter_var < max_iter_var
        iter_var = iter_var + 1;
        
        disp(['   Iteration # ',num2str(iter_var),'...'])        

        %% pack out 
        wtFracGas   = solVec(:,1:nCompGas);
        wtFracLiq   = solVec(:,nCompGas+1:nCompGas+nCompLiq);
        supVelGas   = solVec(:,nCompGas+nCompLiq+1);
        supVelLiq   = solVec(:,nCompGas+nCompLiq+2);
        
        %% mole fractions
        molFracGas = zeros(N,nCompGas);
        molFracLiq = zeros(N,nCompLiq);
        for zPoint = 1:N
            molFracGas(zPoint,:) = wtFracGas(zPoint,:)'./Mw/(sum(wtFracGas(zPoint,:)'./Mw));
            molFracLiq(zPoint,:) = wtFracLiq(zPoint,:)'./Mw/(sum(wtFracLiq(zPoint,:)'./Mw));
        end
        
        parStructAvMolMassGas = struct('molFracGas',molFracGas,'Mw',Mw);
        parStructAvMolMassLiq = struct('molFracLiq',molFracLiq,'Mw',Mw);
        avMolMassGas = get_avMolMassGas(parStructAvMolMassGas);
        avMolMassLiq = get_avMolMassLiq(parStructAvMolMassLiq);
        
        parStructGasDensity = struct('avMolMassGas',avMolMassGas,'pTot',pTot,'gasConst',gasConst,'T',T);
        gasDensity = get_gasDensity(parStructGasDensity);
        liqDensity = liqDensityInit; % constant density        
        
        d_avMolMassGas_dz = D*avMolMassGas;
        parStructGasDensityDeriv = struct('d_avMolMassGas_dz',d_avMolMassGas_dz, ...
            'pTot',pTot,'gasConst',gasConst,'T',T);
        
        d_gasDensity_dz = get_d_gasDensity_dz(parStructGasDensityDeriv);
        parStructSupVelGasDeriv = struct('gasDensity',gasDensity,...
            'liqDensity',liqDensity,'kLa',kLa,'equiConst',equiConst, ...
            'avMolMassGas',avMolMassGas,'avMolMassLiq',avMolMassLiq,...
            'wtFracGas',wtFracGas,'wtFracLiq',wtFracLiq,...
            'd_gasDensity_dz',d_gasDensity_dz,'supVelGas',supVelGas);
        
        d_supVelGas_dz = get_d_supVelGas_dz(parStructSupVelGasDeriv);
%         
       
        %% L        
        
        % gas flux equation
        L1 = diag(ones(N,1));
        L2 = volFracGas*diag(gasDensity)*dispCoefGas*D;
        L3 = zeros(N,N);
        L4 = zeros(N,N);
        
        % gas wt frac equation
        L5 = D;
        L6 = diag(gasDensity)*diag(supVelGas)*D + ...
             diag(gasDensity)*diag(d_supVelGas_dz) + ...
             diag(supVelGas)*diag(d_gasDensity_dz) + ...
             kLa(1)*liqDensity/equiConst(1)*diag(avMolMassGas./avMolMassLiq);
        L7 = zeros(N,N);
        L8 = -kLa(1)*liqDensity*diag(ones(N,1)); %% !! DIGER!
        
        % liquid flux equation
        L9  = zeros(N,N);
        L10 = zeros(N,N);
        L11 = diag(ones(N,1));
        L12 = volFracLiq*liqDensity*dispCoefLiq*D;
        
        % liquid wt frac equation
        L13 = zeros(N,N);
        L14 = -kLa(1)*liqDensity/equiConst(1)*diag(avMolMassGas./avMolMassLiq);
        L15 = D;
        L16 = liqDensity*diag(supVelLiq)*D + ...
              liqDensity*diag(D*supVelLiq) + ...
              kLa(1)*liqDensity*diag(ones(N,1));
          
        L =   [L1 L2 L3 L4; L5 L6 L7 L8; L9 L10 L11 L12; L13 L14 L15 L16];

        
        %% g
        g = zeros(4*N,1);
        parStructReactRate = struct('factor_a',factor_a,'factor_b',factor_b, ...
            'a',a,'b',b,'equiConst',equiConst,'molFracLiq',molFracLiq, ...
            'pTot',pTot);
        
        reactRate    = get_reactRate(parStructReactRate); 
        g(3*N+1:4*N) = volFracLiq.*reactRate.*nu(1).*catDensity.*Mw(1);
      
        %% W
        W = diag([wCOLUMN; wCOLUMN; wCOLUMN; wCOLUMN]);
        
        F = L'*W*g;
        A = L'*W*L;

        B1 = zeros(N,N);
        B2 = zeros(N,N); B2(1,1) = 1; % gas weight fraction of CO set at inlet
        B3 = zeros(N,N);
        B4 = zeros(N,N);
        B5 = zeros(N,N); B5(N,N) = 1; % gas flux of CO set at outlet
        B6 = zeros(N,N); 
        B7 = zeros(N,N);
        B8 = zeros(N,N);
        B9 = zeros(N,N);
        B10 = zeros(N,N);
        B11 = zeros(N,N);
        B12 = zeros(N,N); B12(1,1) = 1; % liquid weight fraction of CO set at inlet
        B13 = zeros(N,N); 
        B14 = zeros(N,N);
        B15 = zeros(N,N); B15(N,N) = 1;
        B16 = zeros(N,N);
        B = [B1 B2 B3 B4; B5 B6 B7 B8; B9 B10 B11 B12; B13 B14 B15 B16];
        
        F_gamma         = zeros(4*N,1);
        F_gamma(1)      = BC(1); % inlet gas weight fraction of CO
        F_gamma(2*N)    = 0;   % outlet gas flux of CO
        F_gamma(2*N+1)  = BC(nCompGas+1); % inlet liquid weight fraction of CO 
        F_gamma(4*N)    = 0;   % outlet liquid flux of CO
        
        %% solve for f
        f_new = (A+B)\(F+F_gamma);
        
        %% use only weight fractions for now - let fluxes go
        disp('CO flux results thrown away.')
        f_new = [f_new(N+1:2*N); f_new(3*N+1:4*N)];
        
        f_new = max(f_new,0);
        
        %% pick up old values for weight fractions
        f_old = [wtFracGas(:,1); wtFracLiq(:,1)];
                     
        %% underrelaxation
        underrelaxation = 0.05;
        f   = f_new*underrelaxation + f_old*(1-underrelaxation);
        
        %% re-shape
        f   = reshape(f,N,2);
        
        %% update solution vector
        solVec(:,1)          = f(:,1); % CO wt fraction gas
        solVec(:,nCompGas+1) = f(:,2); % CO wt fraction liq
        solVec(:,1) = yODE(:,1);
        solVec(:,nCompGas+1) = yODE(:,nCompGas+1);
        
        
        %% calculate L, W, B, A for gas components H2O, C1-C20, C20-C30, C31+
        clear compIndices
        compIndices = [3 4 5 6]; % CO and H2 are omitted; calculated elsewhere
        nCG         = length(compIndices);
                       
        parStructL2wgL3wg = struct('wtFracGas',wtFracGas,'wtFracLiq',wtFracLiq, ...
            'supVelGas',supVelGas, ...
            'gasDensity',gasDensity,'liqDensity',liqDensity, ...
            'kLa',kLa,'equiConst',equiConst, ...
            'avMolMassGas',avMolMassGas,'avMolMassLiq',avMolMassLiq, ...
            'compIndices',compIndices);
        
        L1wg = zeros(N*nCG);
        W    = zeros(N*nCG);
        B    = zeros(N*nCG);
        for gC = 1:nCG
            ind = [(gC-1)*N+1:gC*N];
            L1wg(ind,ind) = D;
            W(ind,ind) = diag(wCOLUMN);
            B(ind(1),ind(1)) = 1;
        end
        
        [L2wg,L3wg] =  get_L2wg_and_L3wg(parStructL2wgL3wg);
                    
        L = L1wg + L2wg + L3wg;

        A = L'*W*L;
        
        %% calculate g
        gMatrix = zeros(N,nCG);
        for zPoint = 1:N
            gMatrix(zPoint,:) = 1./supVelGas(zPoint).*liqDensity./gasDensity(zPoint).* ...
                kLa(compIndices).*wtFracLiq(zPoint,compIndices)';
        end
        g = reshape(gMatrix,N*nCG,1);
        
        %% calculate F
        F = L'*W*g;
        
        %% pick up F_gamma
        F_gamma         = zeros(N,nCG);
        F_gamma(1,:)    = BC(compIndices);
        F_gamma         = reshape(F_gamma,N*nCG,1);
        
        %% solve for f; H2O, C1-C20, C20-C30, C31+
        f_new = (A+B)\(F+F_gamma);
        f_new = max(f_new,0);
        
        f_old = reshape(wtFracGas(:,compIndices),N*nCG,1);
        
        %% underrelaxation
        underrelaxation = 0.05;
        f = f_new*underrelaxation + f_old*(1-underrelaxation);
        
        %% plot weight fractions in gas phase
        figure(1)
        subplot(2,2,1)
        plot(solVec(:,1:nCompGas))
        title('gas weight fractions')
        hold on
        
        %% calculate and evaluate errors and residuals
        res_variable                = L*f - g;
        L2_norm_res_variable        = sqrt((L*f-g)'*W*(L*f-g));
        L2_norm_resBC_variable      = sqrt((B*f-F_gamma)'*W*(B*f-F_gamma));
        L2_norm_tot_var             = L2_norm_res_variable+L2_norm_resBC_variable;
        error_variable              = f-f_old;
        iteration_error_variable    = sqrt((f-f_old)'*W*(f-f_old));
        
        %% plot some diagnostics
        subplot(2,2,2)
        semilogy(iter_var, L2_norm_tot_var,'ko',iter_var,iteration_error_variable,'k*')
        title('L2 residual and iteration error')
        xlabel('iteration_variable number (gass weight fraction variable iteration only)')
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
        
        %% update solVec with four new gas weight fractions
        solVec(:,compIndices) = reshape(f,N,4); % update H2O, C1-C20, C20-C30, C31+ in gas phase
        
        %% calculate the final gas weight fraction (H2)
        wtFracGas = solVec(:,1:nCompGas);
        for zP = 1:N
            wtFracGas(zP,2) = 1-sum(wtFracGas(zP,[1 compIndices]));
        end
        
        solVec(:,2) = wtFracGas(:,2); % update H2 in gas phase
        
        %% --- end of gas weight fractions
        
        %wtFracGas
        pause
        
        %% -- LIQUID FRACTIONS
        %% pack out
        wtFracGas   = solVec(:,1:nCompGas);
%        wtFracLiq   = solVec(:,nCompGas+1:nCompGas+nCompLiq);
        supVelGas   = solVec(:,nCompGas+nCompLiq+1);
        supVelLiq   = solVec(:,nCompGas+nCompLiq+2);
        
        %% CO liquid weight fraction calculated using ADM
        %% C31+ will be calculated at the end
        %% here: components H2, H2O, C1-C20, C20-C30
        clear compIndices
        compIndices = [2 3 4 5];
        nCL = length(compIndices);
        
        %% convert to mole fractions
        molFracGas = zeros(N,nCompGas);
  %      molFracLiq = zeros(N,nCompLiq);
        
        for zPoint = 1:N
            molFracGas(zPoint,:) = wtFracGas(zPoint,:)'./Mw/(sum(wtFracGas(zPoint,:)'./Mw));
%            molFracLiq(zPoint,:) = wtFracLiq(zPoint,:)'./Mw/(sum(wtFracLiq(zPoint,:)'./Mw));
        end
        
        %% update parameters dependent on z or the solution vector
        parStructAvMolMassGas = struct('molFracGas',molFracGas,'Mw',Mw);
 %       parStructAvMolMassLiq = struct('molFracLiq',molFracLiq,'Mw',Mw);
        
        avMolMassGas = get_avMolMassGas(parStructAvMolMassGas);
 %       avMolMassLiq = get_avMolMassLiq(parStructAvMolMassLiq);
        
        %% calculate reaction rate
        parStructReactRate = struct('factor_a',factor_a,'factor_b',factor_b, ...
            'a',a,'b',b,'equiConst',equiConst,'molFracLiq',molFracLiq, ...
            'pTot',pTot);
        
        reactRate = get_reactRate(parStructReactRate);       
        
        parStructWtFracLiqDeriv = struct('wtFracGas',wtFracGas,'wtFracLiq',wtFracLiq, ...
            'supVelLiq',supVelLiq,'liqDensity',liqDensity, ...
            'kLa',kLa,'equiConst',equiConst,'Mw',Mw, ...
            'avMolMassGas',avMolMassGas,'avMolMassLiq',avMolMassLiq, ...
            'nu',nu,'reactRate',reactRate,'catDensity',catDensity, ...
            'compIndices',compIndices);
        
        parStructL2wlL3wl = parStructWtFracLiqDeriv;
         
        %% calculate L, A
        L1wl = zeros(N*nCL);
        W    = zeros(N*nCL);
        B    = zeros(N*nCL);
        for liqComp = 1:nCL
            ind = [(liqComp-1)*N+1:liqComp*N];
            L1wl(ind,ind) = D;
            W(ind,ind) = diag(wCOLUMN);
            B(ind(1),ind(1)) = 1;
        end
        
        [L2wl,L3wl] = get_L2wl_and_L3wl(parStructL2wlL3wl);
        
        L = L1wl + L2wl + L3wl;
        
        A = L'*W*L;
        
        %% calculate g
        gMatrix = zeros(N,nCL);
        for zPoint = 1:N
            gMatrix(zPoint,:) = 1./supVelLiq(zPoint).*kLa(compIndices)*1./equiConst(compIndices).* ...
                avMolMassGas(zPoint)./avMolMassLiq(zPoint).*wtFracGas(zPoint,compIndices)' + ...
                 max(volFracLiq.*nu(compIndices).*reactRate(zPoint)./liqDensity./supVelLiq(zPoint).*catDensity.*Mw(compIndices),0);
        end
        
        g = reshape(gMatrix,N*nCL,1);
        
        %% calculate F
        F = L'*W*g;
        
        %% pick up F_gamma

        F_gamma = zeros(N,nCL);
        F_gamma(1,:) = BC(nCompGas + [compIndices]); %% BC for H2, H2O, C1-C20, C20-C30 wt fraction
        F_gamma = reshape(F_gamma,N*nCL,1);
        
        %% solve for f
        f_new = (A+B)\(F+F_gamma);
        f_new = max(f_new,0);
        
        %% pick up f_old for H2, H2O, C1-C20, C20-C30
        f_old = reshape(wtFracLiq(:,compIndices),N*nCL,1);
        
        %% set underrelaxation here if desired
        underrelaxation = 0.5;
        
        f = f_new*underrelaxation + f_old*(1-underrelaxation);
        
     
        %% plot current f
        figure(2)
        subplot(2,2,1) 
        plot(solVec(:,nCompGas+[1:nCompLiq]))
        title('liquid weight fractions')
        hold on
        
        %% calculate and evaluate errors and residuals
        res_variable                = L*f - g;
        L2_norm_res_variable        = sqrt((L*f-g)'*W*(L*f-g));
        L2_norm_resBC_variable      = sqrt((B*f-F_gamma)'*W*(B*f-F_gamma));
        L2_norm_tot_var             = L2_norm_res_variable+L2_norm_resBC_variable;
        error_variable              = f-f_old;
        iteration_error_variable    = sqrt((f-f_old)'*W*(f-f_old));
        
        %% plot some intermediate diagnostics
        subplot(2,2,2)
        semilogy(iter_var, L2_norm_tot_var,'ko',iter_var,iteration_error_variable,'k*')
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
        
        %% update solution vector with H2, H2O, C1-C20, C20-C30
        solVec(:,nCompGas+[compIndices]) = reshape(f,N,nCL);
        
        
        %% calculate remaining weight fraction: C31+
        wtFracLiq = solVec(:,nCompGas+1:nCompGas+nCompLiq);
        for zP = 1:N
            wtFracLiq(zP,nCompLiq) = 1-sum(wtFracLiq(zP,[1 compIndices]));
        end
        
        solVec(:,nCompGas+nCompLiq) = wtFracLiq(:,nCompLiq); % only C31+ wt fraction changed
        
        if iter_var == max_iter_var
            disp(['Maximal number of iterations (',num2str(max_iter_var),') reached for variable weight fractions'])
        end
        
       % wtFracLiq
      %  pause
        
        
        
    end 
    %% --

    

    disp('***')
    disp('***')
    disp('SUPERFICIAL VELOCITIES')
    disp('***')
    disp('***')
    
    pause
    
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
    F_gamma = zeros(N,2);
    F_gamma(1,:) = BC(nVar-1:nVar);
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
    f_new = max(f_new,0);
    
    %% pick up f_old for plotting
    f_old = [supVelGas; supVelLiq];
    
    %% set underrelaxation here
    underrelaxation = 0.05;
    f = f_new*underrelaxation + f_old*(1-underrelaxation);
    
    %% calculate and evaluate errors and residuals
    res_variable                = L*f - g;
    L2_norm_res_variable        = sqrt((L*f-g)'*W*(L*f-g));
    L2_norm_resBC_variable      = sqrt((B*f-F_gamma)'*W*(B*f-F_gamma));
    L2_norm_tot_var             = L2_norm_res_variable+L2_norm_resBC_variable;
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
    
end

%% unpack and plot results
%% unpack results
wRES_G = solVec(:,1:nCompGas);
wRES_L = solVec(:,nCompGas+1:nCompGas+nCompLiq);
vRES_G = solVec(:,nCompGas+nCompLiq + 1);
vRES_L = solVec(:,nCompGas+nCompLiq + 2);

%% make plots
figure(4)
plot(COLUMN,wRES_G(:,1),...
     COLUMN,wRES_G(:,2),...      
     COLUMN,wRES_G(:,3),...
     COLUMN,wRES_G(:,4),'g--',...
     COLUMN,wRES_G(:,5),'r-.',...
     COLUMN,wRES_G(:,6),'b-' )
title('gas mass fractions')
legend('CO','H_2','H2O','C1-C20','C21-C30','C31+')
ylabel('gas phase mass fractions')
xlabel('reactor length')
grid on

figure(5)
plot(COLUMN,wRES_L(:,1),...
     COLUMN,wRES_L(:,2),...      
     COLUMN,wRES_L(:,3),...
     COLUMN,wRES_L(:,4),'g--',...
     COLUMN,wRES_L(:,5),'r-.',...
     COLUMN,wRES_L(:,6),'b-')
title('liquid mass fractions')
legend('CO','H_2','H2O','C1-C20','C21-C30','C31+')
ylabel('liquid phase mass fractions')
xlabel('reactor length')
grid on

%% calculate conversion of CO (%) on mass basis
wconversion     = 100*wRES_G(:,1);
wcum_conversion = 100*(wtFracGasInit(1) - wRES_G(:,1))./wtFracGasInit(1);

figure(6)
plot(COLUMN,wconversion,COLUMN,wcum_conversion,'r')
title('Conversion on mass basis (%)')
ylabel('conversion (%)')
xlabel('reactor length')
legend('Conversion','Cumulative Conversion')

figure(7)
plot(COLUMN,vRES_G)
title('Superficial gas velocity')

figure(8)
plot(COLUMN,vRES_L)
title('Superficial liquid velocity')
xlabel('reactor length')
ylabel('liquid velocity')
