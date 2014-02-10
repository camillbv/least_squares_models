%% purpose: axial dispersion model for Fischer-Tropsch
%% solved using least squares method
%% weight fractions only

%% camilla.berge.vik@ntnu.no
%% 29.01.2014

%% clear memory
clc
clear all
close all

%% choose long format
format long

%% start figure number counter
figNo = 1;

%% set parameters and initial conditions
parameters

%% calculate points and weights in the reference domain
[z_GLL, wz_GLL]       = GaussLobattoLegendre(N);    % GLL points

%% map points and weights into physical domain
[COLUMN,wCOLUMN] = map(z_GLL,  wz_GLL,  0, totalHeight); % COLUMN length

%% derivative matrix
D_ref  = LagrangeDerivativeMatrix_GLL(N);   % reference domain -1,1
D      = 2/totalHeight*D_ref;               % physical domain   0,totalHeight

%% run ode15s to get initial estimates for weight fractions
y0 = [wtFracGasInit' wtFracLiqInit'];
[zODE,yODE] = ode15s('model_equations',COLUMN,y0);
disp('Cold flow estimate of initial conditions.')

%% THE FOLLOWING VARIABLES ARE CONTAINED IN THE SOLUTION VECTOR

% LSQ
% wtReacG; wtCOgas, wtH2gas
% wtProdG; wtC1gas, wtC2gas, wtC3gas
% wtReacL; wtCOliq, wtH2liq
% wtProdL; wtWAliq, wtC1liq, wtC2liq

% ALGEBRAIC
% wtWateG; wtWAgas
% wtC31pL; wtC3liq

%% number of variables
nVar = (3+nLumps)*2;

%% order of appearance in solution vector:
% SOL(:, 1) = wtCOgas   LSQ
% SOL(:, 2) = wtH2gas   LSQ
% SOL(:, 3) = wtWAgas   ALG
% SOL(:, 4) = wtC1gas   LSQ
% SOL(:, 5) = wtC2gas   LSQ
% SOL(:, 6) = wtC3gas   LSQ

% SOL(:, 7) = wtCOliq   LSQ
% SOL(:, 8) = wtH2liq   LSQ
% SOL(:, 9) = wtWAliq   LSQ
% SOL(:,10) = wtC1liq   LSQ
% SOL(:,11) = wtC2liq   LSQ
% SOL(:,12) = wtC3liq   ALG

%% initial guess of solution vector
SOL = yODE;


%% set those not changing to initial values
liqDensity = liqDensityInit; % constant density
pTot = pTotInit;


%% FIRST: solve for reactants in the liquid phase

%% CO liquid phase loop----------------------------------------------
cNo = 1;

iter_var        = 0;
L2_EQ_var       = 1;
max_iter_var    = 20;
tol_var         = 1e-7;

wtFracGas   = SOL(:,1:nCompGas);

while L2_EQ_var > tol_var && iter_var < max_iter_var
    iter_var = iter_var + 1;
    disp(['   Iteration # ',num2str(iter_var),'...'])
    
    wtFracLiq   = SOL(:,nCompGas+1:nCompGas+nCompLiq);
    f_old       = wtFracLiq(:,cNo); % current CO value 
    
    %% convert to mole fractions
    molFracGas = zeros(N,nCompGas);
    molFracLiq = zeros(N,nCompLiq);
    for zPoint = 1:N
        molFracGas(zPoint,:) = wtFracGas(zPoint,:)'./Mw/(sum(wtFracGas(zPoint,:)'./Mw));
        molFracLiq(zPoint,:) = wtFracLiq(zPoint,:)'./Mw/(sum(wtFracLiq(zPoint,:)'./Mw));
    end
    
    %% update parameters dependent on z or the solution vector
    parAvMolMassGas   = struct('molFracGas',molFracGas,'Mw',Mw);
    parAvMolMassLiq   = struct('molFracLiq',molFracLiq,'Mw',Mw);
    avMolMassGas      = get_avMolMassGas(parAvMolMassGas);
    avMolMassLiq      = get_avMolMassLiq(parAvMolMassLiq);
    parGasDensity     = struct('avMolMassGas',avMolMassGas,'pTot',pTot,...
        'gasConst',gasConst,'T',T);
    gasDensity        = get_gasDensity(parGasDensity);
    parMassTrans      = struct('wtFracGas',wtFracGas,...
        'wtFracLiq',wtFracLiq,'kLa',kLa,'equiConst',equiConst,'Mw',Mw, ...
        'avMolMassGas',avMolMassGas,'avMolMassLiq',avMolMassLiq);
    massTrans         = get_massTrans(parMassTrans);
    
    %% calculate L for CO
    L1 = D;
    L2 = -volFracLiqInit*dispCoef./supVelLiqInit*D.^2;
    L3 = diag(1./supVelLiqInit.*sum(massTrans,2));
    L4 = diag(ones(N,1)*kLa(cNo)/supVelLiqInit);
    L5 = diag(-max(volFracLiqInit.*nu(cNo)*Mw(cNo)*catDensity/liqDensityInit/supVelLiqInit*...
        a*factor_a*pTotInit^2.*wtFracLiq(:,2)*equiConst(1)*equiConst(2)./...
        (1 + b*factor_b*equiConst(1)*wtFracLiq(:,1)*pTotInit).^2,0));
    
    L = L1 + L2 + L3 + L4 + L5;
    
    %% calculate g for CO
    g = kLa(cNo)/supVelLiqInit./equiConst(cNo).*avMolMassGas./avMolMassLiq.*wtFracGas(:,cNo);
    
    %% calculate W for CO
    W = diag(wCOLUMN);
    
    %% calculate F
    F = L'*W*g;
    
    %% calculate A
    A = L'*W*L;
    
    %% calculate F_gamma for CO
    F_gamma = zeros(N,1);
    F_gamma(1) = f_old(1)*supVelLiqInit/dispCoef;
    F_gamma(N) = 0;
    
    %% calculate B for CO
    
    B = zeros(N,N);
    B(1,:) = D(1,:);
    B(N,:) = D(N,:);
    
    f_new = (A+B)\(F+F_gamma);
    
    underrelaxation = 0.2;
    
    f_new = f_new*underrelaxation + f_old*(1-underrelaxation);
    SOL(:,nCompGas+cNo) = f_new;
    
    %% calculate error
    %% plot current f
    figure(figNo)
    
    subplot(2,2,1)
    plot(COLUMN,f_new)
    title('CO liquid weight fractions for each iteration')
    hold on
    
    %% calculate and evaluate errors and residuals
    res_var       = L*f_new - g;
    L2_EQ_var     = sqrt((L*f_new-g)'*W*(L*f_new-g));
    L2_BC_var     = sqrt((B*f_new-F_gamma)'*W*(B*f_new-F_gamma));
    L2_var        = L2_EQ_var+L2_BC_var;
    err_var       = f_new-f_old;
    iter_err_var  = sqrt((f_new-f_old)'*W*(f_new-f_old));
    
    %% plot some intermediate diagnostics
    subplot(2,2,2)
    semilogy(iter_var, L2_EQ_var,'ko',iter_var,iter_err_var,'k*')
    title('L2 residual (equations) and iteration error')
    xlabel('iteration_variable number (liquid weight fraction variable iteration only)')
    legend('L2 norm','iteration error')
    hold on
    
    subplot(2,2,3)
    plot(res_var)
    title('local residual = L*f - g for liquid weight fractions')
    xlabel('z index')
    hold on
    
    subplot(2,2,4)
    plot(err_var)
    title('local iteration error = f-wtFracLiq for liquid weight fractions')
    xlabel('z index')
    hold on
    
        
end
figNo = figNo+1;


%% store final error to use in sum
resCOl = L2_EQ_var;

%% ---------------------------------------------------------------------
cNo = 2;
    
%% H2 liquid phase loop---------------------------------------------------
iter_var        = 0;
L2_EQ_var       = 1;
max_iter_var    = 20;
tol_var         = 1e-7;

while L2_EQ_var > tol_var && iter_var < max_iter_var
    iter_var = iter_var + 1;
    disp(['   Iteration # ',num2str(iter_var),'...'])
    
    wtFracLiq   = SOL(:,nCompGas+1:nCompGas+nCompLiq);
    f_old       = wtFracLiq(:,cNo); % current H2 value 
    
    %% convert to mole fractions
    molFracGas = zeros(N,nCompGas);
    molFracLiq = zeros(N,nCompLiq);
    for zPoint = 1:N
        molFracGas(zPoint,:) = wtFracGas(zPoint,:)'./Mw/(sum(wtFracGas(zPoint,:)'./Mw));
        molFracLiq(zPoint,:) = wtFracLiq(zPoint,:)'./Mw/(sum(wtFracLiq(zPoint,:)'./Mw));
    end
    
    %% update parameters dependent on z or the solution vector
    parAvMolMassGas   = struct('molFracGas',molFracGas,'Mw',Mw);
    parAvMolMassLiq   = struct('molFracLiq',molFracLiq,'Mw',Mw);
    avMolMassGas      = get_avMolMassGas(parAvMolMassGas);
    avMolMassLiq      = get_avMolMassLiq(parAvMolMassLiq);
    parGasDensity     = struct('avMolMassGas',avMolMassGas,'pTot',pTot,...
        'gasConst',gasConst,'T',T);
    gasDensity        = get_gasDensity(parGasDensity);
    parMassTrans      = struct('wtFracGas',wtFracGas,...
        'wtFracLiq',wtFracLiq,'kLa',kLa,'equiConst',equiConst,'Mw',Mw, ...
        'avMolMassGas',avMolMassGas,'avMolMassLiq',avMolMassLiq);
    massTrans         = get_massTrans(parMassTrans);
    
    %% calculate L for H2

    L1 = D;
    L2 = -volFracLiqInit*dispCoef./supVelLiqInit*D.^2;
    L3 = diag(1./supVelLiqInit.*sum(massTrans,2));
    L4 = diag(ones(N,1)*kLa(cNo)/supVelLiqInit);
    L5 = diag(-max(volFracLiqInit.*nu(cNo)*Mw(cNo)*catDensity/liqDensityInit/supVelLiqInit*...
        a*factor_a*pTotInit^2.*wtFracLiq(:,1)*equiConst(1)*equiConst(2)./...
        (1 + b*factor_b*equiConst(1)*wtFracLiq(:,1)*pTotInit).^2,0));
    
    L = L1 + L2 + L3 + L4 + L5;
    
    %% calculate g for H2
    g = kLa(cNo)/supVelLiqInit./equiConst(cNo).*avMolMassGas./avMolMassLiq.*wtFracGas(:,cNo);
    
    %% calculate W for H2
    W = diag(wCOLUMN);
    
    %% calculate F
    F = L'*W*g;
    
    %% calculate A
    A = L'*W*L;
    
    %% calculate F_gamma for H2
    F_gamma = zeros(N,1);
    F_gamma(1) = f_old(1)*supVelLiqInit/dispCoef;
    F_gamma(N) = 0;
    
    %% calculate B for H2
    
    B = zeros(N,N);
    B(1,:) = D(1,:);
    B(N,:) = D(N,:);

    f_new = (A+B)\(F+F_gamma);
    
    underrelaxation = 0.9;
    
    f_new = f_new*underrelaxation + f_old*(1-underrelaxation);
    SOL(:,nCompGas+cNo) = f_new;
    
    %% calculate error
    %% plot current f
    figure(figNo)
    subplot(2,2,1)
    plot(COLUMN,f_new)
    title('H_2 liquid weight fractions for each iteration')
    hold on
    
    %% calculate and evaluate errors and residuals
    res_var       = L*f_new - g;
    L2_EQ_var     = sqrt((L*f_new-g)'*W*(L*f_new-g));
    L2_BC_var     = sqrt((B*f_new-F_gamma)'*W*(B*f_new-F_gamma));
    L2_var        = L2_EQ_var+L2_BC_var;
    err_var       = f_new-f_old;
    iter_err_var  = sqrt((f_new-f_old)'*W*(f_new-f_old));
    
    %% plot some intermediate diagnostics
    subplot(2,2,2)
    semilogy(iter_var, L2_EQ_var,'ko',iter_var,iter_err_var,'k*')
    title('L2 residual (equations) and iteration error')
    xlabel('iteration_variable number (liquid weight fraction variable iteration only)')
    legend('L2 norm','iteration error')
    hold on
    
    subplot(2,2,3)
    plot(res_var)
    title('local residual = L*f - g for liquid weight fractions')
    xlabel('z index')
    hold on
    
    subplot(2,2,4)
    plot(err_var)
    title('local iteration error = f-wtFracLiq for liquid weight fractions')
    xlabel('z index')
    hold on
    
        
end
%% store final error to use in sum
resH2l = L2_EQ_var;

figNo = figNo+1;


%% H2O liquid phase loop
cNo = 3;

iter_var        = 0;
L2_EQ_var       = 1;
max_iter_var    = 20;
tol_var         = 1e-7;

while L2_EQ_var > tol_var && iter_var < max_iter_var
    iter_var = iter_var + 1;
    disp(['   Iteration # ',num2str(iter_var),'...'])
    
    wtFracLiq   = SOL(:,nCompGas+1:nCompGas+nCompLiq);
    f_old       = wtFracLiq(:,cNo); % current H2O value 
    
    %% convert to mole fractions
    molFracGas = zeros(N,nCompGas);
    molFracLiq = zeros(N,nCompLiq);
    for zPoint = 1:N
        molFracGas(zPoint,:) = wtFracGas(zPoint,:)'./Mw/(sum(wtFracGas(zPoint,:)'./Mw));
        molFracLiq(zPoint,:) = wtFracLiq(zPoint,:)'./Mw/(sum(wtFracLiq(zPoint,:)'./Mw));
    end
    
    %% update parameters dependent on z or the solution vector
    parAvMolMassGas   = struct('molFracGas',molFracGas,'Mw',Mw);
    parAvMolMassLiq   = struct('molFracLiq',molFracLiq,'Mw',Mw);
    avMolMassGas      = get_avMolMassGas(parAvMolMassGas);
    avMolMassLiq      = get_avMolMassLiq(parAvMolMassLiq);
    parGasDensity     = struct('avMolMassGas',avMolMassGas,'pTot',pTot,...
        'gasConst',gasConst,'T',T);
    gasDensity        = get_gasDensity(parGasDensity);
    parMassTrans      = struct('wtFracGas',wtFracGas,...
        'wtFracLiq',wtFracLiq,'kLa',kLa,'equiConst',equiConst,'Mw',Mw, ...
        'avMolMassGas',avMolMassGas,'avMolMassLiq',avMolMassLiq);
    massTrans         = get_massTrans(parMassTrans);
    
    %% calculate L for H2
    L1 = D;
    L2 = -volFracLiqInit*dispCoef./supVelLiqInit*D.^2;
    L3 = diag(1./supVelLiqInit.*sum(massTrans,2));
    L4 = diag(ones(N,1)*kLa(cNo)/supVelLiqInit);
    L5 = diag(-max(volFracLiqInit.*nu(cNo)*Mw(cNo)*catDensity/liqDensityInit/supVelLiqInit*...
        a*factor_a*pTotInit^2.*wtFracLiq(:,1)*equiConst(1).*wtFracLiq(:,2).*equiConst(2)./...
        (1 + b*factor_b*equiConst(1)*wtFracLiq(:,1)*pTotInit).^2,0));
    
    L = L1 + L2 + L3 + L4 + L5;
    
    %% calculate g for H2
    g = kLa(cNo)/supVelLiqInit./equiConst(cNo).*avMolMassGas./avMolMassLiq.*wtFracGas(:,cNo);
    
    %% calculate W for H2
    W = diag(wCOLUMN);
    
    %% calculate F
    F = L'*W*g;
    
    %% calculate A
    A = L'*W*L;
    
    %% calculate F_gamma for H2O
    F_gamma = zeros(N,1);
    F_gamma(1) = f_old(1)*supVelLiqInit/dispCoef;
    F_gamma(N) = 0;
    
    %% calculate B for H2O
    
    B = zeros(N,N);
    B(1,:) = D(1,:);
    B(N,:) = D(N,:);

    f_new = (A+B)\(F+F_gamma);
    
    underrelaxation = 0.9;
    
    f_new = f_new*underrelaxation + f_old*(1-underrelaxation);
    SOL(:,nCompGas+cNo) = f_new;
    
    %% calculate error
    %% plot current f
    figure(figNo)
    subplot(2,2,1)
    plot(COLUMN,f_new)
    title('H_2O liquid weight fractions for each iteration')
    hold on
    
    %% calculate and evaluate errors and residuals
    res_var       = L*f_new - g;
    L2_EQ_var     = sqrt((L*f_new-g)'*W*(L*f_new-g));
    L2_BC_var     = sqrt((B*f_new-F_gamma)'*W*(B*f_new-F_gamma));
    L2_var        = L2_EQ_var+L2_BC_var;
    err_var       = f_new-f_old;
    iter_err_var  = sqrt((f_new-f_old)'*W*(f_new-f_old));
    
    %% plot some intermediate diagnostics
    subplot(2,2,2)
    semilogy(iter_var, L2_EQ_var,'ko',iter_var,iter_err_var,'k*')
    title('L2 residual (equations) and iteration error')
    xlabel('iteration_variable number (liquid weight fraction variable iteration only)')
    legend('L2 norm','iteration error')
    hold on
    
    subplot(2,2,3)
    plot(res_var)
    title('local residual = L*f - g for liquid weight fractions')
    xlabel('z index')
    hold on
    
    subplot(2,2,4)
    plot(err_var)
    title('local iteration error = f-wtFracLiq for liquid weight fractions')
    xlabel('z index')
    hold on
    
        
end

resWal = L2_EQ_var;
figNo = figNo + 1;

%% ---------------------------------------------------------------------

%% C1-C20 liquid phase loop---------------------------------------------------
cNo = 4;

iter_var        = 0;
L2_EQ_var       = 1;
max_iter_var    = 20;
tol_var         = 1e-7;

while L2_EQ_var > tol_var && iter_var < max_iter_var
    iter_var = iter_var + 1;
    disp(['   Iteration # ',num2str(iter_var),'...'])
    
    wtFracLiq   = SOL(:,nCompGas+1:nCompGas+nCompLiq);
    f_old       = wtFracLiq(:,cNo); % current C1-C20 value 
    
    %% convert to mole fractions
    molFracGas = zeros(N,nCompGas);
    molFracLiq = zeros(N,nCompLiq);
    for zPoint = 1:N
        molFracGas(zPoint,:) = wtFracGas(zPoint,:)'./Mw/(sum(wtFracGas(zPoint,:)'./Mw));
        molFracLiq(zPoint,:) = wtFracLiq(zPoint,:)'./Mw/(sum(wtFracLiq(zPoint,:)'./Mw));
    end
    
    %% update parameters dependent on z or the solution vector
    parAvMolMassGas   = struct('molFracGas',molFracGas,'Mw',Mw);
    parAvMolMassLiq   = struct('molFracLiq',molFracLiq,'Mw',Mw);
    avMolMassGas      = get_avMolMassGas(parAvMolMassGas);
    avMolMassLiq      = get_avMolMassLiq(parAvMolMassLiq);
    parGasDensity     = struct('avMolMassGas',avMolMassGas,'pTot',pTot,...
        'gasConst',gasConst,'T',T);
    gasDensity        = get_gasDensity(parGasDensity);
    parMassTrans      = struct('wtFracGas',wtFracGas,...
        'wtFracLiq',wtFracLiq,'kLa',kLa,'equiConst',equiConst,'Mw',Mw, ...
        'avMolMassGas',avMolMassGas,'avMolMassLiq',avMolMassLiq);
    massTrans         = get_massTrans(parMassTrans);
    
    %% calculate L for C1-C20

    L1 = D;
    L2 = -volFracLiqInit*dispCoef./supVelLiqInit*D.^2;
    L3 = diag(1./supVelLiqInit.*sum(massTrans,2));
    L4 = diag(ones(N,1)*kLa(cNo)/supVelLiqInit);
    L5 = diag(-max(volFracLiqInit.*nu(cNo)*Mw(cNo)*catDensity/liqDensityInit/supVelLiqInit*...
        a*factor_a*pTotInit^2.*wtFracLiq(:,1)*equiConst(1).*wtFracLiq(:,2)*equiConst(2)./...
        (1 + b*factor_b*equiConst(1)*wtFracLiq(:,1)*pTotInit).^2,0));
    
    L = L1 + L2 + L3 + L4 + L5;
    
    %% calculate g for C1-C20
    g = kLa(cNo)/supVelLiqInit./equiConst(cNo).*avMolMassGas./avMolMassLiq.*wtFracGas(:,cNo);
    
    %% calculate W for C1-C20
    W = diag(wCOLUMN);
    
    %% calculate F
    F = L'*W*g;
    
    %% calculate A
    A = L'*W*L;
    
    %% calculate F_gamma for C1-C20
    F_gamma = zeros(N,1);
    F_gamma(1) = f_old(1)*supVelLiqInit/dispCoef;
    F_gamma(N) = 0;
    
    %% calculate B for C1-C20
    
    B = zeros(N,N);
    B(1,:) = D(1,:);
    B(N,:) = D(N,:);

    f_new = (A+B)\(F+F_gamma);
    
    underrelaxation = 0.9;
    
    f_new = f_new*underrelaxation + f_old*(1-underrelaxation);
    SOL(:,nCompGas+cNo) = f_new;
    
    %% calculate error
    %% plot current f
    figure(figNo)
    subplot(2,2,1)
    plot(COLUMN,f_new)
    title('C1-C20 liquid weight fractions for each iteration')
    hold on
    
    %% calculate and evaluate errors and residuals
    res_var       = L*f_new - g;
    L2_EQ_var     = sqrt((L*f_new-g)'*W*(L*f_new-g));
    L2_BC_var     = sqrt((B*f_new-F_gamma)'*W*(B*f_new-F_gamma));
    L2_var        = L2_EQ_var+L2_BC_var;
    err_var       = f_new-f_old;
    iter_err_var  = sqrt((f_new-f_old)'*W*(f_new-f_old));
    
    %% plot some intermediate diagnostics
    subplot(2,2,2)
    semilogy(iter_var, L2_EQ_var,'ko',iter_var,iter_err_var,'k*')
    title('L2 residual (equations) and iteration error')
    xlabel('iteration_variable number (liquid weight fraction variable iteration only)')
    legend('L2 norm','iteration error')
    hold on
    
    subplot(2,2,3)
    plot(res_var)
    title('local residual = L*f - g for liquid weight fractions')
    xlabel('z index')
    hold on
    
    subplot(2,2,4)
    plot(err_var)
    title('local iteration error = f-wtFracLiq for liquid weight fractions')
    xlabel('z index')
    hold on
            
end

 
%% store final error to use in sum
resC1l = L2_EQ_var;



%% ---------------------------------------------------------------------

%% C20-C30 liquid phase loop
cNo = 5;


%% ---------------------------------------------------------------------
%% algebraic determination of C31+--------------------------------------
for zP = 1:N
    wtFracLiq(zP,nCompLiq)=1-sum(wtFracLiq(zP,1:nCompLiq-1));
end
SOL(:,nCompGas+nCompLiq) = wtFracLiq(:,nCompLiq);
%% --------------------------------------------------------------------- 

%% --------------------------------------------------------------------- 
%% --------------------------------------------------------------------- 
%% GAS MASS FRACTIONS
%% --------------------------------------------------------------------- 
%% --------------------------------------------------------------------- 

wtFracGas   = SOL(:,1:nCompGas);
wtFracLiq   = SOL(:,nCompGas+1:nCompGas+nCompLiq);

%% CO gas phase loop
for cNo = 1:nCompGas
    iter_var        = 0;
    L2_EQ_var       = 1;
    max_iter_var    = 20;
    tol_var         = 1e-7;
    
    figNo = figNo + 1
    
    wtFracGas   = SOL(:,1:nCompGas);
    wtFracLiq   = SOL(:,nCompGas+1:nCompGas+nCompLiq);
    
    while L2_EQ_var > tol_var && iter_var < max_iter_var
        iter_var = iter_var + 1;
        disp(['   Iteration # ',num2str(iter_var),'...'])
        
        wtFracLiq   = SOL(:,nCompGas+1:nCompGas+nCompLiq);
        f_old       = wtFracLiq(:,cNo); % current variable weight fraction
        
        %% convert to mole fractions
        molFracGas = zeros(N,nCompGas);
        molFracLiq = zeros(N,nCompLiq);
        for zPoint = 1:N
            molFracGas(zPoint,:) = wtFracGas(zPoint,:)'./Mw/(sum(wtFracGas(zPoint,:)'./Mw));
            molFracLiq(zPoint,:) = wtFracLiq(zPoint,:)'./Mw/(sum(wtFracLiq(zPoint,:)'./Mw));
        end
        
        %% update parameters dependent on z or the solution vector
        parAvMolMassGas   = struct('molFracGas',molFracGas,'Mw',Mw);
        parAvMolMassLiq   = struct('molFracLiq',molFracLiq,'Mw',Mw);
        avMolMassGas      = get_avMolMassGas(parAvMolMassGas);
        avMolMassLiq      = get_avMolMassLiq(parAvMolMassLiq);
        parGasDensity     = struct('avMolMassGas',avMolMassGas,'pTot',pTot,...
            'gasConst',gasConst,'T',T);
        gasDensity        = get_gasDensity(parGasDensity);
        parMassTrans      = struct('wtFracGas',wtFracGas,...
            'wtFracLiq',wtFracLiq,'kLa',kLa,'equiConst',equiConst,'Mw',Mw, ...
            'avMolMassGas',avMolMassGas,'avMolMassLiq',avMolMassLiq);
        massTrans         = get_massTrans(parMassTrans);
        
        
        
        %% calculate L for each gas phase component
        L1 = D;
               
        L2 = -volFracGasInit*dispCoefGas./supVelGasInit*D.^2;
        L3 = diag(-liqDensityInit./gasDensity./supVelGasInit.*sum(massTrans,2));
        L4 = diag(-liqDensityInit./gasDensity.*kLa(cNo)/supVelGasInit*...
            1./equiConst(cNo).*avMolMassGas./avMolMassLiq);
        
        L = L1 + L2 + L3 + L4;
        
        %% calculate g for each gas phase component
        g = liqDensityInit./gasDensity.*kLa(cNo)/supVelGasInit.*wtFracLiq(:,cNo);
        
        %% calculate W for each gas phase component
        W = diag(wCOLUMN);
        
        %% calculate F for each gas phase component
        F = L'*W*g;
        
        %% calculate A for each gas phase component
        A = L'*W*L;
        
        %% calculate F_gamma for each gas phase component
        F_gamma = zeros(N,1);
        F_gamma(1) = (-gasDensityInit*wtFracGasInit(cNo)+ f_old(1)*gasDensity(1)*supVelGasInit)/(dispCoef*gasDensity(1));
        F_gamma(N) = 0;
        
        %% calculate B for each gas phase component
        
        B = zeros(N,N);
        B(1,:) = D(1,:);
        B(N,:) = D(N,:);
        
        f_new = (A+B)\(F+F_gamma);
        
        SOL(:,cNo) = f_new;
        
        %% calculate error
        %% plot current f
        figure(figNo)
        subplot(2,2,1)
        plot(COLUMN,f_new)
        title(['Component #',num2str(cNo),' gas weight fractions for each iteration'])
        hold on
        
        %% calculate and evaluate errors and residuals
        res_var       = L*f_new - g;
        L2_EQ_var     = sqrt((L*f_new-g)'*W*(L*f_new-g));
        L2_BC_var     = sqrt((B*f_new-F_gamma)'*W*(B*f_new-F_gamma));
        L2_var        = L2_EQ_var+L2_BC_var;
        err_var       = f_new-f_old;
        iter_err_var  = sqrt((f_new-f_old)'*W*(f_new-f_old));
        
        %% plot some intermediate diagnostics
        subplot(2,2,2)
        semilogy(iter_var, L2_EQ_var,'ko',iter_var,iter_err_var,'k*')
        title('L2 residual (equations) and iteration error')
        xlabel('iteration_variable number (liquid weight fraction variable iteration only)')
        legend('L2 norm','iteration error')
        hold on
        
        subplot(2,2,3)
        plot(res_var)
        title('local residual = L*f - g for gas weight fractions')
        xlabel('z index')
        hold on
        
        subplot(2,2,4)
        plot(err_var)
        title('local iteration error = f-wtFracLiq for gas weight fractions')
        xlabel('z index')
        hold on
        
        
    end
    
end






