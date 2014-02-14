%% header
% purpose: solve two-phase axial dispersion model
% for Fischer-Tropsch reactor
% using least-squares method
% camilla.berge.vik@ntnu.no
% 07.02.2014

%% clear memory
clc
clear all
close all

global  gasConst ...
    nCompGas nCompLiq pTot Mw nComp ...
    liqDensityInit nu equiConst dispCoefGas dispCoefLiq 

format long

%% set parameters and initial conditions
parameters

%% set numerical parameters
P = 50;

disp('low number of points - can be increased.')
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
y0 = [wtFracGasInit' wtFracLiqInit' supVelGasInit supVelLiqInit tempGasInit tempSluInit];
[zODE,yODE] = ode15s('model_equations',COLUMN,y0);

%% use ODE estimate as initial guess
solVec = reshape(yODE,N,nVar);

%% boundary conditions
BC = y0';

%% ---- START OF ITERATION LOOP

%% set overall iteration parameters
max_iter            = 50;
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
    max_iter_variable       = 40;
    tol_variable            = 1e-6;
    
    %% make storages
    wtFracGasADM = zeros(N,nComp);
    wtFracLiqADM = zeros(N,nComp);
    
    %% iteration loop for gas weight fractions
    while L2_norm_tot_variable > tol_variable && iter_variable < max_iter_variable
        iter_variable = iter_variable + 1;
        
        %% pack out from solution vector
        wtFracGas   = solVec(:,1:nCompGas);
        wtFracLiq   = solVec(:,nCompGas+1:nCompGas+nCompLiq);
        supVelGas   = solVec(:,nCompGas+nCompLiq+1);
        supVelLiq   = solVec(:,nCompGas+nCompLiq+2);
        tempGas     = solVec(:,nCompGas+nCompLiq+3);
        tempSlu     = solVec(:,nCompGas+nCompLiq+4);
                
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
        
        parStructGasDensity = struct('avMolMassGas',avMolMassGas, ...
            'pTot',pTot,'gasConst',gasConst,'tempGas',tempGas);
        gasDensity = get_gasDensity(parStructGasDensity);
        liqDensity = liqDensityInit; % constant density
               
        %% kinetic constants
        aS = 2.592e-12*exp(-37.3*1000/(gasConst/1000)./tempSlu);
        bS =  1.23e-12*exp( 68.5*1000/(gasConst/1000)./tempSlu);
        
        parStructReactRate = struct('aS',aS,'bS',bS,'equiConst',equiConst, ...
            'molFracLiq',molFracLiq,'pTot',pTot);
        
        reactRate    = get_reactRate(parStructReactRate);
        
        %% mass transfer parameters
        sluDensity  = sluDensityInit;
        heatCapSlu  = heatCapSluInit;
        viscSlu = viscLiq*(1+einsteinK*volFracSol);
        
        parStructkL = struct('liqDensity',liqDensity,'gasDensity',gasDensity, ...
            'diffCoefCO',diffCoefCO,'viscLiq',viscLiq,'graConst',graConst, ...
            'molecularDia',molecularDia,'Mw',Mw,'mWoctacosane',mWoctacosane, ...
            'mDoctacosane',mDoctacosane,'Nparameter',Nparameter, ...
            'tempSlu',tempSlu);
            
        
        kL = getMassTransCoeff(parStructkL);
                
        %% FLUX AND WEIGHT FRACTIONS FOR ALL VARIABLES        
        for cNo = 1:nComp
            %% L
            
            % gas flux equation
            L1 = diag(ones(N,1));
            L2 = volFracGas*diag(gasDensity)*dispCoefGas*D;
            L3 = zeros(N,N);
            L4 = zeros(N,N);
            
            % gas wt frac equation
            L5 = D;
            L6 = diag(supVelGas)*diag(D*gasDensity) + ...
                diag(gasDensity)*diag(supVelGas)*D + ...
                diag(gasDensity)*diag(D*supVelGas) + ...
                diag(kL(:,cNo))*areaDensity*liqDensity/equiConst(cNo)*diag(avMolMassGas./avMolMassLiq);
            L7 = zeros(N,N);
            L8 = -diag(kL(:,cNo))*areaDensity*liqDensity*diag(ones(N,1)); %% !! DIGER!
            
            % liquid flux equation
            L9  = zeros(N,N);
            L10 = zeros(N,N);
            L11 = diag(ones(N,1));
            L12 = volFracLiq*liqDensity*dispCoefLiq*D;
            
            % liquid wt frac equation
            L13 = zeros(N,N);
            L14 = -diag(kL(:,cNo))*areaDensity*liqDensity/equiConst(cNo)*diag(avMolMassGas./avMolMassLiq);
            L15 = D;
            L16 = liqDensity*diag(supVelLiq)*D + ...
                liqDensity*diag(D*supVelLiq) + ...
                diag(kL(:,cNo))*areaDensity*liqDensity*diag(ones(N,1));
            
            L =   [L1 L2 L3 L4; L5 L6 L7 L8; L9 L10 L11 L12; L13 L14 L15 L16];
            
            %% g
            g = zeros(4*N,1);
            
            g(3*N+1:4*N) = volFracLiq.*reactRate.*nu(cNo).*catDensity.*Mw(cNo);
            
            %% W
            W = diag([wCOLUMN; wCOLUMN; wCOLUMN; wCOLUMN]);
            
            %% F
            F = L'*W*g;
            
            %% A
            A = L'*W*L;
            
            %% B
            B1  = zeros(N,N);
            B2  = zeros(N,N); B2(1,1) = 1; % gas weight fraction of CO set at inlet
            B3  = zeros(N,N);
            B4  = zeros(N,N);
            B5  = zeros(N,N); B5(N,N) = 1; % gas flux of CO set at outlet
            B6  = zeros(N,N);
            B7  = zeros(N,N);
            B8  = zeros(N,N);
            B9  = zeros(N,N);
            B10 = zeros(N,N);
            B11 = zeros(N,N);
            B12 = zeros(N,N); B12(1,1) = 1; % liquid weight fraction of CO set at inlet
            B13 = zeros(N,N);
            B14 = zeros(N,N);
            B15 = zeros(N,N); B15(N,N) = 1; % liquid flux of CO set at outlet
            B16 = zeros(N,N);
            B = [B1 B2 B3 B4; B5 B6 B7 B8; B9 B10 B11 B12; B13 B14 B15 B16];
            
            %% F_gamma
            F_gamma         = zeros(4*N,1);
            F_gamma(1)      = BC(cNo);          % inlet gas weight fraction
            F_gamma(2*N)    = 0;                % outlet gas flux
            F_gamma(2*N+1)  = BC(nCompGas+cNo); % inlet liquid weight fraction
            F_gamma(4*N)    = 0;                % outlet liquid flux
            
            %% solve
            f_ADM = (A+B)\(F+F_gamma);
            
            %% assemble old values
            f_oldGAS = wtFracGas(:,cNo);
            f_oldLIQ = wtFracLiq(:,cNo);
        
            %% underrelaxation
            underrelaxation = 0.05;
            
            %% pick out weight fractions and fluxes
            GASFLUX = f_ADM(1:N);
            GASWTFR = max(f_ADM(N+1:2*N),0)*underrelaxation + f_oldGAS*(1-underrelaxation)  ;
            LIQFLUX = f_ADM(2*N+1:3*N);
            LIQWTFR = max(f_ADM(3*N+1:4*N),0)*underrelaxation + f_oldLIQ*(1-underrelaxation);
            
            %% store new values in appropriate matrix
            wtFracGasADM(:,cNo) = GASWTFR;
            wtFracLiqADM(:,cNo) = LIQWTFR;
                
        end
        
        %% update solution vector
        solVec(:,1:nCompGas)                    = wtFracGasADM;
        solVec(:,nCompGas+1:nCompGas+nCompLiq)  = wtFracLiqADM;
        
        %% plot the intermediate results
        figure(300)
        subplot(2,1,1)
        plot(COLUMN,wtFracGasADM)
        ylabel('predicted new gas wt fractions')
        hold on
        subplot(2,1,2)
        plot(COLUMN,wtFracLiqADM)
        ylabel('predicted new liquid wt fractions')
        hold on
        
        
    end 
    
    %% update superficial velocities
    
    %% L
    L1 = diag(gasDensity)*D + diag(D*gasDensity);
    L2 = zeros(N,N);
    L3 = zeros(N,N); 
    L4 = liqDensity*D;
    
    L = [L1 L2; L3 L4];
    
    %% g
    wtFracGas = solVec(:,1:nCompGas);
    wtFracLiq = solVec(:,nCompGas+1:nCompGas+nCompLiq);
    massTransSum = zeros(N,1);
    for zP = 1:N
        massTransSum(zP) = sum(kL(zP,:)'.*areaDensity.*liqDensity.*(1./equiConst.* ...
            avMolMassGas(zP)/avMolMassLiq(zP).*wtFracGas(zP,:)' - ...
           wtFracLiq(zP,:)'));
    end
%     disp('Scaling down mass transfer in velocity predictions')
%     massTransSum = 1e-2.*massTransSum;
    
    g = zeros(2*N,1);
    g(1:N)     = -massTransSum;
    g(N+1:2*N) =  massTransSum;
    
    %% W 
    W = diag([wCOLUMN; wCOLUMN]);
    
    %% F
    F = L'*W*g;
    
    %% A
    A = L'*W*L;
    
    %% B
    B1 = zeros(N,N); B1(1,1) = 1; % set inlet gas superficial velocity
    B2 = zeros(N,N);
    B3 = zeros(N,N); 
    B4 = zeros(N,N); B4(1,1) = 1; % set inlet liquid superficial velocity
    
    B = [B1 B2; B3 B4];
    
    %% F_gamma
    F_gamma      = zeros(2*N,1);
    F_gamma(1)   = BC(nCompGas+nCompLiq+1);
    F_gamma(N+1) = BC(nCompGas+nCompLiq+1);
    
    %% solve
    f_SUPVEL = (A+B)\(F+F_gamma);
    f_newGAS = f_SUPVEL(1:N);
    f_newLIQ = f_SUPVEL(N+1:2*N);
    
    %% pick up old values
    f_oldGAS = supVelGas;
    f_oldLIQ = supVelLiq;
    
    %% plot the intermediate results
    figure(200)
    subplot(2,1,1)
    plot(COLUMN,f_newGAS)
    ylabel('predicted new gas superficial velocity')
    hold on
    subplot(2,1,2)
    plot(COLUMN,f_newLIQ)
    ylabel('predicted new liquid superficial velocity')
    hold on
    
    
    %% underrelaxation
    underrelaxation = 0.5;
    f_GAS = f_newGAS*underrelaxation + f_oldGAS*(1-underrelaxation);
    f_LIQ = f_newLIQ*underrelaxation + f_oldLIQ*(1-underrelaxation);
    
    %% update solution vector
    solVec(:,nCompGas+nCompLiq + 1) = f_GAS;
    solVec(:,nCompGas+nCompLiq + 2) = f_LIQ;
    
    %% update TEMPERATURE for gas and slurry
    %% ---
    disp('temperature')
    
    %% pack out
    %% pack out from solution vector
    wtFracGas   = solVec(:,1:nCompGas);
    wtFracLiq   = solVec(:,nCompGas+1:nCompGas+nCompLiq);
    supVelGas   = solVec(:,nCompGas+nCompLiq+1);
    supVelLiq   = solVec(:,nCompGas+nCompLiq+2);
    tempGas     = solVec(:,nCompGas+nCompLiq+3);
    tempSlu     = solVec(:,nCompGas+nCompLiq+4);
          
    %% heat dispersion parameters
    supVelSlu   = supVelLiq;
   
    lambdaGAS   = dispCoefGas.*gasDensity.*heatCapGas;
    lambdaLIQ   = dispCoefLiq.*sluDensityInit.*heatCapSluInit;
    
    parStructhW = struct('condLiq',condLiq,'condSol',condSol, ...
        'volFracSol',volFracSol,'viscLiq',viscLiq,'graConst',graConst, ...
        'supVelGas',supVelGas,'sluDensity',sluDensity, ...
        'einsteinK',einsteinK,'heatCapSlu',heatCapSlu);
    
    hW = getWallHeatCoeff(parStructhW);

    %% L
    
    % flux equation gas
    L1 = diag(1./volFracGas./lambdaGAS);
    L2 = D;
    L3 = zeros(N,N);
    L4 = zeros(N,N);
    
    % temperature equation gas
    L5 = D*diag(1./gasDensity./heatCapGas./supVelGas);
    L6 = D + diag(hL*areaDensity./gasDensity./heatCapGas./supVelGas);
    L7 = zeros(N,N);
    L8 = -diag(hL*areaDensity./gasDensity./heatCapGas./supVelGas);
    
    % flux equation liquid
    L9  = zeros(N,N);
    L10 = zeros(N,N);
    L11 = diag(ones(N,1)*1./volFracSlu./lambdaLIQ);
    L12 = D;
    
    % temperature equation liquid
    L13 = zeros(N,N);
    L14 = -diag(hL*areaDensity/sluDensityInit/heatCapSluInit./supVelSlu);
    L15 = D*diag(1/sluDensityInit/heatCapSluInit./supVelSlu);
    L16 = D + ...
          diag(hW*perimeter/area/sluDensityInit/heatCapSluInit./supVelSlu)+ ...
          diag(hL*areaDensity/sluDensityInit/heatCapSluInit./supVelSlu);
    
    L = [L1 L2 L3 L4; L5 L6 L7 L8; L9 L10 L11 L12; L13 L14 L15 L16];  
      
    %% g
    g = zeros(4*N,1);
    
    %% ide: putte inn hLa-leddet for gassfasen 
    %% evt. kun den delen av hLa-leddet som er temperaturen i slurryfasen
    
    g(3*N+1:4*N) = volFracLiq*(-deltaHr)*reactRate*catDensity/sluDensityInit/heatCapSluInit./supVelSlu + ...
                   tempSurr*hW*perimeter/area/sluDensityInit/heatCapSluInit./supVelSlu;
    
    %% W
    W = diag([wCOLUMN; wCOLUMN; wCOLUMN; wCOLUMN]);
    
    %% F
    F = L'*W*g;
    
    %% A
    A = L'*W*L;
    
    %% B
    B1  = zeros(N,N); 
    B2  = zeros(N,N); B2(1,1) = 1; % inlet temperature gas
    B3  = zeros(N,N);
    B4  = zeros(N,N);
    B5  = zeros(N,N); B5(N,N) = 1; % outlet heat flux gas
    B6  = zeros(N,N);
    B7  = zeros(N,N);
    B8  = zeros(N,N);
    B9  = zeros(N,N);
    B10 = zeros(N,N);
    B11 = zeros(N,N);
    B12 = zeros(N,N); B12(1,1) = 1; % inlet temperature liquid
    B13 = zeros(N,N);
    B14 = zeros(N,N);
    B15 = zeros(N,N); B15(N,N) = 1; % outlet heat flux liquid
    B16 = zeros(N,N);
    B = [B1 B2 B3 B4; B5 B6 B7 B8; B9 B10 B11 B12; B13 B14 B15 B16];
            
    %% F_gamma
    F_gamma         = zeros(4*N,1);
    F_gamma(1)      = BC(nCompGas+nCompLiq+3);    % inlet gas temperature
    F_gamma(2*N)    = 0;                % outlet gas heat flux
    F_gamma(2*N+1)  = BC(nCompGas+nCompLiq+4); % inlet liquid temperature
    F_gamma(4*N)    = 0;                % outlet liquid heat flux
    
    %% solve
    f = (A+B)\(F+F_gamma);
    f_newGAS = f(N+1:2*N);
    f_newSLU = f(3*N+1:4*N);
    
    %% pick up old values
    f_oldGAS = tempGas;
    f_oldSLU = tempSlu;
    
    %% plot the intermediate results
    figure(100)
    subplot(2,1,1)
    plot(COLUMN,f_newGAS)
    ylabel('gas temperature')
    hold on
    subplot(2,1,2)
    plot(COLUMN,f_newSLU)
    ylabel('liquid temperature')
    hold on
    
    %% underrelaxation
    underrelaxation = 0.5;
    f_GAS = f_newGAS*underrelaxation + f_oldGAS*(1-underrelaxation);
    f_SLU = f_newSLU*underrelaxation + f_oldSLU*(1-underrelaxation);
    
    %% update solution vector
    solVec(:,nCompGas+nCompLiq+3) = f_GAS;
    solVec(:,nCompGas+nCompLiq+4) = f_SLU;
    
    %%
    
end

%% unpack and plot results
%% unpack results
wRES_G = solVec(:,1:nCompGas);
wRES_L = solVec(:,nCompGas+1:nCompGas+nCompLiq);
vRES_G = solVec(:,nCompGas+nCompLiq + 1);
vRES_L = solVec(:,nCompGas+nCompLiq + 2);
tRES_G = solVec(:,nCompGas+nCompLiq + 3);
tRES_L = solVec(:,nCompGas+nCompLiq + 4);

%% calculate mole fractions
mRES_G = zeros(N,nCompGas);
mRES_L = zeros(N,nCompLiq);
for zPoint = 1:N
    mRES_G(zPoint,:) = wRES_G(zPoint,:)'./Mw/(sum(wRES_G(zPoint,:)'./Mw));
    mRES_L(zPoint,:) = wRES_L(zPoint,:)'./Mw/(sum(wRES_L(zPoint,:)'./Mw));
end

%% make plots
figure(1)
subplot(2,1,1)
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

subplot(2,1,2)
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

%% calculate conversion of CO (%) on mole basis
conversion     = 100*mRES_G(:,1);
cum_conversion = 100*(molFracGasInit(1) - mRES_G(:,1))./molFracGasInit(1);

figure(2)
subplot(2,2,1)
plot(COLUMN,conversion,COLUMN,cum_conversion,'r')
ylabel('Conversion (in mole %)')
xlabel('Reactor length')
legend('Conversion','Cumulative Conversion')

subplot(2,2,2)
plot(COLUMN,vRES_G)
xlabel('Reactor length [m]')
ylabel('Superficial gas velocity [m/s]')

subplot(2,2,3)
plot(COLUMN,vRES_L)
xlabel('Reactor length [m]')
ylabel('Superficial liquid velocity [m/s]')

figure(3)
subplot(2,1,1)
plot(COLUMN,tRES_G)
xlabel('Reactor length [m]')
ylabel('Gas Temperature [K]')
subplot(2,1,2)
plot(COLUMN,tRES_L)
xlabel('Reactor length [m]')
ylabel('Slurry Temperature [K]')
