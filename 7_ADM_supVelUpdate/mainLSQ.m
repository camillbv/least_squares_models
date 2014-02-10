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
    nCompGas nCompLiq pTot kLa Mw nComp ...
    liqDensityInit nu equiConst dispCoefGas dispCoefLiq 

format long

%% set parameters and initial conditions
parameters

%% set numerical parameters
P = 40;
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

%% ---- START OF ITERATION LOOP

%% set overall iteration parameters
max_iter            = 10;
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
           
    
        parStructReactRate = struct('factor_a',factor_a,'factor_b',factor_b, ...
            'a',a,'b',b,'equiConst',equiConst,'molFracLiq',molFracLiq, ...
            'pTot',pTot);
        
        reactRate    = get_reactRate(parStructReactRate);
        
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
                kLa(cNo)*liqDensity/equiConst(cNo)*diag(avMolMassGas./avMolMassLiq);
            L7 = zeros(N,N);
            L8 = -kLa(cNo)*liqDensity*diag(ones(N,1)); %% !! DIGER!
            
            % liquid flux equation
            L9  = zeros(N,N);
            L10 = zeros(N,N);
            L11 = diag(ones(N,1));
            L12 = volFracLiq*liqDensity*dispCoefLiq*D;
            
            % liquid wt frac equation
            L13 = zeros(N,N);
            L14 = -kLa(cNo)*liqDensity/equiConst(cNo)*diag(avMolMassGas./avMolMassLiq);
            L15 = D;
            L16 = liqDensity*diag(supVelLiq)*D + ...
                liqDensity*diag(D*supVelLiq) + ...
                kLa(cNo)*liqDensity*diag(ones(N,1));
            
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
            
            %% pick out weight fractions and fluxes
            GASFLUX = f_ADM(1:N);
            GASWTFR = f_ADM(N+1:2*N);
            LIQFLUX = f_ADM(2*N+1:3*N);
            LIQWTFR = f_ADM(3*N+1:4*N);
            
            %% store new values in appropriate matrix
            wtFracGasADM(:,cNo) = GASWTFR;
            wtFracLiqADM(:,cNo) = LIQWTFR;
                
        end
        
%         %% plot the results
%         figure(100)
%         subplot(2,1,1); plot(wtFracGasADM)
%         hold on
%         subplot(2,1,2); plot(wtFracLiqADM)
%         hold on
%         
        %% update solution vector
        solVec(:,1:nCompGas)                    = wtFracGasADM;
        solVec(:,nCompGas+1:nCompGas+nCompLiq)  = wtFracLiqADM;
        
    end 
    
    %% update superficial gas velocity
    
    

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
     COLUMN,wRES_G(:,6),'b-', ...
     COLUMN,wRES_G(:,7),'k:' )
title('gas mass fractions')
legend('CO','H_2','H2O','C1-C10','C11-C20','C21-C30','C31+')
ylabel('gas phase mass fractions')
xlabel('reactor length')
grid on

figure(5)
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
