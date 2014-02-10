%% header
% purpose: solve gas plug flow model in least squares method
% for a Fischer-Tropsch reactor
% using the sequential method
% camilla.berge.vik@ntnu.no
% 21.01.2014

%% clear memory
clc
close all
clear all

%% run file with parameters
parameters

%% set feed gas phase concentrations
x0CO_G     = 0.3;  % kg/kg (unitless) gas feed mole fraction of CO
x0H2_G     = 0.7;  % kg/kg (unitless) gas feed mole fraction of H2
x0H2O_G    = 0;    % kg/kg (unitless) gas feed mole fraction of H2O
x0C16H34_G = 0;    % kg/kg (unitless) gas feed mole fraction of C16H34
x0 = [ x0CO_G  x0H2_G  x0H2O_G  x0C16H34_G]';

%% convert to weight fractions
w0   = x0.*Mw./(x0'*Mw); %  kg/kg (unitless)  gas feed weight fractions
w_init = w0;

%% calculate feed average molar mass
Mw_ave0   = x0'*Mw; % kg/kmol gas feed average molar mass

%% calculate feed gas density
rho_G    = ptot*Mw_ave0/R/T; % kg/m^3  gas feed total density

%% set numerical parameters
P = 100;
N = P+1;

%% number of variables
nVariables = length(w0);

%% calculate points and weights in the reference domain
[z_GLL, wz_GLL]       = GaussLobattoLegendre(N);    % GLL points

%% map points and weights into physical domain
[COLUMN,wCOLUMN] = map(z_GLL,  wz_GLL,  0, L); % COLUMN length

%% derivative matrix
D_ref  = LagrangeDerivativeMatrix_GLL(N);   % reference domain -1,1
D      = 2/L*D_ref;                         % physical domain   0,L

%% calculate sequential variables: L A B W F_gamma f_init
% initialise
globalL = zeros(N,N,nVariables);
globalA = zeros(N,N,nVariables);
globalB = zeros(N,N,nVariables);
globalW = zeros(N,N,nVariables);
globalF_gamma = zeros(N,nVariables);
globalf_init = zeros(N,nVariables);

% set values
for varNum = 1:nVariables
    globalL(:,:,varNum)     = D;
    globalB(1,1,varNum)     = 1;
    globalW(:,:,varNum)     = diag(wCOLUMN);
    globalA(:,:,varNum)     = globalL(:,:,varNum)'*globalW(:,:,varNum)*globalL(:,:,varNum);
    globalf_init(:,varNum)  = w_init(varNum)*ones(N,1);
    globalF_gamma(1,varNum) = w_init(varNum);
end

%% ---- START OF ITERATION LOOP

%% set iteration parameters
max_iter            = 5;
tol                 = 1e-13;
f_old               = globalf_init; % setting the initial value as the old value
iteration_error     = 1; % initialise
L2_norm_tot         = 1; % initialise
iter                = 0; % initialise
relaxation          = 0.1;

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
    
    fB = f_old; % pick up old value to calculate errors at end of each overall iteration
    
    %% -- LOOP EACH VARIABLE IN SEQUENCE
    for varNum = 1:nVariables
        disp('...')
        disp(['Variable # ',num2str(varNum),'...'])
        
        %% re-name for convenience
        weightFracAll = f_old;
        L = globalL(:,:,varNum);
        A = globalA(:,:,varNum);
        B = globalB(:,:,varNum);
        W = globalW(:,:,varNum);
        f_init = globalf_init(:,varNum);
        F_gamma = globalF_gamma(:,varNum);
        
        %% iteration properties for variable iteration
        iter_variable = 0;
        L2_norm_tot_variable = 1;
        max_iter_variable = 100;
        tol_variable = 1e-12;
        
        %% -- VARIABLE ITERATION
        while L2_norm_tot_variable > tol_variable && iter_variable < max_iter_variable
            iter_variable = iter_variable + 1;
  %          disp('   ...')
  %          disp(['   Variable # ',num2str(varNum),', iteration # ',num2str(iter_variable),'...'])
            %% convert to partial pressures for use in kinetic expression
            moleFracAll = zeros(size(weightFracAll));
            partialPAll = zeros(size(weightFracAll));
            for zPoint = 1:N
                moleFracAll(zPoint,:) = weightFracAll(zPoint,:)'./Mw/sum(weightFracAll(zPoint,:)'./Mw);
                partialPAll(zPoint,:) = moleFracAll(zPoint,:).*ptot;
            end
            
            %% calculate reaction rate
            reactRate = zeros(N,1);
            for zPoint = 1:N
                if partialPAll(zPoint,1) < 0.001*10^5 % constraint to avoid negative mass fractions
                    reactRate(zPoint)=0;
                else
                    reactRate(zPoint)= a*factor*factor_a*partialPAll(zPoint,1)*partialPAll(zPoint,2)/...
                        ((1+b*factor_b*partialPAll(zPoint,1))^2);
                end
            end
            
            %% calculate g for this variable only
            g = eps_G*rho_cat/v_GS/rho_G.*nu(varNum)*Mw(varNum).*reactRate;
            
            %% calculate F for this variable only
            F = L'*W*g;
            
            %% solve for f
            f_new = (A+B)\(F+F_gamma);
            
            %% set underrelaxation here if desired
            f = f_new;
            
            %% calculate errors and residuals
         %   disp('   ...')
         %   disp(['   Summary of variable # ',num2str(varNum),', iteration # ',num2str(iter_variable),':']) 
         %   disp('   ...')
            
            res_variable                = L*f - g;
            L2_norm_res_variable        = sqrt((L*f-g)'*W*(L*f-g));
            L2_norm_resBC_variable      = sqrt((B*f-F_gamma)'*W*(B*f-F_gamma));
            L2_norm_tot_variable        = L2_norm_res_variable+L2_norm_resBC_variable;
            iteration_error_variable    = sqrt(f-f_old(:,varNum))'*W*(f-f_old(:,varNum));
            
            f_old(:,varNum) = f; % set new value as old before entering next iteration step
            if iter_variable == max_iter_variable
                disp(['Maximal number of iterations (',num2str(max_iter_variable),') reached for variable # ',num2str(varNum),'.'])
            end 
        end %% end of iteration inside variable 
        globalg(:,varNum)=g;
        globalF(:,varNum)=F;
    end %% end of looping all variables
    fA = f_old; % pick up value of f after iterating over all variables
    
    globalfB = fB;
    globalfA = fA;
    
    %% initialise
    res             = zeros(N,nVariables);
    err             = zeros(N,nVariables);
    L2_norm_res     = zeros(nVariables,1);
    L2_norm_resBC   = zeros(nVariables,1);
    L2_norm_tot     = zeros(nVariables,1);
    iteration_error = zeros(nVariables,1);
    
    %% calculate overall errors and residuals for each overall iteration
    for varNum = 1:nVariables
        
        fA      = globalfA(:,varNum);
        fB      = globalfB(:,varNum);
        L       = globalL(:,:,varNum);
        A       = globalA(:,:,varNum);
        B       = globalB(:,:,varNum);
        W       = globalW(:,:,varNum);
        g       = globalg(:,varNum);
        F       = globalF(:,varNum);
        f_init  = globalf_init(:,varNum);
        F_gamma = globalF_gamma(:,varNum);
        
        res(:,varNum)            = L*fA - g;
        err(:,varNum)            = fA-fB;
        L2_norm_res(varNum)      = sqrt((L*fA-g)'*W*(L*fA-g));
        L2_norm_resBC(varNum)    = sqrt((B*fA-F_gamma)'*W*(B*fA-F_gamma));
        L2_norm_tot(varNum)      = L2_norm_res(varNum)+L2_norm_resBC(varNum);
        iteration_error(varNum)  = sqrt((fA-fB)'*W*(fA-fB));
    end
    
    %% re-shape for plotting
    resPlot = zeros(N*nVariables);
    errPlot = zeros(N*nVariables);
    for varNum = 1:nVariables
        resPlot((varNum-1)*N+1:varNum*N) = res(:,varNum);
        errPlot((varNum-1)*N+1:varNum*N) = err(:,varNum);
    end
    
    %% plot iteration errors during evaluation
    figure(3)
    
    subplot(1,4,1)
    plot(err(:,1))
    title('f-f_old for weight fraction of CO')
    xlabel('z-point index')
    hold on
    subplot(1,4,2)
    plot(err(:,2))
    title('f-f_old for weight fraction of H2')
    xlabel('z-point index')
    hold on
    subplot(1,4,3)
    plot(err(:,3))
    title('f-f_old for weight fraction of H2O')
    xlabel('z-point index')
    hold on
    subplot(1,4,4)
    plot(err(:,4))
    title('f-f_old for weight fraction of C16H34')
    xlabel('z-point index')
    hold on
    
    figure(4)
    semilogy(iter,sum(L2_norm_tot),'ko',iter,sum(iteration_error),'k*')
    title('Total L2 norm residual and iteration error')
    legend('L2 norm residual','Iteration error')
    xlabel('Iteration # (overall iteration)')
    hold on
    
    
end

  

massFractions = f_old;
moleFractions = zeros(size(massFractions));
for zPoint = 1:N
    moleFractions(zPoint,:) = massFractions(zPoint,:)'./Mw/sum(massFractions(zPoint,:)'./Mw);
end

%% calculate conversion of CO (%)
conversion     = 100*moleFractions(:,1);
cum_conversion = 100*(x0CO_G - moleFractions(:,1))./x0CO_G;

%% check the sum of mass fractions
massFractionsIn = sum(massFractions(1,:))

dummy = size(massFractions);
nRows = dummy(1);

massFractionsOut = sum(massFractions(nRows,:))

%% plot the results
figure(1)

subplot(1,2,1)
plot(COLUMN,moleFractions)
xlabel('Reactor length [m]')
ylabel('Gas phase mole fractions (mol/mol)')
legend('CO','H_2','H_2O','C_{16}H_{34}')

subplot(1,2,2)
plot(COLUMN,conversion,'k--',COLUMN,cum_conversion,'r--')
xlabel('Reactor length [m]')
ylabel('Conversion (%)')
legend('Local Conversion of CO','Cumulative Conversion of CO')

