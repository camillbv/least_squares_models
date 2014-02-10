%% header
% purpose: solve gas plug flow model in least squares method
% for a Fischer-Tropsch reactor
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

%% calculate L
L = zeros(N*nVariables);
for varNum = 1:nVariables
    L((varNum-1)*N+1:varNum*N,(varNum-1)*N+1:varNum*N)=D;
end

%% set initial solution vector guess
f_init = zeros(N*nVariables,1);
for varNum = 1:nVariables
    f_init((varNum-1)*N+1:varNum*N) = w_init(varNum)*ones(N,1);
end

%% set weights
for varNum = 1:nVariables
    W((varNum-1)*N+1:varNum*N,(varNum-1)*N+1:varNum*N)=diag(wCOLUMN);
end

%% Calculate A
A = L'*W*L;

%% Calculate B
B           = zeros(nVariables*N);
B1          = zeros(N);
B1(1)       = 1;
for varNum = 1:nVariables
    B((varNum-1)*N+1:varNum*N,(varNum-1)*N+1:varNum*N)=B1;
end

%% Set F_gamma
F_gamma = zeros(nVariables*N,1);
for varNum = 1:nVariables
    dummy = zeros(N,1);
    dummy(1) = w_init(varNum);
    F_gamma((varNum-1)*N+1:varNum*N) = dummy;
end


%% ---- START OF ITERATION LOOP

%% set iteration parameters
max_iter            = 100;
tol                 = 1e-10;
f_old               = f_init; % setting the initial value as the old value
iteration_error     = 1; % initialise
L2_norm_residual    = 1; % initialise
iter                = 0; % initialise
relaxation          = 0.1; 

%% iteration
disp(['Starting ',num2str(max_iter), ' iterations... '])
disp(['Specified tolerance is ',num2str(tol),'.'])
disp('...')
tic
%L2_norm_residual > tol && 

while iter < max_iter % criteria to continue iteration
    %% prepare for a new iteration
    
    %% increase iteration number
    iter                = iter+1;
    
    %% display iteration message
    disp('...')
    disp(['Starting iteration # ',num2str(iter),'...'])
    
    %% re-name weight fractions for convenience
    weightFrac = f_old;
    
    %% convert to mole fractions for use in kinetic expression
    moleFrac = zeros(4*N,1);
    for zPoint = 1:N
        indices = zPoint:N:nVariables*N;
        moleFrac(indices) =  weightFrac(indices)./Mw/(sum(weightFrac(indices)./Mw));
    end
    
    %% calculate partial pressures for use in kinetic expression
    p = moleFrac.*ptot;
    
    %% calculate reaction rate
    reactRate = zeros(N,1);
    for zPoint = 1:N
        if p(zPoint) < 0.1*10^5
            reactRate(zPoint)=0;
        else reactRate(zPoint) = a*factor*factor_a*p(zPoint)*p(N+zPoint)/...
                                     ((1+b*factor_b*p(zPoint))^2);
        end
    end
    
    %% assemble g
    g = zeros(N*nVariables,1);
    for varNum = 1:nVariables
        varIndices = (varNum-1)*N+1:varNum*N;
        g(varIndices)=eps_G*rho_cat/v_GS/rho_G.*nu(varNum)*Mw(varNum).*reactRate;
    end
        
    %% calculate F
    F = L'*W*g;   
    
    %% solve for f_new
    f_new = (A+B)\(F+F_gamma);

    %% under-relaxation may be used here
        
    %% set f value
    f = f_new*relaxation + f_old*(1-relaxation);
    
    %% calculate and store errors
    % iteration error
    iteration_error     = sqrt((f-f_old)'*W*(f-f_old));
    err(iter)           = iteration_error;
    
    % L2 residual
    L2_norm_residual    = sqrt((L*f-g)'*W*(L*f-g));
    res(iter)           = L2_norm_residual;
    
    % L2 residual for boundary conditions
    L2_norm_residualBC  = sqrt((B*f - F_gamma)'*W*(B*f - F_gamma));
    resBC(iter)         = L2_norm_residualBC;
    
    % total L2 residual
    total_L2_norm_residual  =  L2_norm_residual+ L2_norm_residualBC;
    resTOT(iter)            = total_L2_norm_residual;
    
    % iteration number
    iterations(iter)    = iter;
    
    disp(['Iteration summary for iteration # ',num2str(iter),':'])
    disp(['L2-norm residual = ',num2str(L2_norm_residual)])
    disp(['L2-norm residual for boundary conditions = ',num2str(L2_norm_residualBC)])
    disp(['Iteration error = ',num2str(iteration_error)])
    disp('...')
    
    figure(3)
    %% plot error
    plot(f(1:N)-f_old(1:N))
    hold on 
    
    figure(4)
    semilogy(iter,L2_norm_residual,'ko',iter,iteration_error,'k*')
    legend('L2 norm residual','Iteration error')
    hold on
    
    
    %% convergence message to display if convergence is reached
    if iter == max_iter
        disp('...')
        disp('!!')
        disp('NOT CONVERGED - MAXIMUM NUMBER OF ITERATIONS REACHED.')
        disp('!!')
        disp('...')
    end
    
    if L2_norm_residual < tol


        disp('...')
        disp('...')
        disp('CONVERGENCE!')
        disp('...')
        disp('...')
    end
    f_old = f;
end

if L2_norm_residual < tol
    disp(['Overall iteration summary:'])
    disp(['L2-norm residual = ',num2str(L2_norm_residual), ' < tolerance =',num2str(tol),'.'])
    disp(['Convergence reached in ',num2str(iter),' iterations.'])
    disp(['Duration of script execution: ',num2str(toc),' seconds.'])
    
else
    disp('Overall iteration summary:')
    disp(['NOT converged in ',num2str(iter),'iterations.'])
    disp(['L2-norm residual = ',num2str(L2_norm_residual), ' > tolerance =',num2str(tol),'.'])
    disp(['Convergence reached in ',num2str(iter),' iterations.'])
    disp(['Duration of script execution: ',num2str(toc),' seconds.'])
end
%% ---- END OF ITERATION LOOP


%% reshape results vector into results matrix
% nrows: N
% ncols: nVariables
massFractions = zeros(N,nVariables);
for varNum = 1:nVariables
    massFractions(:,varNum)=f((varNum-1)*N+1:varNum*N);
end

%% plot the results
moleFractions = zeros(size(massFractions));

for zstep = 1:length(COLUMN)
    moleFractions(zstep,:) = massFractions(zstep,:)'./Mw/(sum(massFractions(zstep,:)'./Mw));
end


%% print out mole and mass fractions
moleFractions
massFractions

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












