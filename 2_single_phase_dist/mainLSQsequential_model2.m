%% header
% purpose: solve gas plug flow model in least squares method
% using sequential (equation-by-equation) method
% with full product distribution
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
x0CO_G     = 0.3;  % mole/mole (unitless) gas feed mole fraction of CO     
x0H2_G     = 0.7;  % mole/mole (unitless) gas feed mole fraction of H2    
x0H2O_G    = 0;    % mole/mole (unitless) gas feed mole fraction of H2O    

x0 = [ x0CO_G  x0H2_G  x0H2O_G  zeros(Cmax+1,1)']';

%% convert to weight fractions
w0   = x0.*Mw./(x0'*Mw); %  kg/kg (unitless)  gas feed weight fractions
w_init = w0;

%% calculate feed average molar mass
Mw_ave0   = x0'*Mw; % kg/kmol gas feed average molar mass

%% calculate feed gas density
rho_G    = ptot*Mw_ave0/R/T; % kg/m^3  gas feed total density

%% set numerical parameters
P = 50;
N = P+1;

%% number of variables
nVariables = length(x0);

%% calculate points and weights in the reference domain
[z_GLL, wz_GLL]       = GaussLobattoLegendre(N);    % GLL points

%% map points and weights into physical domain
[COLUMN,wCOLUMN] = map(z_GLL,  wz_GLL,  0, L); % COLUMN length

%% derivative matrix
D_ref  = LagrangeDerivativeMatrix_GLL(N);   % reference domain -1,1      
D      = 2/L*D_ref;                         % physical domain   0,L

%% calculate each L
L = zeros(N,N,nVariables);
for varNum = 1:nVariables
    L(:,:,varNum)=D;
end

%% calculate each f_init
for varNum = 1:nVariables
    f_init(:,varNum) = w_init(varNum)*ones(N,1);
end

%% calculate each set of weights
for varNum = 1:nVariables
    W(:,:,varNum)=diag(wCOLUMN);
end

%% calculate each set of A
for varNum = 1:nVariables
    A(:,:,varNum)=L(:,:,varNum)'*W(:,:,varNum)*L(:,:,varNum);
end

%% set each B
B1          = zeros(N);
B1(1,1)       = 1;
for varNum = 1:nVariables
    B(:,:,varNum)=B1;
end

%% set each F_gamma
for varNum = 1:nVariables
    dummy = zeros(N,1);
    dummy(1) = w_init(varNum);
    F_gamma(:,varNum) = dummy;
end

%% calculate constant part of g
constantGTerm = eps_G*rho_cat/v_GS/rho_G*a*factor*factor_a*ptot^2/(Mw(1)*Mw(2));

%% ---- START OF ITERATION LOOP

%% set iteration parameters
max_iter            = 50;
overall_tol                 = 1e-6;
f_old               = f_init; % setting the initial value as the old value
total_iteration_error     = 1; % initialise
total_L2_norm_residual    = 1; % initialise
iteration_error = 1;
iter                = 0; % initialise
relaxation          = 0.5; 

%% iteration
disp(['Starting ',num2str(max_iter), ' iterations... '])
disp(['Specified tolerance is ',num2str(overall_tol),'.'])
disp('...')
tic

while total_L2_norm_residual > overall_tol && iter < max_iter % criteria to continue iteration
    %% prepare for a new iteration
    iter                = iter+1;
    disp('...')
    disp(['Starting iteration # ',num2str(iter),' ...'])
    
    weightFrac = f_old; % for readability
    variable_tol = 1e-10;
    
    for varNum = 1:nVariables  % sequential solution    
        disp(['Solving for variable #',num2str(varNum),' of 34 ...'])
        
        
        while iteration_error > variable_tol % iterate on each variable
            %% calculate each g
            for pNum = 1:N
                sumFracMolWeight    = sum(weightFrac(pNum,:)'./Mw);  % calculate summation term
                wCO = weightFrac(pNum,1); % weight fraction CO
                wH2 = weightFrac(pNum,2); % weight fraction H2
                % calculate point dependent term
                pointDependentGTerm(pNum) = wCO*wH2/(sumFracMolWeight^2)/...
                    ((1+b*factor_b*ptot*wCO/Mw(1)/sumFracMolWeight)^2);
            end
            
            variableDependentGTerm = nu(varNum)*Mw(varNum);
            g(:,varNum) = constantGTerm*variableDependentGTerm*pointDependentGTerm;
            
            %% calculate each F
            F(:,varNum) = L(:,:,varNum)'*W(:,:,varNum)*g(:,varNum);
            
            %% solve each f
            f = (A(:,:,varNum) + B(:,:,varNum))\(F(:,varNum)+F_gamma(:,varNum));
            
            %% calculate error for each variable
            L2_norm_residual    = sqrt( (L(:,:,varNum)*f-g(:,varNum))'*W(:,:,varNum)*...
                (L(:,:,varNum)*f-g(:,varNum)) );
            L2_norm_residualBC  = sqrt( (B(:,:,varNum)*f-F_gamma(:,varNum))'*W(:,:,varNum)*...
                (B(:,:,varNum)*f-F_gamma(:,varNum)) );
            iteration_error     = sqrt((f-f_old(:,varNum))'*W(:,:,varNum)*(f-f_old(:,varNum)));
            
            total_L2_norm_residual =  L2_norm_residual+ L2_norm_residualBC;
            
            %% store error information for each variable and iteration
            err(iter,varNum)           = iteration_error;
            res(iter,varNum)           = L2_norm_residual;
            resBC(iter,varNum)         = L2_norm_residualBC;
            resTOT(iter,varNum)        = total_L2_norm_residual
            
            %% update each f_old
            f_old(:,varNum) = f;

          
        end  
        
        %% -- GO TO NEXT VARIABLE 
        iteration_error = 1;  % re-set iteration error
    end % loop the sequence of variables
    
    total_L2_norm_residual = sum(resTOT(iter,:));
    total_iteration_error = sum(err(iter,:));
           
    disp(['Iteration summary for iteration # ',num2str(iter),':'])
    disp(['Total L2-norm residual = ',num2str(total_L2_norm_residual)])
    disp(['Iteration error = ',num2str(total_iteration_error)])
    disp('...')
    
    %% convergence message to display if convergence is reached
    if iter == max_iter
        disp('...')
        disp('!!')
        disp('Not converged - maximum iterations reached.')
        disp('')
        disp('Please re-run with more points in spatial or phase space.')
        disp('!!')
        disp('...')
    end
    
    if L2_norm_residual < overall_tol

        disp('...')
        disp('...')
        disp('CONVERGENCE!')
        disp('...')
        disp('...')
    end
end % end while

%% plot the results
weightFractions = f_old;

%% plot the results
moleFractions = zeros(size(weightFractions));

for zstep = 1:length(COLUMN)
    moleFractions(zstep,:) = weightFractions(zstep,:)'./Mw/(sum(weightFractions(zstep,:)'./Mw));
end

%% calculate conversion of CO (%)
conversion     = 100*moleFractions(:,1);
cum_conversion = 100*(x0CO_G - moleFractions(:,1))./x0CO_G;

%% plot the results
figure(1)

subplot(1,3,1)
plot(COLUMN,moleFractions)
xlabel('Reactor length [m]')
ylabel('Mole fractions (mol/mol)')
legend('C1','C2','C3',...
    'C4','C5','C6','C7','C8','C9','C10',...
    'C11','C12','C13','C14','C15','C16',...
    'C17','C18','C19','C20','C21','C22',...
    'C23','C24','C25','C26','C27','C28',...
    'C29','C30','C31+')

subplot(1,3,2)
plot(COLUMN,weightFractions)
xlabel('Reactor length [m]')
ylabel('Mole fractions (mol/mol)')
legend('C1','C2','C3',...
    'C4','C5','C6','C7','C8','C9','C10',...
    'C11','C12','C13','C14','C15','C16',...
    'C17','C18','C19','C20','C21','C22',...
    'C23','C24','C25','C26','C27','C28',...
    'C29','C30','C31+')

subplot(1,3,3)
plot(COLUMN,conversion,'k--',COLUMN,cum_conversion,'r--')
xlabel('Reactor length [m]')
ylabel('Conversion (%)')
legend('Local Conversion of CO','Cumulative Conversion of CO')


%% check the sum of mass fractions
massFractionsIn = sum(weightFractions(1,:))

dummy = size(weightFractions);
nRows = dummy(1);

massFractionsOut = sum(weightFractions(nRows,:))