% Name      : GridGLGLL
%               Grid Generator for the 0D trasient PBE
%               spectral in xi
%               spectral in t  
% Author    : Carlos A. Dorao cadorao@nt.ntnu.no   2005
% Location  : <directory in path>/TOOLS
%
% Syntaxis  : struct= LSQ_0D_St_source( N1, ...
%                                       N2, ...
%                                       xi_min, ...
%                                       xi_max, ...
%                                       t_min, ...
%                                       t_max)
% where
%               N1             : Polynomial order
%               N2             : Polynomial order
%               xi_min         : Minimum size 
%               xi_max         : Maximum size 
%               t_min          : Minimum time  
%               t_max          : Maximum time
%              
% Post      : Return a structure with the grid information 
%

function grid=GridGLGLL_1DSt(N1, N2, xi_min, xi_max, t_min, t_max)

   t0=clock;
  %**********************************************************
  % Algorithm 2.  GRID data structure definition         
  
    N=(N1+1)*(N2+1);
    
    grid.x =zeros(N,2);     %  xi= grid.x(:,1) t= grid.x(:,2) 
    grid.w =zeros(N,1);     %  omega  
    
    % quadrature xi  xi_\alpha = grid.xi(:,1)
    %                 w_\alpha = grid.xi(:,2)
    grid.xi=zeros(N1+1,2);   
    
    % quadrature t   t_\alpha = grid.t(:,1)  
    %                 w_\t = grid.t(:,2)
    grid.t=zeros(N2+1,2);   
    
    % 1D Lag. deriv. mat. in (-1, 1)
    grid.dt=zeros(N2+1);    
%    grid.Dt=sparse(N);% 2D Lag. derivative in (t_min, t_max)

    grid.Bxi=zeros(N1+1);   % Matrix Bxi(eta,alpha) 
    
    % Matrix LBxi(beta,eta,alpha)
    grid.BLxi=zeros(N1+1,N1+1,N1+1);    
    
    grid.Cdxi=zeros(N1+1);   % Matrix death (eta,alpha) 

    grid.Cxi=zeros(N1+1);   % Matrix Bxi(eta,alpha) 
    
    % Matrix CLxi(beta,eta,alpha)
    grid.CLxi=zeros(N1+1,N1+1,N1+1);    
    grid.CLFxi=zeros(N1+1,N1+1,N1+1);   
    grid.CdLFxi=zeros(N1+1,N1+1,N1+1);   
    grid.CbLFxi=zeros(N1+1,N1+1,N1+1);   
    grid.CbLF2xi=zeros(N1+1,N1+1,N1+1);   
    
  %**********************************************************
  % Algorithm 3.  Create 2D nodals points         

%   fprintf('t_min= %f  t_max= %f \n ', t_min, t_max); 
%   GLL(N2+1, t_min, t_max)

   %break 
   
  
    %grid.xi = GL 
    [grid.xi(:,1),grid.xi(:,2)]= GL(N1+1, xi_min, xi_max); 
    
    %grid.t  = GLL
    [grid.t(:,1),grid.t(:,2)]= GLL(N2+1, t_min, t_max);  

    c=1;
    deltaXi=xi_max-xi_min;
    deltaT=t_max-t_min;
    for tau=1:N2+1
        for eta=1:N1+1
            grid.x(c,1)=grid.xi(eta,1);   % xi 
            grid.x(c,2)=grid.t(tau,1);    % t
            
            % omega
            grid.w(c,1)=grid.xi(eta,2)*grid.t(tau,2);    
            c=c+1;
        end
    end
 
  %**********************************************************
  % Algorithm 4.  Quadrature rule for breakage birth         
   
     for alpha=1:N1+1
        % algorithm -----------------
        %aux= GL(N1+1, grid.xi(alpha,1) , xi_max); 
        %for k=1:N1+1
        %    grid.BXi(alpha,k)=aux();   % xi 
        %end
        % mathlab short form 
        aux=GL(N1+1, grid.xi(alpha,1) , xi_max); % GL
        grid.Bxi(alpha,:)=aux';
    end

     
  %**********************************************************
  % Algorithm 5. Lagrangean Polynomial for breakage birth      
  deltaXi=xi_max-xi_min;
  for alpha=1:N1+1
      for beta=1:N1+1
         for k=1:N1+1
             % mapping (xi_min, xi_max) -> (-1,1)  where is 
             % given the Lagrangean polynomial
             pt= 2/deltaXi*(grid.Bxi(alpha,k)-xi_min ) -1;
             
             grid.BLxi(alpha,beta, k)=LagrangeGL(beta, pt, N1+1);        
        end
    end
 end

 
  %**********************************************************
  % Algorithm .  Quadrature rule for Coalesce birth         
   
     for alpha=1:N1+1
        aux=GL(N1+1,xi_min, grid.xi(alpha,1) ); % GL
        grid.Cxi(alpha,:)=aux';
    end
 
%    grid.Cxi
  %**********************************************************
  % Algorithm . Lagrangean Polynomial for coalescence birth      
  deltaXi=xi_max-xi_min;
  for alpha=1:N1+1
%      delta=(xi_max+xi_min-grid.xi(alpha,1))/(xi_max-xi_min);
      delta=(grid.xi(alpha,1)-xi_min)/(xi_max-xi_min);
      for beta=1:N1+1
         for k=1:N1+1
             % mapping (xi_min, xi_max) -> (-1,1)  where is 
             % given the Lagrangean polynomial
             ptxi= 2/deltaXi*(grid.xi(alpha,1)-xi_min ) -1;
             pt= 2/deltaXi*(grid.Cxi(alpha,k)-xi_min ) -1;
      
             
             grid.CLxi(alpha,beta, k)=LagrangeGL(beta, ptxi- pt, N1+1)*grid.xi(k,2)*delta;        
        end
    end
 end

  %**********************************************************
  % Algorithm . Lagrangean Polynomial for coalescence birth F*     
  deltaXi=xi_max-xi_min;
  for alpha=1:N1+1
      for beta=1:N1+1
         for k=1:N1+1
             % mapping (xi_min, xi_max) -> (-1,1)  where is 
             % given the Lagrangean polynomial
             pt= 2/deltaXi*(grid.Cxi(alpha,k)-xi_min ) -1;
             
             grid.CLFxi(alpha,beta, k)=LagrangeGL(beta, pt, N1+1);        
        end
    end
 end

 %***********************************************************
 
  %**********************************************************
  % Algorithm .  Quadrature rule for Coalesce death          
   
     for alpha=1:N1+1
        aux=GL(N1+1,xi_min,  xi_max+xi_min-grid.xi(alpha,1) ); % GL
        grid.Cdxi(alpha,:)=aux';
    end

    %grid.xi
    %grid.Cdxi
    
  %**********************************************************
  % Algorithm . Lagrangean Polynomial for coalescence death F*     
  deltaXi=xi_max-xi_min;
  for alpha=1:N1+1
      for beta=1:N1+1
         for k=1:N1+1
             % mapping (xi_min, xi_max) -> (-1,1)  where is 
             % given the Lagrangean polynomial
             pt= 2/deltaXi*(grid.Cdxi(alpha,k)-xi_min ) -1;
             
             grid.CdLFxi(alpha,beta, k)=LagrangeGL(beta, pt, N1+1);        
        end
    end
 end


  % BIRTH *************************************************** 
  %**********************************************************
  % Algorithm . Lagrangean Polynomial for coalescence death F*     
  deltaXi=xi_max-xi_min;
  for alpha=1:N1+1
      for beta=1:N1+1
         for k=1:N1+1
             % mapping (xi_min, xi_max) -> (-1,1)  where is 
             % given the Lagrangean polynomial
             pt= 2/deltaXi*(grid.Cxi(alpha,k)-xi_min ) -1;
             
             grid.CbLFxi(alpha,beta, k)=LagrangeGL(beta, pt, N1+1);        
        end
    end
 end

 
 %**********************************************************
  % Algorithm . Lagrangean Polynomial for coalescence death F*     
  deltaXi=xi_max-xi_min;
  for alpha=1:N1+1
      for beta=1:N1+1
         for k=1:N1+1
             % mapping (xi_min, xi_max) -> (-1,1)  where is 
             % given the Lagrangean polynomial
             pt  = 2/deltaXi*(grid.Cxi(alpha,k)-xi_min ) -1;
             ptxi= 2/deltaXi*(grid.xi(alpha,1)-xi_min ) -1;
             deltapt  = 2/deltaXi*(grid.xi(alpha,1)-grid.Cxi(alpha,k)-xi_min ) -1;
           
             %grid.CbLF2xi(alpha,beta, k)=LagrangeGL(beta, ptxi-pt, N1+1);        
             grid.CbLF2xi(alpha,beta, k)=LagrangeGL(beta, deltapt, N1+1);        
        end
    end
 end

 
 
 
 
 
 
 
 
 
    
  %**********************************************************
%  % Algorithm 5. Compute Lagrangean derivative matrix Dt   
    
    % 1D Lag. derivative in (-1,1)
    grid.dt=LagrangeDerivativeMatrix_GLL(N2+1);
    deltaT=t_max-t_min;
    
%    grid.dt 
    
    %1D Lag. derivative in (t_min, t_max)
    grid.dt=grid.dt* 2/deltaT;

 %   grid.dt
    
%    for gamma=1:N2+1
%        for eta=1:N1+1
%            for tau=1:N2+1
%                i=eta + tau *(N1+1);
%                j=eta + gamma *(N1+1);
%                grid.Dt(i,j)=grid.dt(tau,gamma) * 2/deltaT;
%            end
%        end
%    end

grid.time= etime(clock,t0);

