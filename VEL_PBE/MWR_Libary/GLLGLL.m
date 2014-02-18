% Name      : GLLGLL
%               Grid Generator for the 0D trasient PBE
%               spectral in x
%               spectral in y  
% Author    : Carlos A. Dorao cadorao@nt.ntnu.no   2005
% Location  : <directory in path>/TOOLS
%
% Syntaxis  : struct= GLLGLL(N1, N2, x_min, x_max, y_min, y_max)
% where
%               N1             : Polynomial order
%               N2             : Polynomial order
%               x_min          : Minimum size 
%               x_max          : Maximum size 
%               y_min          : Minimum time  
%               y_max          : Maximum time
%              
% Post      : Return a structure with the grid information 
%

function grid=GridGLLGLL(N1, N2, x_min, x_max, y_min, y_max)

   t0=clock;
  %**********************************************************
  % Algorithm 2.  GRID data structure definition         
  
    N=(N1+1)*(N2+1);
    grid.N=N;
    
    grid.h=(x_max-x_min)*(y_max-y_min)/N^2;
    
    grid.Id=sparse(N,N);
    for i=1:N
        grid.Id(i,i)=1;
    end
    
    grid.xt =zeros(N,2);     %  xi= grid.x(:,1) t= grid.x(:,2) 
    grid.wt =zeros(N,1);     %  omega  
    
    % quadrature xi  xi_\alpha = grid.xi(:,1)
    %                 w_\alpha = grid.xi(:,2)
    grid.x=zeros(N1+1,2);   
    grid.y=zeros(N2+1,2);   
    
    % quadrature t   t_\alpha = grid.t(:,1)  
    %                 w_\t = grid.t(:,2)
    grid.dx=zeros(N1+1);    
    grid.dy=zeros(N2+1);    
%    grid.Dt=sparse(N);% 2D Lag. derivative in (t_min, t_max)


  %**********************************************************
  % Algorithm 3.  Create 2D nodals points         
  
    %grid.x = GLL 
    [grid.x(:,1),grid.x(:,2)]= GLL(N1+1, x_min, x_max); 
    
    %grid.y  = GLL
    [grid.y(:,1),grid.y(:,2)]= GLL(N2+1, y_min, y_max);  

    c=1;
    deltaX=x_max-x_min;
    deltaY=y_max-y_min;
    for tau=1:N2+1
        for eta=1:N1+1
            grid.xt(c,1)=grid.x(eta,1);   % xi 
            grid.xt(c,2)=grid.y(tau,1);    % t
            
            % omega
            grid.wt(c,1)=grid.x(eta,2)*grid.y(tau,2);    
            c=c+1;
        end
    end
 
    grid.W=sparse(N,N);
    for i=1:N
        grid.W(i,i)=grid.wt(i,1);
    end
  %**********************************************************
%  % Algorithm 5. Compute Lagrangean derivative matrix Dt   
    
    % 1D Lag. derivative in (-1,1)
    grid.dx=LagrangeDerivativeMatrix_GLL(N1+1);
    grid.dy=LagrangeDerivativeMatrix_GLL(N2+1);
    grid.Dx=zeros(N,N);
    grid.Dy=zeros(N,N);
    deltaX=x_max-x_min;
    deltaY=y_max-y_min;
    
    %1D Lag. derivative in (t_min, t_max)
    %grid.dt=grid.dt* 2/deltaT;

for q=1:N2+1
    s=q;
    %for s=1:N2+1
          
                for p=1:N1+1
                    for r=1:N1+1
                        
                        i=p+(q-1)*(N1+1);
                        j=r+(s-1)*(N1+1);
                        grid.Dx(i,j)=grid.dx(p,r);%*2/deltaX; %*delta(s,q)
                        %grid.Dy(i,j)=grid.dy(q,s)*delta(r,p)*2/deltaY;
                        %D2(i,j)=D(i2,j2)*delta(i1,j1);
                    end
                end
            
                
         %end
    end

    for q=1:N2+1
         for s=1:N2+1
          
                for p=1:N1+1
                    r=p;
                    %for r=1:N1+1
                        
                        i=p+(q-1)*(N1+1);
                        j=r+(s-1)*(N1+1);
                        %grid.Dx(i,j)=grid.dx(p,r)*2/deltaX; %*delta(s,q)
                        grid.Dy(i,j)=grid.dy(q,s);%*2/deltaY; % *delta(r,p)
                        %D2(i,j)=D(i2,j2)*delta(i1,j1);
                    %end
                end
            
                
         end
    end

    
    
%    figure(60)
%    subplot(1,2,1)
%    spy(grid.Dx)
%    
%    subplot(1,2,2)
%    spy(grid.Dy)


grid.time= etime(clock,t0);

%-----------------------------------------------------------------------------------------
% delta function 
function y=delta (i,j)
    if i==j 
        y=1;
    else
        y=0;
    end
