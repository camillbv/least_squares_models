% Name      :   cong_gradient 
% Author    :   Carlos A. Dorao cadorao@nt.ntnu.no      2005
% Location  :   <directory in path>/TOOLS
%
% Purpose   :   Solve the system A x = F using a congujate gradient solver
% Systaxis  :   [U,iteration,delta0]=cong_gradient(A,F,epsilon)
%                   A = SPD matrix N x N 
%                   F = vector N x 1
%               return
%                   U = solution N x 1
%                   iteration = number of iteration used
%                   delta0    = accuracy 

function[U,iteration,delta0]=cong_gradient(A,F,epsilon)

    iteration=1;
    %U=sparse(length(F),1);
    U=zeros(length(F),1);
    R=F;
    R=A*U-R;
    delta0=R'*R;
    if (delta0>epsilon) 
        D=-R;
        while(1)
            iteration=iteration+1;
            B=A*D;
            tau=delta0/(D'*B);
            U=U+tau*D;
            R=R+tau*B;
            delta1=R'*R;
            if delta1 <= epsilon
                clear R;
                clear D;
                clear B;
                break;    
            end
            beta=delta1/delta0;
            delta0=delta1;
            D=-R+beta*D;
        end
        
    else
    end
    
    
  