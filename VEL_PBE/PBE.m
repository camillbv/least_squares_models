function [fd] = PBE(vZ,fdOLD,RHOd)


global  nZ nXI nGLOB GM derCOLUMN derBUBBLE COLUMN LAMBDAglob B
global  BUBBLEmax BUBBLEmin BUBBLE
global  xGLL wGLL xGL wGL
global  Crate_Cd Crate_Cb LagGLL_Cd LagGLL_Cb LagGLL_Cb2 
global  L_B fdBC fdIN
global  LagXI LagZ
global  ITERmax ITERtol


%==========================================================================
%L-MATRIX - CONVECTION IN PHYSICAL SPACE
%==========================================================================
dervd = zeros(nZ,nXI);
for iXI=1:nXI
 	dervd(:,iXI) = derCOLUMN*vZ(:,iXI);
end


L_Con = zeros(nGLOB,nGLOB);

for jZ = 1:nZ
    for jXI = 1:nXI
        jGLOB = GM(jZ,jXI);
            for iqZ = 1:nZ
                for iqXI = 1:nXI
                    iqGLOB = GM(iqZ,iqXI);
                
                        L1 = derCOLUMN(iqZ,jZ)*LagXI(iqXI,jXI)*vZ(iqZ,iqXI);
                        L2 = LagZ(iqZ,jZ)*LagXI(iqXI,jXI)*dervd(iqZ,iqXI);
                        L_Con(iqGLOB,jGLOB) = L_Con(iqGLOB,jGLOB) + L1 + L2;
                
                end%END iqXI
            end%END iqZ
    end%END jXI
end%END jZ


%==========================================================================
%L-MATRIX - GROWTH IN PROPERTY SPACE
%==========================================================================
derRHOd = derCOLUMN*RHOd;

%PRECALCULATION - GROWTH TERM
L_Growth = zeros(nGLOB,nGLOB);
for jZ = 1:nZ
    iqZ = jZ;
        for jXI = 1:nXI
            jGLOB = GM(jZ,jXI);
                for iqXI = 1:nXI
                    iqGLOB = GM(iqZ,iqXI);
                                
                        L1 = -BUBBLE(iqXI)/(3*RHOd(iqZ))*derRHOd(iqZ)*derBUBBLE(iqXI,jXI)*LagZ(iqZ,jZ)*vZ(iqZ,iqXI);
                        L2 = -1/(3*RHOd(iqZ))*derRHOd(iqZ)*LagXI(iqXI,jXI)*LagZ(iqZ,jZ)*vZ(iqZ,iqXI);
                        L3 = -1/(3*RHOd(iqZ))*derRHOd(iqZ)*LagXI(iqXI,jXI)*LagZ(iqZ,jZ)*BUBBLE(iqXI)*(derBUBBLE(iqXI,:)*vZ(iqZ,:)'); 
                        
                        L_Growth(iqGLOB,jGLOB) = L_Growth(iqGLOB,jGLOB) + L1 + L2 + L3;
                        
                end%END iqXI
        end%END jXI
end%END jZ


%==========================================================================
%PRECALCULATED L-MATRIX AHEAD OF ITERATION LOOP 
%==========================================================================
L_PRE = L_Con + L_B + L_Growth;  



%**************************************************************************
%**************************************************************************
%ITERATION LOOP
%**************************************************************************
%**************************************************************************

g        = zeros(nGLOB,1);
fdMATRIX = zeros(nZ,nXI);


ITER     = 0;
RELAX    = 0.7;
ITERerr  = 1;

%*********************************************
while  ITER < ITERmax  &&  ITERerr > ITERtol 
%*********************************************
 
ITER = ITER + 1; 
 


%==========================================================================
%L-MATRIX - COALESCENCE DEATH
%==========================================================================
L_Cd = zeros(nGLOB,nGLOB);

for iqZ = 1:nZ
    jZ = iqZ;

    C    = zeros(nXI,nXI);
    for iq = 1:nXI
        
        iqGLOB = GM(iqZ,iq);
        
        if BUBBLE(iq) > ( BUBBLEmax^3-BUBBLEmin^3 )^(1/3)
            L_Cd(iqGLOB,:) = 0;
        else
            [ZETA,wZETA] = GLL_2(nXI,BUBBLEmin, (BUBBLEmax^3-BUBBLE(iq)^3)^(1/3) ,xGLL,wGLL);
            VolZETA      = pi/6*ZETA.^3;
            C(:,iq)      = Crate_Cd(:,iq);
                      
            fZETA = zeros(nXI,1);
            for iqZETA = 1:nXI
                Lag = LagGLL_Cd(iq,:,iqZETA);
                
                DUM = 0;
                for k = 1:nXI
                    DUM = DUM + fdOLD(iqZ,k)*Lag(k);
                end%END FOR k-LOOP 1
                fZETA(iqZETA) = DUM;
            end%END FOR iqZETA-LOOP
            
            INT = 0;
            for k = 1:nXI
                INT = INT + C(k,iq)/VolZETA(k)/RHOd(iqZ)*fZETA(k)*wZETA(k);
            end%END FOR k-LOOP 2
            
            j = iq;
            jGLOB = GM(jZ,j);
            
            L_Cd(iqGLOB,jGLOB) = INT;
  
        end%END IF-LOOP
    end%END FOR iq-LOOP
end%END iqZ-LOOP


%==========================================================================
%L-MATRIX - COALESCENCE BIRTH
%==========================================================================
L_Cb = zeros(nGLOB,nGLOB);

for iqZ = 1:nZ
    jZ = iqZ;
    
    fZETA = zeros(nXI,nXI);
    for iq = 1:nXI
        Lag = zeros(nXI,nXI);
        Lag(:,:) = LagGLL_Cb(iq,:,:);
        
        for iqZETA = 1:nXI
            DUMMY = 0;
            for k = 1:nXI
                DUMMY = DUMMY + fdOLD(iqZ,k)*Lag(k,iqZETA);
            end%END k-LOOP
            fZETA(iq,iqZETA) = DUMMY;
        end%END iqZETA-LOOP
    end%END iq-LOOP
    
    
    for iq = 1:nXI
        
        iqGLOB = GM(iqZ,iq);
        
        if( BUBBLE(iq) <= (2*pi/6*BUBBLEmin^3)^(1/3) )
            L_Cb(iqGLOB,:) = 0;
        else
            
            [ZETA,wZETA] = GL_2(nXI,BUBBLEmin, (BUBBLE(iq)^3-BUBBLEmin^3)^(1/3) ,xGL,wGL);%NB GL-POINTS
            VolXI        = pi/6*BUBBLE(iq)^3;
            VolZETA      = pi/6*ZETA.^3;
            delVol       = pi/6*( BUBBLE(iq)^3 - ZETA.^3 );
            XIZETA       = ( BUBBLE(iq)^3 - ZETA.^3 ).^(2/3);
            
            for j = 1:nXI
                
                Lag    = zeros(nXI,1);
                Lag(:) = LagGLL_Cb2(iq,j,:);
                C      = Crate_Cb(:,iq);
                
                fZETA2    = zeros(nXI,1);
                fZETA2(:) = fZETA(iq,:);
                
                DUM = C./XIZETA./delVol.*Lag.*fZETA2./VolZETA./RHOd(iqZ);
                
                INT = 0;
                for k = 1:nXI
                    INT = INT + DUM(k)*wZETA(k);
                end%END k-LOOP
                
                jGLOB = GM(jZ,j);
                L_Cb(iqGLOB,jGLOB) = 0.5*BUBBLE(iq)^2*VolXI*INT;
            end%END j-LOOP
        end%END IF
    end%END iq-LOOP
end%END iqZ-LOOP


%==========================================================================
%PRECALCULATED L-MATRIX
%==========================================================================
L =  L_PRE + L_Cd - L_Cb;


%==========================================================================
%GENERATE Af=F
%==========================================================================

% LTW = L'*LAMBDAglob;
% BTW = B'*LAMBDAglob;
% 
% 
% A = LTW*L + BTW*B;
% F = LTW*g + BTW*fdBC;


A = L;
F = fdBC;
for iZ = 1:nZ
    for iXI = 1:nXI
        iGLOB = GM(iZ,iXI);
        if iZ==1
            A(iGLOB,:)=0;
            A(iGLOB,iGLOB)=1;
        end%END IF-LOOP
    end%END iXI-LOOP
end%END iZ-LOOP

%==========================================================================
%SOLVE EQUATION SYSTEM
%==========================================================================
f = A\F;

%==========================================================================
%SOLUTION IN MATRIX FORM
%==========================================================================
fdMATRIX = fdMATRIX*0;
for iXI = 1:nXI
    for iZ = 1:nZ
        iGLOB = GM(iZ,iXI);
            fdMATRIX(iZ,iXI) = f(iGLOB);
    end
end

%==========================================================================
%EVALUATE THE APPROXIMATED SOLUTION
%==========================================================================
Dum1 = A*f-F;%L*f-g; 
% Dum2 = B*f-fdBC;
% ErrorL2 = sqrt( Dum2'*Dum2 + Dum1'*Dum1 );
ErrorL2 = sqrt( Dum1'*Dum1 );
% disp(['ErrorL2 = ',num2str(ErrorL2)])

ITERerr = max(max(abs(fdOLD - fdMATRIX)))/max(fdIN);
% disp(['ITERerr = ',num2str(ITERerr)])

%==========================================================================
%UPDATE, RELAX
%==========================================================================
fdOLD = RELAX*fdMATRIX + (1-RELAX)*fdOLD;


%******************
end%END WHILE LOOP
%******************

fd = fdOLD;

disp(['ErrorL2 = ',num2str(ErrorL2)])
disp(['ITERerr = ',num2str(ITERerr)])

end%END FUNCTION