%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%MULTIFLUID-PBE MODEL FOR BUBBLE COLUMN
%ORTHOGONAL COLLOCATION METHOD
%JANNIKE SOLSVIK (jannike.solsvik@chemeng.ntnu.no)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
format compact
format short
clc
addpath('MWR_Libary','MATLAB')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global  nZ nXI nGLOB GM derCOLUMN derBUBBLE COLUMN LAMBDAglob B
global  BUBBLEmax BUBBLEmin BUBBLE
global  xGLL wGLL xGL wGL
global  Crate_Cd Crate_Cb LagGLL_Cd LagGLL_Cb LagGLL_Cb2 
global  L_B fdBC fdIN
global  LagXI LagZ
global  ITERmax ITERtol


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PHYSICAL PARAMETERS
GASCONST          = 8.3145;
Mw                = 29E-3;
P0                = 101500;                     %[Pa]
T                 = 298;                        %[K]
RHOc              = 998.01;                     %[kg/m3]
GRAVITY           = 9.8;
VolFLOWd          = 0.14;                       %Dispersed phase [m3/(m2s)]
VolFLOWc          = 1.017;                      %Continuous phase [m3/(m3s)]
ALPHA             = 0.121;                      %Void fraction at inlet (-)
SIGMA             = 0.072;
kBREAKAGE         = 6E-2; k1 = 0.336*kBREAKAGE
kCOALESCENCE      = 20E-3; kCOALESCENCE
EPSILON           = VolFLOWd*GRAVITY;
RHOd0             = Mw*P0/(GASCONST*T);         %Dispersed phase density [kg/m3], Z=Zmax
MYc               = 9.7754E-4;                  %Dynamic viscosity [kg/(ms)]
dPIPE             = 5E-2;                       %Diameter of pipe [m] 
ITERmax           = 55;
ITERtol           = 1E-10;


%Z-COORDINATE (PHYSICAL SPACE)
nZ                = 22;
COLUMNmin         = 0;
COLUMNmax         = 4;
dCOLUMN           = COLUMNmax-COLUMNmin;
Zmin              = COLUMNmin/dCOLUMN;
Zmax              = COLUMNmax/dCOLUMN;
LagZ              = eye(nZ,nZ);
[COLUMN,wCOLUMN]  = GLL_(nZ,COLUMNmin,COLUMNmax);


%XI-COORDINATE (PROPERTY SPACE)
nXI               = 45;
BUBBLEmin         = 0.5E-3;
BUBBLEmax         = 15E-3;
dBUBBLE           = BUBBLEmax-BUBBLEmin;
LagXI             = eye(nXI,nXI);
[BUBBLE,wBUBBLE]  = GLL_(nXI,BUBBLEmin,BUBBLEmax);
[xGLL,wGLL]       = GaussLobattoLegendre(nXI);
[xGL,wGL]         = GaussLegendre(nXI);


%DERIVATIVE MATRIX
derCOLUMN = LagrangeDerivativeMatrix_GLL(nZ);
derCOLUMN = derCOLUMN*2/dCOLUMN;
derBUBBLE = LagrangeDerivativeMatrix_GLL(nXI);
derBUBBLE = derBUBBLE*2/dBUBBLE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%GLOBAL INDEX NOTATION
nGLOB = nZ*nXI;

%GLOBAL MATRIX
GM = zeros(nZ,nXI);
for iZ = 1:nZ
    for iXI = 1:nXI
        iGLOB = iXI + nXI*(iZ-1);
        GM(iZ,iXI) = iGLOB;
    end
end

%GLOBAL WEIGHT
wGLOB = zeros(nGLOB,1);
for iqZ = 1:nZ
    for iqXI = 1:nXI
        iqGLOB = GM(iqZ,iqXI);
            wGLOB(iqGLOB) = wCOLUMN(iqZ)*wBUBBLE(iqXI);
    end
end
LAMBDAglob = diag(wGLOB);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INITIAL DISTRIBUTION 
% [xxx,notuse] = GLL_(nXI,1E-3,2.7);
% a = 0.7;
% b = 4;
% fdIN = a*b*xxx.^(b-1).*exp(-a*xxx.^b)*36;

%Gaussian distribution
Amp=0.195;
mean=5E-3;
sig=0.0009;
fdIN = Amp/(sig*sqrt(2*pi))*exp(-(BUBBLE-mean).^2/(2*sig^2));

plot(BUBBLE*1E3,fdIN)

fdOLD = zeros(nZ,nXI);
for k = 1:nZ
    fdOLD(k,:) = fdIN(:);
end

%BC
fdBC        = zeros(nGLOB,1);
fdBC(1:nXI) = fdIN(:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PRECALCULATIONS OF THE PBE - BOUNDARY MATRIX B AND LINEAR TERMS IN L-MATRIX  

%BOUNDARY MATRIX
B = zeros(nGLOB,nGLOB);
for iZ = 1:nZ
    for iXI = 1:nXI
        iGLOB = GM(iZ,iXI);
        if iZ==1
            B(iGLOB,iGLOB)=1;
        end%END IF-LOOP
    end%END iXI-LOOP
end%END iZ-LOOP


%PRECALCULATION - BREAKAGE FREQUENCY ---> BREAKAGE DEATH
BREAKAGEd = zeros(nXI,1);
for iq = 1:nXI
    BREAKAGEd(iq) = BREAKAGE_RATE_CT(BUBBLE(iq),BUBBLEmin,SIGMA,RHOc,EPSILON,kBREAKAGE);
end%END iq-LOOP
plot(BUBBLE*1e3,BREAKAGEd)

%PRECALCULATION - L-MATRIX ---> BREAKAGE DEATH
L_Bd = zeros(nGLOB,nGLOB);
for jZ = 1:nZ
    for jXI = 1:nXI
        jGLOB = GM(jZ,jXI);
            for iqZ = 1:nZ
                for iqXI = 1:nXI
                    iqGLOB = GM(iqZ,iqXI);
                
                        b = BREAKAGEd(iqXI);
                        L = b*LagZ(iqZ,jZ)*LagXI(iqXI,jXI);
                        
                        L_Bd(iqGLOB,jGLOB) = L_Bd(iqGLOB,jGLOB) + L;
                
                end%END iqXI
            end%END iqZ
    end%END jXI
end%END jZ


%PRECALCULATIONS - BREAKAGE FREQUENCY - REDISTRIBUTION - LAGRANGE --> BREAKAGE BIRTH
REDISTb   = zeros(nXI,nXI);
BREAKAGEb = zeros(nXI,nXI);
LagGLL_Bb = zeros(nXI,nXI,nXI);
for iq = 1:nXI
    [ZETA,wZETA]    = GLL_2(nXI,BUBBLE(iq),BUBBLEmax,xGLL,wGLL);
    
    BREAKAGEb(iq,:) = BREAKAGE_RATE_CT(ZETA,BUBBLEmin,SIGMA,RHOc,EPSILON,kBREAKAGE);
    
    %REDISTb(iq,:)   = BREAKAGE_REDIST_CT(BUBBLE(iq),ZETA);
    %REDISTb(iq,:)   = BREAKAGE_REDIST_VALENTAS(BUBBLE(iq),ZETA);
    %REDISTb(iq,:)   = BREAKAGE_REDIST_HSIA_TAVLARIDES(BUBBLE(iq),ZETA);
    REDISTb(iq,:)   = BREAKAGE_REDIST_DIEMER_OLSON(BUBBLE(iq),ZETA);

    Coord           = 2/dBUBBLE * (ZETA-BUBBLEmin) -1;
    for j = 1:nXI
        LagGLL_Bb(iq,j,:) = LagrangeGLL2(j,Coord,nXI,xGLL);
    end%END j-LOOP  
end%END iq-LOOP



%PRECALCUALTIONS - L-MATRIX ---> BREAKAGE BIRTH
L_Bb = zeros(nGLOB,nGLOB);
for iqZ = 1:nZ
    jZ = iqZ;
    
    Lag  = zeros(nXI,1);
    for iqXI = 1:nXI
    
        iqGLOB = GM(iqZ,iqXI);
    
        [ZETA,wZETA] = GLL_2(nXI,BUBBLE(iqXI),BUBBLEmax,xGLL,wGLL);
    
        VolXI = BUBBLE(iqXI)^3; %pi/6 is also in the denominator
    
        hb = REDISTb(iqXI,:)';
        b  = BREAKAGEb(iqXI,:)';
    
        DUMMY = b.*hb./ZETA.^3.*wZETA;
    
        for jXI = 1:nXI
            jGLOB = GM(jZ,jXI);
            Lag(:) = LagGLL_Bb(iqXI,jXI,:);
        
            INT = 0;
            for k = 1:nXI
                INT = INT + DUMMY(k)*Lag(k);
            end%END k
        
            L_Bb(iqGLOB,jGLOB) = VolXI*INT;
            
        end%END jXI
    end%END iqXI
end%END iqZ


%PRECALCULATION - COALESCENCE RATE --> COALESCENCE DEATH TERM
LagGLL_Cd = zeros(nXI,nXI,nXI);
Crate_Cd  = zeros(nXI,nXI);

for iq = 1:nXI
   
    [ZETA,wZETA] = GLL_2(nXI,BUBBLEmin, (BUBBLEmax^3-BUBBLE(iq)^3)^(1/3) ,xGLL,wGLL);
    Coord        = 2/dBUBBLE * (ZETA-BUBBLEmin) -1;
    
    for j = 1:nXI
        LagGLL_Cd(iq,j,:) = LagrangeGLL2(j,Coord,nXI,xGLL);
    end%END j-LOOP
    
    Crate_Cd(:,iq) = COALESCENCE_FREQUENCY(ZETA,BUBBLE(iq),RHOc,EPSILON,SIGMA,kCOALESCENCE);
end%END iq-LOOP


%PRECALCULATION - COALESCENCE RATE - LAGRANGE --> COALESCENCE BIRTH TERM
LagGLL_Cb  = zeros(nXI,nXI,nXI);
LagGLL_Cb2 = zeros(nXI,nXI,nXI);
Crate_Cb   = zeros(nXI,nXI);

for iq = 1:nXI

    [ZETA,wZETA] = GL_2(nXI,BUBBLEmin, (BUBBLE(iq)^3-BUBBLEmin^3)^(1.0/3.0) ,xGL,wGL);%NB! GL-POINTS
    BUBBLE2      = (BUBBLE(iq)^3 - ZETA.^3).^(1.0/3.0);
   
   Coord  = 2/dBUBBLE * (ZETA-BUBBLEmin) -1;
   Coord2 = 2/dBUBBLE * (BUBBLE2 - BUBBLEmin) -1;
   
   for j = 1:nXI
        LagGLL_Cb(iq,j,:)  = LagrangeGLL2(j,Coord,nXI,xGLL);
        LagGLL_Cb2(iq,j,:) = LagrangeGLL2(j,Coord2,nXI,xGLL);
   end%END j-LOOP
   
   Crate_Cb(:,iq) = COALESCENCE_FREQUENCY(ZETA(:),BUBBLE2(:),RHOc,EPSILON,SIGMA,kCOALESCENCE);
    
end%END iq-LOOP


%PROBLEM OPERATOR FOR BREAKAGE
L_B = L_Bd - L_Bb;



%**************************************************************************
%**************************************************************************
%START SOLUTION ALGORITHM OF THE MULTIFLUID-PBE MODEL
%**************************************************************************
%**************************************************************************


%==========================================================================    
%EQUATION: CONTINUOUS PHASE MOMENTUM (SIMPLIFIED EQ)
%VARIABLE: PRESURE
%INDEPENDENT COODINATE: Z
%==========================================================================
P = P0 + RHOc*GRAVITY*(COLUMNmax-COLUMN);

figure(100),subplot(3,3,1),plot(COLUMN,P),ylabel('P'),xlabel('Z')
xlim([COLUMNmin-0.2,COLUMNmax+0.2]),drawnow


%==========================================================================    
%EQUATION: IDEAL GAS LAW
%VARIABLE: DISPERSED PHASE DENSITY
%INDEPENDENT COODINATE: Z
%==========================================================================
RHOd = P*RHOd0/P0;

figure(100),subplot(3,3,2),plot(COLUMN,RHOd),ylabel('\rho_d'),xlabel('Z')
xlim([COLUMNmin-0.2,COLUMNmax+0.2]),drawnow


%==========================================================================    
%VARIABLE: DISPERSED PHASE VOLUME FRACTION
%Z=Zmin
%==========================================================================
ALPHAdIN = fdIN'/RHOd(1)*wBUBBLE;


%==========================================================================    
%VARIABLE: DISPERSED PHASE VELOCITY
%==========================================================================
vdIN = VolFLOWd/ALPHAdIN;
vcIN = VolFLOWc/(1-ALPHAdIN);
vd   = ones(nZ,nXI)*1.35;
vd(1,:) = vdIN;
figure(100),subplot(3,3,9),mesh(BUBBLE,COLUMN,vd),ylabel('Z'),xlabel('\xi'),zlabel('v_d'),drawnow


disp('----------------------------------------------')
disp('Precalulations done - enter iteration loop')
disp(' ')


%***********************************************************************************
%GLOBAL ITERATION LOOP 
ITERerrP   = 1;
ITERerrV   = 1;
ITERerrPBE = 1;
ITER1      = 0;
while  ITER1 < 5E2  &&  ITERerrP > 1E-10  &&  ITERerrPBE > 1E-10  %&&  ITERerrV > 1E-10     
%**************************************************************************************
ITER1 = ITER1+1;    
   

%==========================================================================    
%EQUATION: IDEAL GAS LAW
%VARIABLE: DISPERSED PHASE DENSITY
%INDEPENDENT COODINATE: Z
%==========================================================================
RHOd = P*RHOd0/P0;

figure(100),subplot(3,3,2),plot(COLUMN,RHOd),ylabel('\rho_d'),xlabel('Z')
xlim([COLUMNmin-0.2,COLUMNmax+0.2]),drawnow


%==========================================================================
%EQUATION: PBE
%VARIABLE: fd
%INDEPENDENT COORDINATES: Z & XI
%==========================================================================
disp('----------------------------------------------')
disp('PBE')

if ITER1==1
    fd  = PBE(vd,fdOLD,RHOd);
    fdOLD = fd;
else
    fd         = PBE(vd,fdOLD,RHOd);
    ITERerrPBE = max(max(abs(fd-fdOLD)))/max(fdIN);
    fdOLD      = 0.1*fd+0.9*fdOLD;
    fd         = fdOLD;
end

figure(100),subplot(3,3,4),mesh(BUBBLE*1E3,COLUMN,fd),ylabel('Z'),xlabel('\xi [mm]'),zlabel('f_d'),drawnow
figure(100),subplot(3,3,5),plot(BUBBLE*1E3,fd(1,:),'r--',BUBBLE*1E3,fd(nZ,:),'b'),ylabel('f_d'),xlabel('\xi [mm]'),legend('Z=0','Z=L'),drawnow

disp(' ')

%==========================================================================    
%EQUATION: DISPERSED PHASE VOLUME FRACTION
%VARIABLE: ALPHAd
%INDEPENDENT COODINATE: z
%==========================================================================
ALPHAd = zeros(nZ,1);
for iZ = 1:nZ
   ALPHAd(iZ) = (fd(iZ,:)*wBUBBLE)/RHOd(iZ);
end
ALPHAdIN = ALPHAd(1);

figure(100),subplot(3,3,3),plot(COLUMN,ALPHAd)
ylabel('\alpha_d'),xlabel('Z')
xlim([COLUMNmin-0.2,COLUMNmax+0.2]),drawnow


%==========================================================================
%VARIABLE: VELOCITY
%Z=Zmin
%==========================================================================
vcIN = VolFLOWc/(1-ALPHAdIN);
vdIN = VolFLOWd/ALPHAdIN;


%==========================================================================
%EQUATION: CONTINUOUS PHASE CONTINUITY
%VARIABLE: LIQUID PHASE VELOCITY
%INDENPENDENT COORDINATE: Z
%==========================================================================
vc = (1-ALPHAdIN)*vcIN./(1-ALPHAd);

figure(100),subplot(3,3,7),plot(COLUMN,vc),ylabel('v_c'),xlabel('Z')
xlim([COLUMNmin-0.2,COLUMNmax+0.2]),drawnow


%==========================================================================
%EQUATION: DISPERSED PHASE MOMENTUM EQUATION
%VARIABLE: DISPERSED PHASE VELOCITY
%INDEPENDENT VARIABLES: Z,XI
%==========================================================================
disp('----------------------------------------------')
disp('Dispersed velocity')
 
vdGUESS = vd;
vdNEW   = zeros(nZ,nXI);

dPdZ  = derCOLUMN*P;

Eo = zeros(nZ,nXI);
Re = zeros(nZ,nXI);
CD = zeros(nZ,nXI);

ITER2 = 0;
ITERerr2 = 1;
while ITERerr2 > 1E-10 %&& ITER2 < 50
ITER2 = ITER2+1;    
    
for iZ=1:nZ  
    for iXI=1:nXI
        dV = vdGUESS(iZ,iXI)-vc(iZ);
        Eo(iZ,iXI) = GRAVITY*(RHOc-RHOd(iZ))*BUBBLE(iXI)*BUBBLE(iXI)/SIGMA;
        Re(iZ,iXI) = RHOc*abs( dV )*BUBBLE(iXI)/MYc;
        DUM1 = 16/Re(iZ,iXI)*(1+0.15*Re(iZ,iXI)^0.687);
        DUM2 = 48/Re(iZ,iXI);
        DUM3 = 8/3*Eo(iZ,iXI)/(Eo(iZ,iXI)+4);
        %CD(iZ,iXI) = max( min(DUM1,DUM2) , DUM3 );                     %Tomiyama (1998)
        %CD(iZ,iXI) = ( 0.63 + 4.8/sqrt(Re(iZ,iXI)) )^2;                %Dalle Ville (1948)
        %CD(iZ,iXI) = 2/3*sqrt(Eo(iZ,iXI));                             %Ishii & Zuber (1979)
        %CD(iZ,iXI) = 5.645/(Eo(iZ,iXI)^-1+2.835);                      %Grevskott et al. (1996)
%         if Re(iZ,iXI) <= 1E3                                          %Schiller & Naumaan (1935)
%             CD(iZ,iXI) = 24*(1+0.15*Re(iZ,iXI)^0.687)/Re(iZ,iXI);
%         else
%             CD(iZ,iXI) = 0.44;
%         end
        CD(iZ,iXI) = 24*(1+0.10*Re(iZ,iXI)^0.75)/Re(iZ,iXI);               %Ma & Ahmadi (1990)
        
        DUM(iZ,iXI) = 3*RHOc*CD(iZ,iXI)*abs( dV )/(4*BUBBLE(iXI)*RHOd(iZ));

        vdNEW(iZ,iXI) = vc(iZ,1) - dPdZ(iZ)/(RHOd(iZ)*DUM(iZ,iXI)) - GRAVITY/DUM(iZ,iXI);         
    end
end
        vdNEW(1,:) = vdIN;
        vdNEW(nZ,:) = vdNEW(nZ-1,:); %NB!!!!!!!!!!!!!!!!
        

        ITERerr2 = max(max(abs(vdNEW-vdGUESS)))./max(max(vd));
        RELAX = 1E-1; 
        vdGUESS  = RELAX*vdNEW + (1-RELAX)*vdGUESS;
end

ITERerrV = max(max(abs(vdNEW-vd)))/max(max(vd));
RELAX = 1E-1; 
vd = RELAX*vdNEW + (1-RELAX)*vd;
figure(100),subplot(3,3,9),mesh(BUBBLE*1E3,COLUMN,vd),ylabel('Z'),xlabel('\xi [mm]'),zlabel('v_d'),drawnow 

disp(' ')


%==========================================================================
%VARIABLE: AVERAGE DISPERSED PHASE VELOCITY
%INDEPENDENT VARIABLES: Z
%==========================================================================        
vdZ = zeros(nZ,1);
for iZ = 1:nZ
    vdZ(iZ) = ( fd(iZ,:).*vd(iZ,:) )*wBUBBLE / (ALPHAd(iZ)*RHOd(iZ));
end
figure(100),subplot(3,3,8),plot(COLUMN,vdZ),ylabel('v_d(Z) [m/s]'),xlabel('Z [m]')
xlim([COLUMNmin-0.2,COLUMNmax+0.2]),drawnow


%For xi-uniform velocity:
% for iZ = 2:nZ
%     vd(iZ,:) = vdZ(iZ);
% % %     vd(iZ,:) = mean(vdZ);
% end


%==========================================================================
%EQUATION: CONSTITUTIVE EQUATION FOR THE DRAG FORCE 
%VARIABLE: fDRAG
%==========================================================================
%THE DRAG COEFFICIENT
fDRAG = zeros(nZ,nXI);
for iZ=1:nZ
    for iXI=1:nXI       
        
        dV = vd(iZ,iXI)-vc(iZ);
        Eo = GRAVITY*(RHOc-RHOd(iZ))*BUBBLE(iXI)*BUBBLE(iXI)/SIGMA;
        Re = RHOc*abs(dV)*BUBBLE(iXI)/MYc;
        
        DUM1 = 16/Re*(1+0.15*Re^0.687);
        DUM2 = 48/Re;
        DUM3 = 8/3*Eo/(Eo+4);
        %CD   = max( min(DUM1,DUM2) , DUM3 );           %Tomiyama (1998)
        %CD = ( 0.63 + 4.8/sqrt(Re) )^2;                %Dalle Ville (1948)
        %CD = 2/3*sqrt(Eo);                             %Ishii & Zuber (1979)
        %CD = 5.645/(Eo^-1+2.835);                      %Grevskott et al. (1996)
%         if Re <= 1E3                                  %Schiller & Naumaan (1935)
%             CD = 24*(1+0.15*Re^0.687)/Re;
%         else
%             CD = 0.44;
%         end 
        CD = 24*(1+0.10*Re^0.75)/Re;                    %Ma & Ahmadi (1990)

        fDRAG(iZ,iXI) = -3*RHOc*CD*fd(iZ,iXI)*abs(dV)*dV/(4*BUBBLE(iXI)*RHOd(iZ));
    end
end


%==========================================================================
%EQUATION: CONTINUOUS PHASE MOMENTUM EQUATION
%VARIABLE: LIQUID PHASE PRESSURE
%INDEPENDENT COORDINATE: Z
%==========================================================================
Re = RHOc*vc*dPIPE/MYc;
fw = (0.79*log(Re)-1.64).^-2;

dvcdZ = derCOLUMN*vc;

L = derCOLUMN;
g = -RHOc*GRAVITY*ones(nZ,1) ...
    -RHOc*vc.*dvcdZ ...
    -0.5*RHOc*fw.*vc.*vc/dPIPE ...
    -(fDRAG*wBUBBLE)./(1-ALPHAd)*1;
  
L(nZ,:)  = 0;
L(nZ,nZ) = 1;
g(nZ,1)  = 0;
g(nZ,1)  = P0;
DUM = P;
P = L\g;
RELAX = 9E-2;
P = RELAX*P + (1-RELAX)*DUM;
ITERerrP = max(abs(DUM-P))/P0;



%==========================================================================
%VARIABLE: SAUTER MEAN DIAMETER
%INDEPENDENT VARIABLES: Z
%==========================================================================
dS = zeros(nZ,1);
for iZ=1:nZ
    dS(iZ) = ( max(1E-5,fd(iZ,:))*wBUBBLE) / ((max(1E-5,fd(iZ,:))./BUBBLE')*wBUBBLE) ;
end

figure(100),subplot(3,3,6),plot(COLUMN,dS*1E3)
ylabel('dS [mm]'),xlabel('Z [m]')
xlim([COLUMNmin-0.2,COLUMNmax+0.2]),drawnow

disp('------------------------------------')
disp(['GLOBAL ITER: ',num2str(ITER1)])
disp(['ITERerrV: ',num2str(ITERerrV)])
disp(['ITERerrP: ',num2str(ITERerrP)])
disp(['ITERerrPBE: ',num2str(ITERerrPBE)])
disp(' ')

%*******************************
end%END GLOBAL ITERATION-LOOP 
%*******************************

jd = vdZ(1)*ALPHAdIN;
jc = vc(1)*(1-ALPHAdIN);

disp('------------------------------------')
disp(['Computed VolFLOWd: ',num2str(jd)]) 
disp(['Set VolFLOWd: ',num2str(VolFLOWd)]) 
disp(' ')
disp('------------------------------------')
disp(['Computed VolFLOWc: ',num2str(jc)]) 
disp(['Set VolFLOWc: ',num2str(VolFLOWc)])
disp(' ')
disp('------------------------------------')
disp(['Computed alpha: ',num2str(ALPHAdIN)]) 
disp(['Set alpha: ',num2str(ALPHA)])

%**************************************************************************
%**************************************************************************
%END SOLVING THE MULTIFLUID-PBE MODEL
%**************************************************************************
%**************************************************************************