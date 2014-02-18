%Drag coefficient
clear all
close all
clc

POS=[500 500 320 240];


f1=figure(1);
set(f1,'DefaultAxesLineStyleOrder','-|-.|:|-o|-d','DefaultAxesColorOrder',[0 0 0],'DefaultLineLineWidth',1.4,'Position',POS)
hold all,box on
    load('Case1','dS','COLUMN'),plot(COLUMN,dS*1E3,'k')
    %load('Case2','dS','COLUMN'),plot(COLUMN,dS*1E3,'k')
    %load('Case3','dS','COLUMN'),plot(COLUMN,dS*1E3,'k')
    %load('Case4','dS','COLUMN'),plot(COLUMN,dS*1E3,'k')
    %load('Case5','dS','COLUMN'),plot(COLUMN,dS*1E3,'k')
    load('Case7','dS','COLUMN'),plot(COLUMN,dS*1E3,'k')
xlabel('z [m]')
ylabel('\xi_{sdm} [mm]')
legend('Tomiyama','Ville','Ishii')



f2=figure(2);
set(f2,'DefaultAxesLineStyleOrder','-|-.|:|-o|-d','DefaultAxesColorOrder',[0 0 0],'DefaultLineLineWidth',1.4,'Position',POS)
hold all,box on
    load('Case1','vdZ','COLUMN'),plot(COLUMN,vdZ,'k')
    %load('Case2','vdZ','COLUMN'),plot(COLUMN,vdZ,'k')
    %load('Case3','vdZ','COLUMN'),plot(COLUMN,vdZ,'k')
    %load('Case4','vdZ','COLUMN'),plot(COLUMN,vdZ,'k')
    %load('Case5','vdZ','COLUMN'),plot(COLUMN,vdZ,'k')
    load('Case7','vdZ','COLUMN'),plot(COLUMN,vdZ,'k')
xlabel('z [m]')
ylabel('\xi_{sdm} [mm]')
legend('Tomiyama','Ville','Ishii')