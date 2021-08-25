%% 
clc;
%clear all;
clearvars -except ban_C ban_D ban_I
%close all;

%% Data Input

time1 = datetime('3-Mar-2020','Format','dd-MMM-yy'):1:datetime('31-Dec-2021','Format','dd-MMM-yy');
time1=time1';
%load('D:\COVID_19\Model\regRatio_E_I.mat');

Bangladesh=dlmread('data_bangladeshR.txt');
%ratio=dlmread('infectDiagnoseRatio.txt');
day=Bangladesh(:,1);
I=Bangladesh(:,2);
D=Bangladesh(:,3);
C = Bangladesh (:,4);
E=Bangladesh(:,5);

degPoly =1;  % Degree of Polynomial

%% Global Parameters
N=length (day);   %Input Data Length
P=1.353e9; % India
%P=60.36e6; %Italy
%P=47e6;      %Spain    %Total Population of the Region of Interest
eN =length(day);          % Size of Data Vector will be used for parameter estimation
sN =length(time1);          %Total Simulation Time
alpha= 1/7;      % Value of Alpha
mfactor=5;

timeStep=1;
time=1:timeStep:sN;
pStep=1/timeStep;
lengTime=(length(time));

S=zeros(eN,1); %E=zeros(eN,1);
gama=zeros(eN,1);   %Gama from Real Data
lamda=zeros(eN,1);  %Lamda from Real Data
beta1=zeros (eN,1); %beta1 from real Data
eC=zeros(eN,1); eD=zeros(eN,1); dE=zeros(eN,1); dI=zeros(eN,1);
BR=zeros(eN,1); CFR=zeros(eN,1); CFRS=zeros(lengTime,1); 
DR=zeros(eN,1); DRS=zeros(lengTime,1);
CR=zeros(eN,1); CRS=zeros(lengTime,1);


%[ratioS,R_bita1]=fun_leastSquare(day(1:7),ratio(1:7), degPoly, sN+1);


for i=1:eN
 I(i)=I(i)-C(i)-D(i);
 %E(i)=abs(I(i)*regRatio.predictFcn(i));
 S(i)=P-I(i)-E(i)-C(i)-D(i);   %Suspectible Population
 CFR(i)=D(i)/I(i);
end

% for i=4:eN-3
% I(i)=(I(i)+I(i-1)+I(i-2)+I(i-3)+I(i+1)+I(i+2)+I(i+3))/7;
% E(i)=(E(i)+E(i-1)+E(i-2)+E(i-3)+E(i+1)+E(i+2)+E(i+3))/7;
% C(i)=(C(i)+C(i-1)+C(i-2)+C(i-3)+C(i+1)+C(i+2)+C(i+3))/7;
% D(i)=(D(i)+D(i-1)+D(i-2)+D(i-3)+D(i+1)+D(i+2)+D(i+3))/7;
% S(i)=(S(i)+S(i-1)+S(i-2)+S(i-3)+S(i+1)+S(i+2)+S(i+3))/7;
% end


for i=2:1:eN
    beta1(i)=abs((((E(i)-E(i-1)+alpha*E(i))*P)/((I(i)*S(i))+(mfactor*E(i)*S(i))/timeStep)));   
    gama(i)=abs(((C(i)-C(i-1))/(timeStep*I(i))));
    lamda(i)=abs((((D(i)-D(i-1))/(I(i)*timeStep))));
    BR(i)=(beta1(i)+5*beta1(i))/(gama(i)+lamda(i));
end

for i=9:eN
DR(i)=D(i)/(D(i)+C(i));
CR(i)=C(i)/(D(i)+C(i));
end
BR(1)=BR(2);

eS=2;

%%  REGRESSION APP CALLING

%[regData,RSM]=parameterEstimationOthers();
%[regData,RSM]=parameterEstimationOthers(lengTime,2,eN-1, beta1, gama, lamda);
%[regData,RSM]=parameterEstimationOthers(lengTime,day,eS,eN, beta1, gama, lamda);
eS=5;
[regData,RSM]=parameterEstimationBangladeshR(lengTime,day,eS,eN, beta1, gama, lamda);
%%
%regData=dlmread('D:\COVID_19\Model\RegressionLearner\regData.txt');

ICR=zeros(sN,1);
betaU=zeros(sN,1); betaD=zeros(sN,1); 
gamaU=zeros(sN,1); gamaD=zeros(sN,1);
lamdaU=zeros(sN,1); lamdaD=zeros(sN,1);
betaN = (regData(1:sN,1)*1.038);%1.027 %1.208
a=0.04; b=0.0035;
gamaN = (regData(1:sN,2));  %*1.485
lamdaN = (regData(1:sN,3)); %*1.7   

betaN = (regData(1:sN,1)*1.08); %For Re
a=0.020; b=0.0012;
SEMbeta = std(betaN)/sqrt(length(betaN));tsbeta = tinv([0.025  0.975],length(betaN)-1);      
SEMgama = std(gamaN)/sqrt(length(gamaN));tsgama = tinv([0.025  0.975],length(gamaN)-1);
SEMlamda = std(lamdaN)/sqrt(length(lamdaN));tslamda = tinv([0.025  0.975],length(lamdaN)-1);
for i=1:sN
betaU(i) =betaN(i) + tsbeta(2)*SEMbeta; 
betaD(i)=betaN(i)+tsbeta(1)*SEMbeta;
gamaU(i) =gamaN(i) + tsgama(2)*SEMgama; 
gamaD(i)=gamaN(i)+tsgama(1)*SEMgama;
lamdaU(i)=lamdaN(i) + tslamda(2)*SEMlamda; 
lamdaD(i)=lamdaN(i)+tslamda(1)*SEMlamda;
end
%close all;

SS=zeros (length(time),1); ES=zeros (length(time),1);
IS=zeros (length(time),1); CS=zeros (length(time),1); DS=zeros (length(time),1);

SSU=zeros (length(time),1); ESU=zeros (length(time),1);
ISU=zeros (length(time),1); CSU=zeros (length(time),1); DSU=zeros (length(time),1);

SSD=zeros (length(time),1); ESD=zeros (length(time),1); ISD=zeros (length(time),1);
CSD=zeros (length(time),1); DSD=zeros (length(time),1);
ID=zeros(length(time),1); TI=zeros(length(time),1); 

A=zeros(5,5);
derivatives=zeros(length(time),5);
IS(1)=I(1); ISD(1)=I(1); ISU(1)=I(1);
ES(1)=E(1); ESD(1)=E(1); ESU(1)=E(1);
SS(1)=P-E(1)-I(1); SSD(1)=P-E(1)-I(1); SSU(1)=P-E(1)-I(1); 

CS(1)=0; CSD(1)=0; CSU(1)=0;
DS(1)=0; DSD(1)=0; DSU(1)=0;

%RegData=zeros(length(time),3);
R0=zeros(length(time),1); R1=zeros(length(time),1); R2=zeros(length(time),1);
%lamdaSS=zeros(length(time),1);

for j=1:3
for idx=1:sN-1
  
    if j==1
    A(1,1)=-(1/P)*ISD(idx)*SSD(idx);  A(1,2)=-(1/P)*ESD(idx)*SSD(idx);
    A(2,1)=(1/P)*ISD(idx)*SSD(idx); A(2,2)=(1/P)*ESD(idx)*SSD(idx);
    A(2,3)=-ESD(idx); A(3,3)=ESD(idx); A(3,4)=-ISD(idx);
    A(3,5)=-ISD(idx);A(4,4)=ISD(idx); A(5,5)=ISD(idx);
    elseif j==2
    A(1,1)=-(1/P)*IS(idx)*SS(idx);  A(1,2)=-(1/P)*ES(idx)*SS(idx);
    A(2,1)=(1/P)*IS(idx)*SS(idx); A(2,2)=(1/P)*ES(idx)*SS(idx);
    A(2,3)=-ES(idx); A(3,3)=ES(idx); A(3,4)=-IS(idx);
    A(3,5)=-IS(idx);A(4,4)=IS(idx); A(5,5)=IS(idx);
    else
    A(1,1)=-(1/P)*ISU(idx)*SSU(idx);  A(1,2)=-(1/P)*ESU(idx)*SSU(idx);
    A(2,1)=(1/P)*ISU(idx)*SSU(idx); A(2,2)=(1/P)*ESU(idx)*SSU(idx);
    A(2,3)=-ESU(idx); A(3,3)=ESU(idx); A(3,4)=-ISU(idx);
    A(3,5)=-ISU(idx);A(4,4)=ISU(idx); A(5,5)=ISU(idx);
   end


if j==1     
    X=[betaD(idx); betaD(idx)*5; 1/7 ;a; b]; %LCL
    R0(idx+1)=(X(1)/(X(4)+X(5)))+X(2)*7;
elseif j==2
    %X=[regData(idx,1); regData(idx,1)*0; 1/7 ;regData(2); regData(3)];
    %X=[beta(idx); beta(idx)*0; 1/7 ;gama(idx); lamda(idx)];
    X=[betaN(idx); betaN(idx)*5; 1/7 ;a; b];
    %X=[betaN(idx); betaN(idx)*5; 1/7 ;gamaN(idx); lamdaN(idx)];
    R1(idx+1)=(X(1)/(X(4)+X(5)))+X(2)*7;
    ip(idx+1)=1/(X(4));
    %CFRS(idx)=DS(idx)/IS(idx);
    %CCRS(idx)=X;
    ICR(idx+1)=X(1)/(P*(X(4)+X(5)));
else
    X=[betaU(idx); betaU(idx)*5; 1/7 ; a; b]; %UCL
    R2(idx+1)=(X(1)/(X(4)+X(5)))+X(2)*7;
end
   
    
      derivatives=derCalc(A,X);
  
    if j==1
    SSD(idx+1)=derivatives(1)*timeStep+SSD(idx);
    ESD(idx+1)=derivatives(2)*timeStep+ESD(idx);
    ISD(idx+1)=derivatives(3)*timeStep+ISD(idx);
    CSD(idx+1)=derivatives(4)*timeStep+CSD(idx);
    DSD(idx+1)=derivatives(5)*timeStep+DSD(idx);
    elseif j==2
    SS(idx+1)=derivatives(1)*timeStep+SS(idx);
    ES(idx+1)=derivatives(2)*timeStep+ES(idx);
    IS(idx+1)=derivatives(3)*timeStep+IS(idx);
    CS(idx+1)=derivatives(4)*timeStep+CS(idx);
    DS(idx+1)=derivatives(5)*timeStep+DS(idx);
%     DRS(idx)=DS(idx)/(DS(idx)+CS(idx));
%     CRS(idx)=CS(idx)/(DS(idx)+CS(idx));
    TI(idx+1)=IS(idx+1)+CS(idx+1)+DS(idx+1);    
    else
    SSU(idx+1)=derivatives(1)*timeStep+SSU(idx);
    ESU(idx+1)=derivatives(2)*timeStep+ESU(idx);
    ISU(idx+1)=derivatives(3)*timeStep+ISU(idx);
    CSU(idx+1)=derivatives(4)*timeStep+CSU(idx);
    DSU(idx+1)=derivatives(5)*timeStep+DSU(idx);
    end
    if idx>3
      if IS(idx)-IS(idx-1)>0
      ID(idx+1)=IS(idx)-IS(idx-1);
      end
    end
end
end


ind_I=IS;
ind_C=CS;
ind_D=DS;

for i=1:sN
CFR(i)=DS(i)/(IS(i)+CS(i)+DS(i));
DRatio(i)=DS(i)/(DS(i)+CS(i));
CRatio(i)=CS(i)/(DS(i)+CS(i));
end



fig=figure();
set(fig,'color','white')
grid on
hold on
plot(IS(1:N),'r','DisplayName','Infected');
% plot(TI,'-.','MarkerIndices', 1:1:eN,'DisplayName','Total Infected');
% plot(ID,'--','MarkerIndices', 1:1:eN,'DisplayName','Total Infected');
plot(I,'.','MarkerIndices', 1:1:eN,'DisplayName','Real');
legend('show');
% 
fig=figure();
set(fig,'color','white')
grid on
hold on
plot(DS(1:N),'black-.','DisplayName','Death');
plot(D,'.','MarkerIndices', 1:1:eN,'DisplayName','Real');
legend('show');
% 
fig=figure();
set(fig,'color','white')
grid on
hold on
plot(CS(1:N),'m--','DisplayName','CURED');
plot(C,'.','MarkerIndices', 1:1:eN,'DisplayName','Real');
legend('show');
%%

Range='B1:B300';
% xlswrite('dInf.xlsx',IS,'Infected',Range);
% xlswrite('dInf.xlsx',CS,'Cured',Range);
% xlswrite('dInf.xlsx',DS,'Death',Range);
% xlswrite('dInf.xlsx',IS+CS+DS,'Total',Range);
xlswrite('South Asian.xlsx',IS,'Infected',Range);
xlswrite('South Asian.xlsx',CS,'Cure',Range);
xlswrite('South Asian.xlsx',DS,'Death',Range);
xlswrite('South Asian.xlsx',R0,'R0',Range);
xlswrite('South Asian.xlsx',R1,'R1',Range);
xlswrite('South Asian.xlsx',R2,'R2',Range);


%%
fig=figure();
set(fig,'color','white')
grid on
hold on
plot(time1,IS(1:pStep:lengTime),'r','DisplayName','Infected');
% plot(TI,'-.','MarkerIndices', 1:1:eN,'DisplayName','Total Infected');
% plot(ID,'--','MarkerIndices', 1:1:eN,'DisplayName','Total Infected');
plot(I,'.','MarkerIndices', 1:1:eN,'DisplayName','Real');
legend('show');
% 
fig=figure();
set(fig,'color','white')
grid on
hold on
plot(time1,DS(1:pStep:lengTime),'black-.','DisplayName','Death');
plot(D,'.','MarkerIndices', 1:1:eN,'DisplayName','Real');
legend('show');
% 
fig=figure();
set(fig,'color','white')
grid on
hold on
plot(time1,CS(1:pStep:lengTime),'m--','DisplayName','CURED');
plot(C,'.','MarkerIndices', 1:1:eN,'DisplayName','Real');
legend('show');

%% Visualization Infected
fig=figure();
set(fig,'color','white')
set(gca,'XMinorTick','on','YMinorTick','on')
grid on
hold on
plot(ISD(1:pStep:lengTime),'m-.','DisplayName','LB');
plot(IS(1:pStep:lengTime),'r','DisplayName','Simulation');
plot(ISU(1:pStep:lengTime),'black--','DisplayName','UB');
plot(I,'o','MarkerIndices', 1:10:eN,'DisplayName','Reported');
ylabel('Infected (In thousand)');
xlabel('Day');
yticks(0:5e4:max(ISU));
yticklabels(0:5e1:max(ISU)/1000);
xticklabels('auto')
xticks(0:20:lengTime);
xtickangle(90)
xlim([0 lengTime]);
ylim([0 max(ISU)]);
box('on');
legend('show','Location','northwest');

set(gca,'GridColor',[0.1 0.2 0.9])
%xticklabels('manual')

%% CURED

fig=figure();
set(fig,'color','white')
set(gca,'XMinorTick','on','YMinorTick','on')
grid on
hold on
plot(CSD(1:pStep:lengTime),'m-.','DisplayName','LB');
plot(CS(1:pStep:lengTime),'r','DisplayName','Simulation');
plot(CSU(1:pStep:lengTime),'black--','DisplayName','UB');
plot(C,'o','MarkerIndices', 1:5:eN,'DisplayName','Reported');
xlabel('Day');
ylabel('Cure ( In Thousand)');
yticks(0:5e5:max(CSU));
yticklabels(0:5e5:max(CSU));
xticklabels('auto')
xticks(0:20:lengTime);
xtickangle(90)
xlim([0 lengTime]);
ylim([0 max(CSU)]);
box('on');
legend('show','Location','northwest');
set(gca,'GridColor',[0.1 0.2 0.9])
%xticklabels('manual')
%% Death
fig=figure();
set(fig,'color','white')
grid on
hold on
plot(DSD(1:pStep:lengTime),'m-.','DisplayName','LB');
plot(DS(1:pStep:lengTime),'r','DisplayName','Simulation');
plot(DSU(1:pStep:lengTime),'black--','DisplayName','UB');
plot(D,'o','MarkerIndices', 1:5:eN,'DisplayName','Reported');
ylabel('Death ( In Thousand)');
xlabel('Day');
yticks(0:100:max(DSU));
yticklabels(0:1:max(DSU));
xticklabels('auto')
xticks(0:10:lengTime);
xtickangle(90)
xlim([0 lengTime]);
ylim([0 max(DSU)]);
xtickangle(90)
box('on');
legend('show','Location','northwest');
set(gca,'GridColor',[0.1 0.2 0.9])
%xticklabels('manual')

%% Reproduction Ratio
fig=figure();
set(fig,'color','white')
grid on
hold on
plot(R0(1:pStep:lengTime),'m-.','DisplayName','LB');
plot(R1(1:pStep:lengTime),'r','DisplayName','Simulation');
plot(R2(1:pStep:lengTime),'black--','DisplayName','UB');

ylabel('R(t)');
xlabel('Day');

xticklabels('auto')
xticks(0:5:lengTime);
yticks(0:0.5:max(R1));
xtickangle(90)
xlim([2 lengTime]);

xtickangle(90)
box('on');
legend('show','Location','northwest');
set(gca,'GridColor',[0.1 0.2 0.9])
%xticklabels('manual')

%% %% ALL
fig=figure();
set(fig,'color','white')
grid on
hold on
%plot(ES(1:pStep:lengTime),'black','DisplayName','Exposed');
plot(IS(1:pStep:lengTime),'m','DisplayName','Infected');
plot(CS(1:pStep:lengTime),'--','DisplayName','Cured');
plot(DS(1:pStep:lengTime),'black','DisplayName','Death');
ylabel('Population (Exposed, Infected, Cured, Death)');
xlabel('Day');
yticks(0:5e4:max(CS+2000));
xticklabels('auto')
xticks(0:15:lengTime);
xtickangle(90)
xlim([0 lengTime]);
ylim([0 max(CS+2000)]);
xtickangle(90)
box('on');
legend('show','Location','northwest');
set(gca,'GridColor',[0.1 0.2 0.9])
hold on
% yyaxis('right')
% ylabel('Population (Suceptible');
% plot(SS(1:pStep:lengTime),'m-.','DisplayName','Suseptible');
%xticklabels('manual')


%% 

% clc
% close all;
fig=figure();
set(fig,'color','white')
grid on
hold on
plot(time1,IS(1:pStep:lengTime),'r','DisplayName','Infected');
plot(I,'.','MarkerIndices', 1:1:eN,'DisplayName','Real');
legend('show');

fig=figure();
set(fig,'color','white')
grid on
hold on
plot(time1,DS(1:pStep:lengTime),'black-.','DisplayName','Death');
plot(D,'.','MarkerIndices', 1:1:eN,'DisplayName','Real');
legend('show');

fig=figure();
set(fig,'color','white')
grid on
hold on
plot(time1,CS(1:pStep:lengTime),'m--','DisplayName','CURED');
plot(C,'.','MarkerIndices', 1:1:eN,'DisplayName','Real');
legend('show');


%% Visualization Infected
fig=figure();
set(fig,'color','white')
grid on
hold on
plot(time1,ISD(1:pStep:lengTime),'m-.','DisplayName','LB');
plot(time1,IS(1:pStep:lengTime),'r','DisplayName','SIER');
plot(time1,ISU(1:pStep:lengTime),'black--','DisplayName','UB');
plot(I,'o','MarkerIndices', 1:3:eN,'DisplayName','Real');
ylabel('Infected(In thousand)');
yticks(0:5e3:max(ISU));
yticklabels(0:5:max(ISU/1000));
xticklabels('auto')
xticks([datetime(min(time1)):3:datetime(max(time1))]);
xtickangle(90)
xlim([min(time1) max(time1)]);
ylim([0 max(ISU)]);
box('on');
legend('show','Location','northwest');

set(gca,'GridColor',[0.1 0.2 0.9])
%xticklabels('manual')


%% CURED

fig=figure();
set(fig,'color','white')
grid on
hold on
plot(time1,CSD(1:pStep:lengTime),'m-.','DisplayName','LB');
plot(time1,CS(1:pStep:lengTime),'r','DisplayName','SIER');
plot(time1,CSU(1:pStep:lengTime),'black--','DisplayName','UB');
plot(C,'o','MarkerIndices', 1:3:eN,'DisplayName','Real');
ylabel('Cure');
yticks(0:250:max(CSU));
xticklabels('auto')
xticks([datetime(min(time1)):3:datetime(max(time1))]);
xlim([min(time1) max(time1)]);
ylim([0 max(CSU)]);
xtickangle(90)
xticklabels('Manual')
box('on');
legend('show','Location','northwest');
set(gca,'GridColor',[0.1 0.2 0.9])
%xticklabels('manual')
%% Death
fig=figure();
set(fig,'color','white')
grid on
hold on
plot(time1,DSD(1:pStep:lengTime),'m-.','DisplayName','LB');
plot(time1,DS(1:pStep:lengTime),'r','DisplayName','SIER');
plot(time1,DSU(1:pStep:lengTime),'black--','DisplayName','UB');
plot(D,'o','MarkerIndices', 1:3:eN,'DisplayName','Real');
ylabel('Death');
yticks(0:25:max(DSU));
xticklabels('auto')
xticks([datetime(min(time1)):3:datetime(max(time1))]);
xlim([min(time1) max(time1)]);
ylim([0 max(DSU)]);
xtickangle(90)
xticklabels('Manual')
box('on');
legend('show','Location','northwest');
set(gca,'GridColor',[0.1 0.2 0.9])
%xticklabels('manual')

%% 
fig=figure();
set(fig,'color','white')
grid on
hold on
%plot(time1,SSD(1:pStep:lengTime),'m-.','DisplayName','LB');
hold on
yyaxis right
yticks([1.7e7:1.8e8])
yticklabels([170:1:180])
plot(time1,SS(1:pStep:lengTime),'black','DisplayName','Suspeectible');
hold on
% yyaxis left
% hold on
% plot(time1,ES(1:pStep:lengTime),'-.','DisplayName','Exposed');
% plot(time1,IS(1:pStep:lengTime),'b','DisplayName','Infected');
% plot(time1,CS(1:pStep:lengTime),'m','DisplayName','Cured');
% plot(time1,DS(1:pStep:lengTime),'--','DisplayName','Death');
%plot(time1,SSU(1:pStep:lengTime),'black--','DisplayName','UB');
%plot(S,'o','MarkerIndices', 1:5:eN,'DisplayName','Real');
%xlabel('Day');
xticklabels('manual')
xtickangle(45)
ylabel('Susceptible');
ylabel('Exposed, Infected, Cured and Death)');
yticks(0:3e4:15e4);
box('on');
legend('show');
%% %% Visualization All 
fig=figure();
set(fig,'color','white')
grid on
hold on
plot(time1,ES(1:pStep:lengTime),'m-.','DisplayName','Affected');
plot(time1,IS(1:pStep:lengTime),'b-.','DisplayName','Infected');
plot(time1,CS(1:pStep:lengTime),'r','DisplayName','Cured');
plot(time1,DS(1:pStep:lengTime),'black--','DisplayName','Death');
plot(I,'r.','MarkerIndices', 1:3:eN,'DisplayName','Reported(Infected)');
% plot(C,'g.','MarkerIndices', 1:3:eN,'DisplayName','Reported(Cured)');
% plot(D,'b.','MarkerIndices', 1:3:eN,'DisplayName','Reported(Death)');
ylabel('Infected(In thousand)');
yticks(0:5e3:max(ISU));
yticklabels(0:5:max(ISU/1000));
xticklabels('auto')
xticks([datetime(min(time1)):15:datetime(max(time1))]);
xtickangle(90)
xlim([min(time1) max(time1)]);
ylim([0 max(ISU)]);
box('on');
legend('show','Location','northwest');

set(gca,'GridColor',[0.1 0.2 0.9])
%xticklabels('manual')

%% Exposed
fig=figure();
set(fig,'color','white')
grid on
hold on
plot(time1,ESD(1:pStep:lengTime),'m-.','DisplayName','LB');
plot(time1,ES(1:pStep:lengTime),'r','DisplayName','SIER');
plot(time1,ESU(1:pStep:lengTime),'black--','DisplayName','UB');
%plot(E,'o','MarkerIndices', 1:5:eN,'DisplayName','Real');
xlabel('Day');
ylabel('Exposed (In thousand)');
yticks(0:2e4:8e4);
box('on');
legend('show');
hold on
xticks(b,[datetime('Mar-08','Format','MMM-dd'):3:datetime('May-31','Format','MMM-dd')])

xtickangle(45)
yticks([0:5000:75000])
yticklabels([0:5:75])
grid on
set(gca,'GridColor',[0.1 0.2 0.9])
xticklabels('manual')

%% Fertility and CUred Rate
fig=figure();
set(fig,'color','white')
grid on
hold on
plot(time1,CFRS,'m-.','DisplayName','CFR(Simulation');
plot(CFR,'-.','MarkerIndices', 1:5:eN,'DisplayName','CFR(Real)');

xlabel('Day');
ylabel('Case Fertility Rate');
legend('show');

%% Death and Cure Rate
fig=figure();
set(fig,'color','white')
grid on
hold on
plot(time1,CRS,'m-.','DisplayName','CR(Simulation');
plot(CR,'-.','MarkerIndices', 1:5:eN,'DisplayName','CR(Real)');
plot(time1,DRS,'black','DisplayName','CR(Simulation');
plot(DR,'+','MarkerIndices', 1:5:eN,'DisplayName','CR(Real)');
xlabel('Day');
ylabel('Case Fertility Rate');
legend('show');
xim()
%% Reproduction Rate
fig=figure();
set(fig,'color','white')
grid on
hold on
%plot(time1,R0(1:pStep:lengTime),'b','DisplayName','LB');
plot(time1,R1(1:pStep:lengTime),'r','DisplayName','SIER');
%plot(time1,R2(1:lengTime),'black--','DisplayName','LB');
% plot(BR,'o','DisplayName','Real');

xlabel('Time (Day)');
ylabel('Basic Reproduction Number (R_0)');
legend('show');

xlim(datetime('Mar-08','Format','MMM-dd'):3: datetime('May-31','Format','MMM-dd'));

%% 
xticks([datetime('Mar-08','Format','MMM-dd'):datetime('May-31','Format','MMM-dd')])
yticks([0:1000:8000])
yticklabels([0:1:8])
grid on
set(gca,'GridColor',[0.1 0.2 0.9])

%% CONDITIONS

q=zeros(sN,1);
R0=zeros(sN,1);
L0=zeros(sN,1);

for i=1:sN
  
  q(i)=(X(1))/(P*(X(4)+X(4)));
  R0(i)=R1(i)*P;  
  L0(i)=IS(i)/(IS(i)+5*ES(i));
  
  if  R0(i)<=L0(i)
     disp("Inequality Denies, EPIDEMIC IS UNDER CONTROL ")
     disp(i)
  end    
  
end
plot(time1,R1,'DisplayName','Reproduction Number');
hold on
plot(L0,'DisplayName','Contact Rate');
ylim([-0.5 12])

xticklabels('auto')
xticks([datetime(min(time1)):10:datetime(max(time1))]);
xtickangle(90)
xlim([min(time1)+1 max(time1)-1]);

box('on');
legend('show','Location','northwest');

set(gca,'GridColor',[0.1 0.2 0.9])

%% % for i=2:eN
%     beta1(i)=(((E(i)-E(i-1)+alpha*E(i))*P)/((I(i)*S(i))+(mfactor*E(i)*S(i))/timeStep));   
%     %beta1(i)=((S(i)-S(i-1))/timeStep)/(-((I(i)*S(i))/P)-(mfactor*E(i)*S(i)/P));  
%     gama(i)=abs((C(i)-C(i-1))/(timeStep*I(i)));
%     lamda(i)=((D(i)-D(i-1))/(I(i)*timeStep));
%     
% end
% beta2=mfactor*beta1; % Calculation of Beta2

% %% %%  LEAST SQUARE FUNCTION CALLING
% 
% % Infacter rate bita1 estimation
% [beta1S,R_bita1]=fun_leastSquare(day(1:eN),beta1(1:eN), degPoly, sN);
% %Death Rate lamda estimation
% [gamaS,R_lamda]=fun_leastSquare(day(1:eN),gama(1:eN), degPoly, sN);
% %Cure Rate gama estimation
% [lamdaS,R_gama]=fun_leastSquare(day(1:eN),lamda(1:eN),degPoly, sN);
% 
% 
% beta2S=mfactor*beta1S;   % beta2 Calculation

%% 
%S(t)=((t/Pop)*bita1*t*(I(t)-I(1))*S(1)+(t/Pop)*bita2*(E(t)-E(1))*S(1))/(1+(t/Pop)*bita1*(I(t)-I(1))+(t/Pop)*bita2*(E(t)-E(1)));
%E(t)= (E(1)+(1/Pop)*(bita1*t*(I(t)-I(1))*(S(t)-S(1)))-(1/Pop)*(bita2*t*E(1)*(S(t)-S(1)))+alpha*t*E(1))/(1-(1/Pop)*t*(S(t)-S(1))+alpha);
%I(t)=(I(1)+alpha*t*(E(t)-E(1))+gama(t)*t*I(1)+lamda(i)*t*I(1))/(1+gama(t)+lamda(t));
%I(t)=(I(1)+alpha*t*(E(t)-E(1))+gama(t)*t*I(1)+lamda(i)*t*I(1))/(1+gama(t)+lamda(t));
%C(t)=C(1)+gama*t*(I(t)-I(1));
%D(t)=D(1)+lamda*t*(I(t)-I(1));

%dS(i)=-(1/Pop)*bita1(t)*(I(t)-I(t0))*(S(t)-S(t0))-(1/Pop)*bita2(t)*(E(t)-E(t0))*(S(t)-S(t0));
%dE(i)=(1/Pop)*bita1(t)*(I(t)-I(t0))*(S(t)-S(t0))+(1/Pop)*bita2(t)*(E(t)-E(t0))*(S(t)-S(t0))-alpha*(E(t)-E(t0));
%dI(i)=alpha*(E(t)-E(t0))-gama(t)*(I(t)-I(t0))-lamda(t)*(I(t)-I(t0));
%dC(i)=gama(t)*(I(t)-I(t0));
%dD(i)=lamda(t)*(I(t)-I(t0));
 
%% Not Necessary

% for j=1:last-1
%     st=(j-1)*timeWindow+1;
%     en= j*timeWindow;
%   for i= st : en 
%     beta1(i)=((E(i+1)-E(i)+alpha*E(i+1))*P)/((I(i+1)*S(i+1))+(5*E(i+1)*S(i+1)));   
%     gama(i)=(C(i+1)-C(i))/I(i+1);
%     lamda(i)=(D(i+1)-D(i))/I(i+1);
%   end
%   abeta1(j)=mean(beta1(st:end));
%   agama(j)=mean(gama(st:end));
%   alamda(j)=mean(lamda(st:end));
% end
