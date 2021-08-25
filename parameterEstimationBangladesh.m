function [regData,RSM]=parameterEstimationBangladesh(lengTime,day,eS,eN, beta1, gama, lamda)

gamaR(:,1)=day(eS:eN); gamaR(:,2)=gama(eS:eN);
%betaR(:,1)=1:(eN-eS)+1; betaR(:,2)=beta1(eS:eN);
betaR(:,1)=day(eS:eN); betaR(:,2)=beta1(eS:eN);
lamdaR(:,1)=day(eS:eN); lamdaR(:,2)=lamda(eS:eN);

% [betaReg, RMSE2] = trainRegressionModelFineG(betaR);
% [gamaReg, RMSE1] = trainRegressionModelFineG(gamaR);
% [lamdaReg,RMSE3] =trainRegressionModelFineG(lamdaR);

% [betaReg, RMSE2] = trainRegressionModelCoarse(betaR);
 [gamaReg, RMSE1] = trainRegressionModelCoarse(gamaR);
% [lamdaReg,RMSE3] =trainRegressionModelCoarse(lamdaR);
%[gamaReg, RMSE1] = trainRegressionModelBeta(gamaR);
%[betaReg, RMSE2] = trainRegressionModelLinearSMV(betaR);
%[betaReg, RMSE2] = trainRegressionModelLinearSMV(betaR);
%[betaReg, RMSE2] = trainRegressionModelCoarse(betaR);
%[betaReg, RMSE2] = trainRegressionModelCubicSMV(betaR);
[betaReg, RMSE2] =trainRegressionModelPak(betaR);
%[betaReg, RMSE2] =trainRegressionModelPExpo(betaR);
%[betaReg, RMSE2] =trainRegressionModelFG(betaR);
%[betaReg, RMSE2] =trainRegressionModelLinear(betaR);
[lamdaReg,RMSE3] =trainRegressionModelCoarse(lamdaR);
clc
regData=zeros(lengTime,3);
for i=1:lengTime
  if betaReg.predictFcn(i)>0
     regData(i,1)=betaReg.predictFcn(i);
 else
     regData(i,1)=0;
  end 

 if gamaReg.predictFcn(i)>0
     regData(i,2)=gamaReg.predictFcn(i);
 else
     regData(i,2)=0;
 end
 
  if lamdaReg.predictFcn(i)>0
      regData(i,3)=lamdaReg.predictFcn(i);
  else
     regData(i,3)=0;
  end
 RSM=[RMSE1,RMSE2,RMSE3];
end