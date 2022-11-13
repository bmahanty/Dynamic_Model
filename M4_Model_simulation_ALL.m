function M4_Model_simulation_ALL % Solving the system, check
%finding mumax,ks,kn
clear, clc, format short , format compact
close all
clear global
profile on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=importdata('Fabri_data.txt'); %% Three dataset- time-[Biomass-Carbon-Nitrogen]
Data=a.data; 
Data=Data(1:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Constants   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tsi=Data(:,1);I0=Data(1,2:end);
Para=[0.077738 0.89448 0.35444 0.018287 0.14383 1.629 0.54386 0.15836];
option1=odeset('NonNegative',4); %% Biomass C N
%% Fitness%%
[~,PredictY]=ode23(@ODEfun,tsi,I0,option1,Para);
[AllR2, avgR2, adR2, RSS,AIC, AICc] = fitness (Data(:,2:end), PredictY,Para);
fprintf('The Final average R2  : %f, Adj R2 : %f\n',avgR2, adR2);

%% Long term predict
tlong=0:0.5:71.5;
[tpred,PredictLong]=ode23(@ODEfun,tlong,I0,option1,Para);
disp('');



%% Prediction interval
Paralb=[0.0709      0.75871       0.2896     0.089082       1.2612];
Paraub=[0.0836      0.96684       0.4894      0.17262        1.763];

Xi=[];Ni=[];Ci=[];Pi=[];

for i=1:50000
Parati=Paralb+(Paraub-Paralb).*rand(1,length(Paralb));

tpi=0:0.25:71.6; %%Act_Data(:,1);
[~,PredictYi]=ode23(@ODEfun,tpi,I0,option1,Parati);
Xi=[Xi PredictYi(:,1)];
Ni=[Ni PredictYi(:,2)];
Ci=[Ci PredictYi(:,3)];
Pi=[Pi PredictYi(:,4)];

end
Xii=[min(Xi,[],2) max(Xi,[],2)];
Nii=[min(Ni,[],2) max(Ni,[],2)];
Cii=[min(Ci,[],2) max(Ci,[],2)];
Pii=[min(Pi,[],2) max(Pi,[],2)];

plot(tpi,Xii,Nii,Cii,Pii)



function dYfundt = ODEfun (~,Yfun,Param)
% S,Ef,P
X=Yfun(1);N=Yfun(2);G=Yfun(3);P=Yfun(4);
mumax=Param(1);Xmax=Param(2);alfa=Param(3);beta=Param(4);Y_XG=Param(5);Y_XN=Param(6);Y_PG=Param(7);m=Param(8);%Ki=Param(9);%tl=Param(9);%
%mumax=Param(1);ks=Param(2);kn=Param(3);alfa=Param(4);beta=Param(5);Y_XG=Param(6);Y_XN=Param(7);Y_PG=Param(8);m=Param(9);%Pmax=Param(10);%
mu=mumax*(1-(X/Xmax));%*(1-(P/Pmax))
%mu=mumax*(G/(ks+G))*(N/(kn+N));%*(1-exp(-t/tl))
%mu=mumax*(G/((ks*X)+G))*(N/((kn*X)+N));

dXdt=mu*X;
dPdt=alfa*dXdt+beta*X;
dGdt=-(1/Y_XG)*dXdt-(1/Y_PG)*dPdt-m*X;%*(Kig/(Kig+G))*(Ki/(Ki+N))
dNdt=-(1/Y_XN)*dXdt;

if N<=0
dXdt=0;
dNdt=0;
end

if G<=0 
dGdt=0;
dXdt=0;
dPdt=0;
dNdt=0;
end

% if G<=0 || N<=0
% dXdt=0;
% dGdt=0;
% dNdt=0;
% end
% 
% if G<=0 
% dPdt=0;
% end


dYfundt = [dXdt;dNdt;dGdt;dPdt]; 

function [AllR2, avgR2, adR2, RSS,AIC, AICc] = fitness (Exp, Mod,Para)
AA=sum((Exp-Mod).^2,1);
BB=sum((Exp-mean(Exp)).^2);
RSS=sum(AA);
err=AA./BB;
AllR2=1-err;
% SSE=sum(err);%([1,2,4,5,7,8])
avgR2=(sum(AllR2)/4);
N=length(Exp)-1;P=length(Para);
adR2=1-((1-avgR2)*(N*4-1))/(N*4-P-1);
fprintf('The current average R2  : %f, Adj R2 : %f\n',avgR2,adR2);

nAIC=length(Exp)*4;
K=length(Para);
AIC=nAIC*log(RSS/nAIC)+2*K;
AICc=nAIC*log(RSS/nAIC)+((2*K*nAIC)/(nAIC-K-1));
disp('')
