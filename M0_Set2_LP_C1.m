% Luedeking-piret; only sucrose-product data with 3 parameter
%
function M0_Set2_LP_C1 % Solving the system, check
clear, clc, format short g, format compact
close all
profile on
global avgR2 AdjR2 Biodata run OriPredict PredictY ParaOpt fna

fna='C1_CI_data.txt';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load TRS data');disp (' ');
a=importdata('C1_data.txt'); %% time, Biomass, Glucose, Nitrogen, PHA
Oridata=a.data;
Biodata=Oridata(1:end,:);%%% If considering all time pints
run=0;
%Oridata=Oridata(2:end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = optimset('LargeScale','on','Algorithm','interior-point','display', 'off');
lbp=[0 0 0];ubp=[1 Inf Inf];Para=[0.31  0.25016 0.0044]; %%%alfa;beta;Y_PG;m;kg;kn
%mumax= Para (1);alfa=Para(2);beta=Para (3);Y_XG= Para (4);Y_XN= Para (5);Y_PG=Para(6);m=Para(7);
[ParaOpt,Fval,~,~]=fmincon(@SumSqr,Para,[],[],[],[],lbp,ubp,[],options,Biodata);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Constants   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp ('Optimal Para '); 
disp(ParaOpt);
disp (['Final Sum of Squares: ' num2str(Fval)]); 
disp (' ');
disp([avgR2 AdjR2])
OriPredict=PredictY;

[AllR2, avgR2, adR2, RSS,AIC, AICc] = fitness (Biodata(:,2:end), OriPredict,ParaOpt);
fprintf('The Final average R2  : %f, Adj R2 : %f\n',avgR2, adR2);

%cross()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp ('Now going into CI routine for all parameters');disp(' ')
bootn=input('What many bootstrap you want to evaluate : ');disp(' ')
errM=Biodata(2:end,2:end)-OriPredict(2:end,:);
%errM=normalize(errM,'center','mean');
[ci1,~]=bootci(bootn,{@bootmodel,errM},'alpha',0.05,'type','per');%,
disp('Confidence interval for estimates');disp (' ');
disp(ci1);

Param=importdata(fna); 
if isstruct(Param)
    Param=Param.data;
end
boxplot(Param)
disp('Confidence interval on entire dataset');disp (' ');
prcmuci=prctile(Param,[2.5 97.5]);
disp(prcmuci);
disp('Hello');
hist(Param);

function SSE=SumSqr(Para,data)
global avgR2 AdjR2 PredictY CP

tspan=data(:,1);
CP=data(:,2:end); %% Three sets of Biomass/C/P
Ini=CP(1,:);
PredictY=mdl(tspan,Ini,Para);

AA=sum((CP-PredictY).^2,1);
BB=sum((CP-mean(CP)).^2);
err=AA./BB;
SSE=sum(err(:,:));
avgR2=1-(SSE/2);
N=length(tspan)-1;P=length(Para);
AdjR2=1-((1-avgR2)*(N*2-1))/(N*2-P-1);
fprintf('The current average R2  : %f, Adj R2 : %f\n',avgR2,AdjR2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y=mdl(tsi,Ini,Par)
%% p independent expt dataset
%% q variables i.e. biomass, prod
%%
global CP
p=1;
q=2;
Y=zeros(length(tsi),length(Ini));
for j=1:p %%% 3 independent expt dataset
    a=q*(j-1)+1;b=q*(j-1)+q;
Ini=CP(1,a:b);
option1=odeset('NonNegative',q); %% Biomass C N
[~,Yi]=ode23(@ODEfun,tsi,Ini,option1,Par);
Y(:,a:b)=Yi;
end

function dYfundt = ODEfun (~,Yfun,Para)
% S,Ef,P
%X=Yfun(1);N=Yfun(3);
G=Yfun(1);PHA=Yfun(2);
%mumax=0.155;
%mumax= 0.08;Kg=0.1;Kn=1.51E-05;
%Y_XG=0.6876;Y_XN=0.55973;

%mumax= 0.08;Kg=0.1;Kn=1.51E-05;Y_XG= 0.0048847;Y_XN= 0.59109;
%alfa=Para(1);;m=Para(4);
X=0.497;
%m=0.25016;
Y_PG=Para (1);beta=Para (2);m=Para (3);
%m=Para (3);
%mu=mumax*(G/(Kg+G))*(N/(Kn+N));
%dXdt=mu*X;

if G>=0
dPdt=(beta*X);
else
dPdt=0;
end

dGdt=-(1/Y_PG)*dPdt-m*X;
%dNdt=-(1/Y_XN)*dXdt;


dYfundt = [dGdt;dPdt]; 

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

function [outFrend]=bootmodel (er)
%%%%%
% [a, b]=size(er);
% for jj=1:b
%     c=er(:,jj);
%     er(:,jj)=c(randperm(a));
% end
%%%%%
global OriPredict ParaOpt avgR2 AdjR2 
global run Biodata fna

yft=OriPredict(2:end,:)+er;
er(yft<0)=0;
yf=OriPredict(2:end,:)+er;
fabriy=Biodata;
fabriy(2:end,2:end)=yf;

options = optimset('LargeScale','on','Algorithm','interior-point','display', 'off');
lbp=0.50*ParaOpt;ubp=1.5*ParaOpt;
srtPara=[0.31 0.0044 0.25016];%%[0.1 0.9 0.8 0.9];%%
[Paraboot,~,~,~]=fmincon(@SumSqr,srtPara,[],[],[],[],lbp,ubp,[],options,fabriy);
%%Scaling**
mult=[1 100 1];
Paraboot=mult.*Paraboot;
outFrend=[Paraboot avgR2 AdjR2];
run=run+1;
fprintf('This is the %d th bootstrap sample running.\n',run);

savefile(outFrend,run,fna);

function savefile(datalist,batchrun,filename)
n=length(datalist);
Rest=0;
if batchrun==1
Rest=input('would you like me to (1) start fresh table (deletes existing)  or (2) addup to existing say 1, 2: ');disp (' ');  
%Rest=2;
end


if batchrun==1 && Rest==1
%crit=Para(1);alfa=Para(2);beta=Para(3);gama=Para(4); delta=Para(5);
fid = fopen(filename, 'w');
fprintf(fid, '||- P1 P2 P3 Avg_R2 Adj_R2||-- \r\n\r\n'); 
fprintf(fid, '%10.4f  %10.4f %10.4f %10.4f %10.4f\r\n', datalist);
fclose(fid);
end

if batchrun>1 || Rest==2
fid = fopen(filename, 'a');
fprintf(fid, '%10.4f  %10.4f %10.4f %10.4f %10.4f\r\n', datalist);
fclose(fid);
end
function[]=cross()
global Biodata count CVheader
%%% cross validation
CVheader=Biodata(1,:);
datacross=Biodata(2:end,:);
count=1;

f=@(Xtrain,Xtest)cvsqrerr(Xtrain,Xtest);
PCRmsep=crossval(f,datacross,'kfold',3,'mcreps',10);%%,'mcreps',10
SSE=sum(PCRmsep,1)/10;


Rdata=Biodata(:,2:end);
SSA=sum((Rdata-mean(Rdata)).^2,1);

RMSE=(SSE/size(PCRmsep,1)).^0.5;
R2CV=1-(SSE./SSA);disp('');
fprintf('The RMSE  : %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f\n',RMSE);
fprintf('The R2CV  : %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f\n',R2CV);
function sumSqerr = cvsqrerr(Xtrain,Xtest)
global ParaOpt count CVheader
Xtraine=[CVheader;Xtrain];

lb=0.1*ParaOpt;ub=10*ParaOpt;
Param=ParaOpt;
options = optimset('LargeScale','on','Algorithm','interior-point','display', 'off','MaxFunEvals',1000);
[Paracv,~,~,~]=fmincon(@SumSqr,Param,[],[],[],[],lb,ub,[],options,Xtraine);

datacv=[CVheader;Xtest];
tspan=datacv(:,1);
cvY=mdl(tspan,CVheader(2:end),Paracv);

if size(tspan,1)==2
cvY=cvY(end,:);
else
  cvY(1,:)=[];  
end

Sqerr=(Xtest(:,2:end)-cvY).^2; %% Removing the time entry
sumSqerr=sum(Sqerr,1);
fprintf('The %d, fold with %d samples  \n',count,size(Xtest,1));
count=count+1;
