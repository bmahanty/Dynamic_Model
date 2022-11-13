function M7_Model_corre % Solving the system, check
%finding mumax,ks,kn
clear, clc, format short g, format compact
close all
profile on
global fna OriPredict ParaOpt run Biodata

fna='1CN4_corre.txt';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=importdata('CN4_data_DWCA.txt'); %% Three dataset- time-[Biomass-Carbon-Nitrogen]
Oridata=a.data; 
Biodata=Oridata(1:end,:);
Biodata([4 9 12],:)=[];%(1:18,:)
run=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = optimset('LargeScale','on','Algorithm','interior-point','display', 'off');

Para=[0.0907 3 1 0.25016 0.5 0.6 0.31 0.0044];
%Para=[0.065484        4.375     0.056288       1.7894       0.1322      0.11405];
lbp=[0 0 0 0 0 0 0 0];ubp=[Inf Inf Inf Inf Inf Inf 1 Inf];   
%lbp=[0 0 0 0 0 0 0];ubp=[Inf Inf Inf Inf Inf Inf Inf]; %%%mumax;Kg Kn; Y_GX=Para(6);Y_NX=Para(8);%[initial guess values0.0025847 1.8335 0.9 0.676 0.463 0.001]1.8335 0.9 0.676 0.463 0.001 [0.005847 15.335 2 0.276 0.463 0.001]
[ParaOpt,~,~,~]=fmincon(@SumSqr,Para,[],[],[],[],lbp,ubp,[],options,Biodata);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Constants   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp (ParaOpt);
rng('default')
[~,OriPredict]=ode45(@ODEfun,Biodata(:,1),Biodata(1,2:end),odeset('NonNegative',4),ParaOpt);
% nerr=[normrnd(0,0.07,[8,1]) normrnd(0,0.01,[8,1]) normrnd(0,0.4,[8,1]) normrnd(0,0.02,[8,1])];
% Mod=OriPredict(2:end,:)+nerr;
% OriPredict(2:end,:)=Mod;


[AllR2, avgR2, adR2, RSS,AIC, AICc] = fitness (Biodata(:,2:end), OriPredict,ParaOpt);
fprintf('The Final average R2  : %f, Adj R2 : %f\n',avgR2, adR2);

%cross()
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp ('Now going into CI routine for all parameters');disp(' ')
bootn=input('What many bootstrap you want to evaluate : ');disp(' ')

errM=Biodata(2:end,2:end)-OriPredict(2:end,:);
% M=mean(errM);
% Stdev=std(errM);
% nerr=[normrnd(0,Stdev(1),[8,1]) normrnd(0,Stdev(2),[8,1]) normrnd(0,Stdev(3),[8,1]) normrnd(0,Stdev(4),[8,1])];
%errM=errM/10;
%errM=normalize(errM,'center','mean');
%Error=((errM-M)./Stdev);
[ci1,~]=bootci(bootn,{@bootmodel,errM},'alpha',0.05,'type','per');%,

Param=importdata(fna); 
if isstruct(Param)
    Param=Param.data;
 end
[R,Pval]= corrcoef(tb1(:,1:8));

% disp('Confidence interval for estimates');disp (' ');
% disp(ci1);
% disp('');
% 
%  figure();
%  boxplot(Param)
%  Param=rmoutliers(Param);
% disp('Confidence interval on entire dataset');disp (' ');
% prcmuci=prctile(Param,[2.5 97.5]);
% disp(prcmuci);
disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function SSE=SumSqr(Para,data)


tspan=data(:,1);
XCN=data(:,2:end); %% Three sets of Biomass/N/C
Ini=XCN(1,:);
option1=odeset('NonNegative',4);
[~,PredictY]=ode45(@ODEfun,tspan,Ini,option1,Para);

AA=sum((XCN-PredictY).^2,1);
BB=sum((XCN-mean(XCN)).^2);
err=AA./BB;
SSE=sum(err);%([1,2,4,5,7,8])
avgR2=1-(SSE/4);
N=length(tspan)-1;P=length(Para);
AdjR2=1-((1-avgR2)*(N*4-1))/(N*4-P-1);
fprintf('The current average R2  : %f, Adj R2 : %f\n',avgR2,AdjR2);



function [outFrend]=bootmodel (er)
%%%%%
% [a, b]=size(er);
% for jj=1:b
%     c=er(:,jj);
%     er(:,jj)=c(randperm(a));
% end
%%%%%
global OriPredict ParaOpt run Biodata fna


yft=OriPredict(2:end,:)+er;
er(yft<0)=0;
yf=OriPredict(2:end,:)+er;
fabriy=Biodata;
fabriy(2:end,2:end)=yf;

options = optimset('LargeScale','on','Algorithm','interior-point','display', 'off');
%lbp=[0 0 0 0 0 0 0 0];ubp=[Inf Inf Inf Inf Inf Inf 1 Inf]; 
lbp=0.5*ParaOpt;ubp=1.5*ParaOpt;%ubp(4)=1;

%srtPara=ParaOpt;
srtPara=[0.0907 3 1 0.25016 0.5 0.6 0.31 0.0044];
[Paraboot,~,~,~]=fmincon(@SumSqr,srtPara,[],[],[],[],lbp,ubp,[],options,fabriy);

mult=[1 1 1 1 1 1 1 1];
[~,FebPred]=ode45(@ODEfun,fabriy(:,1),fabriy(1,2:end),odeset('NonNegative',4),Paraboot);

[~,BavgR2, BadjR2,~,~,~] = fitness (fabriy(:,2:end), FebPred,Paraboot);
Parabooti=mult.*Paraboot;
outFrend=[Parabooti BavgR2 BadjR2];
run=run+1;
fprintf('This is the %d th bootstrap sample running.\n',run);

savefile(outFrend,run,fna);
disp('');

function Y=mdl(tsi,Ini,Par)
global p q

Y=zeros(length(tsi),length(Ini));
for j=1:p %%% 3 independent expt dataset
    a=q*(j-1)+1;b=q*(j-1)+q;
    I0=Ini(a:b);
    option1=odeset('NonNegative',q); %% Biomass C N
    [~,Yi]=ode23(@ODEfun,tsi,I0,option1,Par);
    Y(:,a:b)=Yi(:,:);
%     if length(tsi)==2
%         Y(:,a:b)=Yi([1 end],:);
%     else
%         Y(:,a:b)=Yi(:,:);
%     end
end

function dYfundt = ODEfun (~,Yfun,Param)
% S,Ef,P
X=Yfun(1);N=Yfun(2);G=Yfun(3);P=Yfun(4);
%mumax=Param(1);ks=Param(2);kn=Param(3);alfa=Param(4);Y_XG=Param(5);Y_XN=Param(6);%m=Param(7);%Pmax=Param(8);
%mumax=Param(1);Xmax=Param(2);alfa=Param(3);Y_XG=Param(4);Y_XN=Param(5);%m=Param(6);%Ki=Param(7);
mumax=Param(1);Xmax=Param(2);alfa=Param(3);beta=Param(4);Y_XG=Param(5);Y_XN=Param(6);Y_PG=Param(7);m=Param(8);
%Y_PG=0.54386;beta=0.018287;m=0.15836;% 0.97044Y_PG=Param(6);
%Y_PG=0.55049;
%ks=Para(1);kn=Para(2);Y_XG=Para(3);Y_XN=Para(4);m=Para(5);
%mu=mumax*(1-exp(-G/kt));
%mu=mumax*(G/(ks+G+(G^2/Ki)))*(N/(kn+N));
%mu=mumax*(1-exp(-ks*G))*(1-exp(-kn*N));
mu=mumax*(1-(X/Xmax));
%mu=mumax*(G/(ks+G))*(N/(kn+N));%*(1-(P/Pmax));
%mu=mumax*(G/((ks*X)+G))*(N/((kn*X)+N));

dXdt=mu*X;
dPdt=alfa*dXdt+beta*X;
%dPdt=qmax*(G/(ks+G))*X;
dGdt=-(1/Y_XG)*dXdt-(1/Y_PG)*dPdt-m*X;
dNdt=-(1/Y_XN)*dXdt;

% if G<=0 || N<=0
% dXdt=0;
% dGdt=0;
% dNdt=0;
% end
% 
% if G<=0 
% dPdt=0;
% end

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


function savefile(datalist,batchrun,filename)
n=length(datalist);
C = {'%10.4f ';'\r\n'};
spec=[C{[ones(1,n) 2]}];
Rest=0;
if batchrun==1
Rest=input('would you like me to (1) start fresh table (deletes existing)  or (2) addup to existing say 1, 2: ');disp (' ');  
%Rest=2;
end


if batchrun==1 && Rest==1
%crit=Para(1);alfa=Para(2);beta=Para(3);gama=Para(4); delta=Para(5);
fid = fopen(filename, 'w');
fprintf(fid, '||- P1 P2 P3 P4 Avg_R2 Adj_R2||-- \r\n\r\n'); 
fprintf(fid, spec, datalist);
fclose(fid);
end

if batchrun>1 || Rest==2
fid = fopen(filename, 'a');
fprintf(fid, spec, datalist);
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
n=length(RMSE);
C = {'%2.4f ';'\n'};
spec=[C{[ones(1,n) 2]}];

fprintf(spec,RMSE);
fprintf(spec,R2CV);

function sumSqerr = cvsqrerr(Xtrain,Xtest)
global ParaOpt count CVheader
Xtraine=[CVheader;Xtrain];

lb=0.5*ParaOpt;ub=1.5*ParaOpt;
%Param=ParaOpt;
Param=[0.0907 3 1 0.5 0.6];
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

