function M4_Model_Sensitivity_Plot % Solving the system, check
%%% here biomass profile not required and adopted as parameter
%%http://www.scielo.br/scielo.php?pid=S0104-66322014000200007&script=sci_arttext#taba1
clear, clc, format short g, format compact
close all
profile on
global ParaOpt Oridata% %Bio0 Sub0 Prod0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load experimental data file name COLUMN WISE: time, Biomass, Glucose, TTC, Formazan');disp (' ');
rawdata=importdata('CN4_data_DWCA.txt'); 
Oridata=rawdata.data;
Oridata([4 9 12],:)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ParaOpt=[0.077738 0.89448 0.35444 0.018287 0.14383 1.629 0.54386 0.15836];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Constants   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pro=input('number of experimentally measued profiles ? ');
Yn=input('number of profiles in modelled Y vector? ');
sental=[];
parsiz=size(ParaOpt,2);
for i=1:parsiz
    
    sen=relsen(i,pro,Yn,ParaOpt);
    sental=[sental;sen]; %#ok<*AGROW>
end
disp (' ');
%% data storage
fid = fopen('Sensitivity_data.txt', 'w');
fprintf(fid, '%10.4f %10.4f %10.4f %10.4f \r\n', sental'); %% on three response
fclose(fid);
%% Sensitivity plot
a=importdata('Sensitivity_data.txt');
sensiplot(a,parsiz,pro)


function sensiplot(data,par,profi)
disp('we are reshaping parameter sensitivities for each of three output')
resen=reshape(data,[],par*profi); %#ok<*NASGU>
resen=resen(1:8,:);%%% only first ten time values/data points
figure();
set(gcf,'color','w')
s='abcd';
a=1;b=par;
for i=1:profi
subplot(2,2,i)
c=i*b;
a=((i-1)*b)+1;
bar3(resen(:,a:c));

set(gca, 'XTickLabel', {'\mu_{\it max}','X_{max}', '\alpha', '\beta','Y_{X,S}','Y_{X,N}','Y_{PHA,S}','\it m'},'FontSize',12);  % Modify the y axis tick labels
set(gca, 'YTickLabel', {'15' '25' '41' '43' '45' '47' '63' '65'},'FontSize',12);  % Modify the y axis tick labels 15    25    41    43    45    47    63    65
xlabel('Parameters','rotation',25, 'FontSize',12)
ylabel('Time (h)','rotation',330, 'FontSize',12)
zlabel('Rel sensitivity')
title(['(',s(i),')'])
end

function ZZavg=relsen(pi,nYexe,nYmod,FixPara)

%%% pi - Parameter number for analysis
%%% nYexe - num of profiles experimentally measured
%%% nYmod - number of profiles in modelled Y vector
global Oridata 
time=Oridata(:,1);Initial=Oridata(1,2:5);
ndp=size(time,1);
sensi=[FixPara(pi) FixPara(pi)*1.01 FixPara(pi)*0.99]; %% with first parameter sensitivity
%There are 6 different measured profile, so "z" will have 6x3=18 columns
z=zeros(ndp,nYexe*3);mm=1;pp=1;ZZ=zeros(ndp,nYexe*2);ZZavg=zeros(ndp,nYexe);

for i=1:3 %% original and two modified para
FixPara(pi)=sensi(i);
%disp ('********************************************************************** '); disp (' ');
option1=odeset('NonNegative',1:nYmod);
[~,y]=ode23(@ODEfun,time,Initial,option1,FixPara);
%%% Biomass Sucrose Glucose Fructose Viability Product

z(:,[i,i+3,i+6,i+9])=[y(:,1) y(:,2) y(:,3) y(:,4)]; 
end

for k=1:nYexe  %%% beacuse fot UREA, DIN Model PH Ca
ZZ(:,mm)=((z(:,pp)-z(:,pp+1))./z(:,pp)).*(sensi(1)/(sensi(1)-sensi(2)));
ZZ(:,mm+1)=((z(:,pp)-z(:,pp+2))./z(:,pp)).*(sensi(1)/(sensi(1)-sensi(3)));
ZZavg(:,k)=(ZZ(:,mm)+ZZ(:,mm+1))./2;
mm=mm+2;
pp=pp+3;
end
ZZavg(1,:)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dYfundt = ODEfun (~,Yfun,Param)
%% critF/X, alfa beta
%% Time, Biomass, Glucose, formazan
mumax=Param(1);Xmax=Param(2);alfa=Param(3);beta=Param(4);Y_XG=Param(5);Y_XN=Param(6);Y_PG=Param(7);m=Param(8);
X=Yfun(1);N=Yfun(2);G=Yfun(3);P=Yfun(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu=mumax*(1-(X/Xmax));

dXdt=mu*X;
dPdt=alfa*dXdt+beta*X;
dGdt=-(1/Y_XG)*dXdt-(1/Y_PG)*dPdt-m*X;
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

dYfundt = [dXdt;dNdt;dGdt;dPdt]; 

