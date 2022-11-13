function M3_Bootstrap_Plot
clear, clc, format short g, format compact
close all
profile on
%%%%%%%%

figure();
set(gcf,'color','w')
set(gcf,'units','centimeters','position',[5,5,20,12])
name='a':'z';

fd1='C1_CI_data.txt';
Param1=importdata(fd1); 
if isstruct(Param1)
    Param1=Param1.data;
end
Param1=rmoutliers(Param1);


fd='1CN4_LM_newbiil0.5_u1.5.txt';
Param=importdata(fd); 
if isstruct(Param)
    Param=Param.data;
end
Param=rmoutliers(Param);

Xlabel1={'\mu_{\it max}','X_{\it max}', '\alpha', '\it Y_{\it X,S}','\it Y_{\it X,N}'};
Xlabel2={'\it Y_{\it PHA,S}','\beta','\it m'};

%Xlabel2={'\mu_{\it max}','k_s', 'k_n', '\alpha', 'Y_{X/G}','Y_{X/N}', 'm', 'R^2','Adj R^2'};

for k=1:5
subplot(2,4,k);
histogram(Param(:,k),100);
title(['(' name(k) ')'])
xlabel(Xlabel1{k})
ylabel('Frequency')
end

for k=1:3
subplot(2,4,k+5);
histogram(Param1(:,k),100);
title(['(' name(k+5) ')'])
xlabel(Xlabel2{k})
ylabel('Frequency')
end

