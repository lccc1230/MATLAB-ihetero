%  Generate theoretical data, and then fit the data with our model
%
%
%  Author: Chang Lu (lu-c12@mails.tsinghua.edu.cn)
%  Paper: Quantitative analysis of ligand induced Hetero-dimerization
%
%
% NOTE: alpha represents the reverse of coorperativity factor used in the paper. 

clear;
clc;
%% generate data
l0=1.0e+03.*[0.0084;0.0169;0.0337;0.0562;0.0843;0.2528;0.5617;1.1235;1.6852;3.3704;6.7407; 12];
X1=l0';
% total_conc = [39.6 39; 55.4 39; 71.3 39];
total_conc = [50 2000; 50 4000; 50 6000];
data_length = [12 12 12]; % Number of measurements in each data set
n=length(data_length); % three sets of data

K1=12;
K2=1160;
%K3=60
K4=6000;
C=10;
for i=1:n
    p10=total_conc(i,1);
    p20=total_conc(i,2);
    YMatrix1(:,i)=generate_data(K1,K2,K4,l0,p10,p20);
%     YMatrix1(:,i)=generate_data_Douglass(K1,K2,K4,l0,p10,p20);
end
%% for fitting
% construct s0 & p1lp2 for fitting (p1lp2 denots FRET signal)
% s0(m,:)=[p10, p20, l0], m is the number of total data points for all experiments
a=length(l0);
s0=zeros(a*3,3); p1lp2=zeros(a*3,1);
count=0;
for i=1:n
    s0(count+1:count+a,1)=total_conc(i,1); 
    s0(count+1:count+a,2)=total_conc(i,2); 
    s0(count+1:count+a,3)=l0;
    p1lp2(count+1:count+a)=YMatrix1(:,i).*C;
    count=count+a;
end
% createfigure2(data_length, [s0(:,3), p1lp2], [], [], [],[])
%% fitting
% initial parameters
K1_0=10;  % D+L->DL
K2_0=1000; % A+L->LA
K4_0=1000;  % DL+A->DLA
C_0=10; % ratio between intensity and complex conc., proportional to FRET efficiency
coeff0=[K1_0 K2_0 K4_0 C_0];
% coeff0=[K1 K2 K4 C];
% nonlinear regression
options = statset(...
    'MaxIter',1000,...
    'TolFun',1e-26,'DerivStep',1,...
    'TolX',1e-20);
%     'Display','iter',...

mdl = NonLinearModel.fit(s0,p1lp2,@ihetero_nlinfit,coeff0,'Options',options,'CoefficientNames',{'K1','K2','K4','C'});
mdl
%     [c,r,J] = runlsq([c0,A], [p10,p20], l0, p1lp2);
    
%% figure generation
X2=10.^linspace(-1,5,500)'; % for fitted line generation
coeff=mdl.Coefficients.Estimate;
for i=1:n
    YMatrix2(:,i)=ihetero_nlinfit(total_conc(i,1),total_conc(i,2),coeff,X2);
end     
createfigure2(data_length, [s0(:,3), p1lp2], X2, YMatrix2, [],[])