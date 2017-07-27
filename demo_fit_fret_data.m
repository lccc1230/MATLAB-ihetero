%% Script to process FRET experimental data
%   Nonlinear fit
%   Suitable for multiple sets of data (different protein concentration)
%   K1,K2,K4,C are parameters to be estimated;
%
%  Author: Chang Lu (lu-c12@mails.tsinghua.edu.cn)
%  Paper: Quantitative analysis of ligand induced Hetero-dimerization
%
%
% NOTE: alpha represents the reverse of coorperativity factor used in the paper. 

function demo_fit_fret_data
%% data input

% data:  
% data(:,1) -- l0 conc.; 
% data(:,2~n) -- FRET signal for different conc. of P1 (which is FKBP-Cypet)
data = importdata('data/FRET_data.txt');

% total_conc: 
% total_conc(n,:) -- Conc. of total P1(FKBP) and P2(FRB) for each experiment
total_conc = [39.6 39; 55.4 39; 71.3 39];

% data_length: # of measured data points in each experiment
data_length = [12 12 12];

% construct s0 & p1lp2 for fitting (p1lp2 denots FRET signal)
% s0(m,:)=[p10, p20, l0], m is the number of total data points for all experiments
[a,b]=size(data);
s0=zeros(a*(b-1),3);
p1lp2=zeros(a*(b-1),1);
count=0;
n=b-1;
for i=1:n
    s0(count+1:count+a,1)=total_conc(i,1); 
    s0(count+1:count+a,2)=total_conc(i,2); 
    s0(count+1:count+a,3)=data(:,1);
    p1lp2(count+1:count+a)=data(:,i+1);
    count=count+a;
end

%% initial parameters
K1=1;  % D+L->DL
K2=1000; % A+L->LA
K4=10;  % DL+A->DLA
C0=10; % ratio between intensity and complex conc., proportional to FRET efficiency
coeff0=[K1 K2 K4 C0];

%% nonlinear regression
options = statset(...
    'MaxIter',1000000,...
    'TolFun',1e-26,'DerivStep',1,...
    'TolX',1e-20);
%     'Display','iter',...

mdl = NonLinearModel.fit(s0,p1lp2,@ihetero_nlinfit,coeff0,'Options',options,'CoefficientNames',{'K1','K2','K4','C'});
mdl

%% figure generation
X2=10.^linspace(-1,4,500)'; % for fitted line generation
coeff=mdl.Coefficients.Estimate;
for i=1:n
    YMatrix2(:,i)=ihetero_nlinfit(total_conc(i,1),total_conc(i,2),coeff,X2);
end 
createfigure2(data_length, [s0(:,3), p1lp2], X2, YMatrix2,[],[])
end
