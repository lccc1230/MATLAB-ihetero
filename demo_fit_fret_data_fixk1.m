%% Script to process FRET experimental data
%   Nonlinear fit
%   Suitable for multiple sets of data (different protein concentration)
%   K1 is fixed to be 0.3
%   K2,K4,C are parameters to be estimated;
%
%  Author: Chang Lu (lu-c12@mails.tsinghua.edu.cn)
%  Paper: Quantitative analysis of ligand induced Hetero-dimerization
%
%
% NOTE: alpha represents the reverse of coorperativity factor used in the paper. 

function demo_fit_fret_data_fixK1 
%% data input
%   data: Experiment data
%   data(:,1) is ligand concentration
%   data(:,2~n) is fluorescence intensity for each experiment
data = importdata('data/FRET_data.txt');
total_conc = [39.6 39; 55.4 39; 71.3 39];
data_length = [12 12 12];

% s0 = n*3: s0(1,:)=[p10, p20, l0]
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
% K1=1;  % D+L->DL
K2=2000; % A+L->LA
K4=10;  % DL+A->DLA
C0=10; % ratio between intensity and complex conc., proportional to FRET efficiency
coeff0=[K2 K4 C0];
%% nonlinear regression using 'nlinfit.m'
% ihetero_nlinfit(u);
options = statset('MaxIter',1000000,...
    'TolFun',1e-16,'DerivStep',1,...
    'TolX',1e-20);
mdl = NonLinearModel.fit(s0,p1lp2,@ihetero_nlinfit_fixK1,coeff0,'Options',options,'CoefficientNames',{'K2','K4','C'});
mdl

%% figure generation
coeff=mdl.Coefficients.Estimate;
X2=10.^linspace(-1,5,500)'; % for fitted line generation
for i=1:n
    YMatrix2(:,i)=ihetero_nlinfit_fixK1(total_conc(i,1),total_conc(i,2),coeff,X2);
end 

createfigure2(data_length, [s0(:,3), p1lp2], X2, YMatrix2,[],[])
end

%%
function y = ihetero_nlinfit_fixK1 (varargin)
%     K_2=1; K_1/K_2=n; 
%     c=[K1,K2,K4]
%     c=lsqcurvefit('ihetero',c0,l0,p1lp2);
%     K1, K2, K4 always >0 0<A<10
 persistent u;
switch length(varargin)
    case 1
        u=varargin{1};
        return;
    case 2
        c1 = varargin{1}; 
        s0= varargin{2};
        p10=s0(:,1);    p20=s0(:,2);    l0=s0(:,3);
    case 4
        l0 =varargin{4};
        p10=ones(size(l0)); p20=ones(size(l0));
        p10=p10.*varargin{1};
        p20=p20.*varargin{2};
        c1  =varargin{3};
end

% coefficients transformation 
% c=i_logistictransform(c1,u);
c=c1;
    K1=0.3;
    K2=c(1);
    K4=c(2);
    A=c(3);
alpha = K4./K2;

y=zeros(size(l0));
for i=1:length(l0)

    c=coeff_ihetero(alpha, K1, K2, l0(i), p10(i), p20(i));
    % check if c contains NaN or Inf
    for j=1:length(c)
        if isnan(c(j)) || c(j)==Inf;
            break;
        end
    end
    
    r=roots(c);
    l=0;
    for j=1:length(r)
        temp=p1lp2_ihetero(K1, K2, r(j), l0(i), p10(i), p20(i));
        if ((r(j)>0) && (r(j)<l0(i)) && (temp<p10(i)) && (temp<p20(i)) && (temp>0))
            l=r(j);
            y(i)=temp;
        end
    end
    if l==0
        fprintf('no real root \n');
    end
end

y=A*y;
end