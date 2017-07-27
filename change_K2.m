%% generate computational experimental data - change K2
% use the following functions:
%   ihetero_equil_conc(P10,P20,c,L0)
%   maxp1lp2_ihetero(P10,P20, [K1,K2,k3/K2,1])
%   createfigure1
%
%
%  Author: Chang Lu (lu-c12@mails.tsinghua.edu.cn)
%  Paper: Quantitative analysis of ligand induced Hetero-dimerization
%
%
% NOTE: alpha represents the reverse of coorperativity factor used in the paper. 

clc;
clear;
%% computation generation of data
% assuming kinetic constants & total conc.
% K1    % D+L->DL
% K2    % A+L->LA
% k3    % D+LA->DLA
% k4;   % DL+A->DLA
% alpha=k3/K1=k4/K2
K1=10;
% K2=10.^linspace(0,4,500)';
K2=[1,10,100,1000];
alpha=1;
P10=10;
P20=200;

n=length(K2);
% n=1;

%% produce equilibrium data
% X2=10.^linspace(-1,4,500)'; % for fitted line generation
L0=10.^linspace(-2,5,500)'; % for fitted line generation
YMatrix2=zeros(500,n);
maxL0=zeros(n,1);maxp1lp2=zeros(n,1);
EC50=zeros(n,1);IC50=zeros(n,1);EDR=zeros(n,1);
for i=1:n
    c=[K1,K2(i),alpha,1];
    % calculate equilibrium conc. p1lp2 & L
    YMatrix2(:,i)=ihetero_equil_conc(P10,P20,c,L0);
%     YMatrix2(:,i)=ihetero_equil_conc_ligand_in_excess(P10,P20,c,L0);

    % calculate maximum p1lp2
    [maxL0(i), maxp1lp2(i)]  = maxp1lp2_ihetero(P10,P20,[K1,K2(i),alpha,1]);
    % calculate EDR
    [EC50(i), IC50(i), EDR(i)] = EDR_ihetero (maxp1lp2(i), K1, K2(i), alpha, P10, P20);
    % gui yi hua
%     YMatrix2(:,i)=YMatrix2(:,i)/maxp1lp2(i);
    YMatrix2(:,i)=YMatrix2(:,i)/P10;
end 

createfigure1(n,[], L0, YMatrix2,[],[]);
hold on;
plot(maxL0,maxp1lp2/P10,'MarkerSize',12,'Marker','x','LineWidth',2,...
    'LineStyle','none',...
    'Color',[0 0 0]);
