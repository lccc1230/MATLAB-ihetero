%% generate computational experimental data
%  generate all lines
%  use the following functions:
%   ihetero_nlinfit(P10,P20,c,L0)
%   p1lp2_ihetero(K1,K2,L,L0,P10,P20)
%   createfigure3
%
%  Author: Chang Lu (lu-c12@mails.tsinghua.edu.cn)
%  Paper: Quantitative analysis of ligand induced Hetero-dimerization
%
%
% NOTE: alpha represents the reverse of coorperativity factor used in the paper. 

clc;
clear;
%% computative generation of data, 
% assuming kineitc parameters & total conc. of receptors
% D for Donor, A for Acceptor
K1=0.3;  % D+L->DL
K2=1434.8; % A+L->LA
alpha=1/117;  % DL+A->DLA
P10=55.4;
P20=39;
pmin=min([P10 P20]);
% P0_average=(P10+P20)/2;
n=length(P10);

%% produce equilibirum data
% X2=10.^linspace(-1,4,500)'; % for fitted line generation
L0=10.^linspace(-2,7,500)';
YMatrix2=zeros(500,5);
maxL0=zeros(n,1);maxp1Lp2=zeros(n,1);
    c=[K1,K2,alpha,1];
    % calculate equilibrium conc. p1lp2 & L
    [YMatrix2(:,1),L]=ihetero_equil_conc(P10,P20,c,L0);
    YMatrix2(:,1)=YMatrix2(:,1)./pmin;
    % 
    [p1Lp2,p1,p2,p1L,Lp2] =p1lp2_ihetero(K1, K2, L, L0, P10, P20);
    YMatrix2(:,2)=p1/P10;
    YMatrix2(:,3)=p2/P20;
    YMatrix2(:,4)=p1L/P10;
    YMatrix2(:,5)=Lp2/P20;
%     YMatrix2(:,2)=0;
%     YMatrix2(:,3)=0;
%     YMatrix2(:,4)=0;
%     YMatrix2(:,5)=0;
    % [maxL0(i), maxp1Lp2(i)]  = maxp1Lp2_ihetero(P10,P20,c);

%% figure
% createfigure1(n,[], X2, YMatrix2,maxL0,maxp1Lp2);
createfigure3(L0, YMatrix2,[],[]);
%% print result
% print -r300 -depsc /Users/luc/Documents/Projects/sample_figures/