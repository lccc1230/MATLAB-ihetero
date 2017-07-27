% Calculate maximum P1LP2 concentration
%
%
%  Author: Chang Lu (lu-c12@mails.tsinghua.edu.cn)
%  Paper: Quantitative analysis of ligand induced Hetero-dimerization


function [maxl0, maxp1lp2] = maxp1lp2_ihetero( p10, p20,c)
% u=[10,30000,1000,20]; % upper bound on parameters

% c=i_logistictransform(c,u);
K1=c(1);
K2=c(2);
alpha=c(3);

maxl0=sqrt(K1.*K2)+(sqrt(K1).*p20+sqrt(K2).*p10)./(sqrt(K1)+sqrt(K2));

% alpha=K4./K2;

tmp = p10+p20+alpha.*(sqrt(K1)+sqrt(K2)).^2;

delta =  tmp.^2 - 4.*p10.*p20;

maxp1lp2 = (tmp - sqrt(delta))./2.*c(4);