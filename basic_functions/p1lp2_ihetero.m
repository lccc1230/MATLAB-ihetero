% Calculate equilibrium concentrations using coefficients and knowledge of
% free ligand [L] at equilibrium
%
%
%  Author: Chang Lu (lu-c12@mails.tsinghua.edu.cn)
%  Paper: Quantitative analysis of ligand induced Hetero-dimerization

function [p1lp2,p1,p2,p1l,lp2] =p1lp2_ihetero(K1, K2, l, l0, p10, p20)

p1lp2 = p20+(K2+l)./(K1.*K2-l.^2).*(l.*(l0-l-p10)+K1.*(l0-l-p20));

p1 = (l.*(l-l0+p20)+K2.*(l-l0+p10)).*K1./(K1.*K2-l.^2);

p2 = (l.*(l-l0+p10)+K1.*(l-l0+p20)).*K2./(K1.*K2-l.^2);

p1l = p1.*l./K1;

lp2 = p2.*l./K2;