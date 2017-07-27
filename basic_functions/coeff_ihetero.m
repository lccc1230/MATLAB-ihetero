%  Calculate coefficents of quintic equation
%
%  Author: Chang Lu (lu-c12@mails.tsinghua.edu.cn)
%  Paper: Quantitative analysis of ligand induced Hetero-dimerization

function c = coeff_ihetero(alpha, K1, K2, l0, p10, p20)

% NOTE: alpha represents the reverse of coorperativity factor used in the paper. 

a0 = 1-alpha;
a1 = -l0+(1-alpha).*(K1+K2-l0+p10+p20);
a2 = (K1-l0+p10).*(K2-l0+p20)-(1-alpha).*(K2.*(l0-p10)+K1.*(l0-p20));
a3 = K2.*(K1-l0+p10).*(K1.*alpha-l0+p10)+K1.*(K2-l0+p20).*(K2.*alpha-l0+p20)+2.*K1.*K2.*alpha.*l0;
a4 = K1.*K2.*(K1.*alpha-l0+p10).*(K2.*alpha-l0+p20)+(K1.*K2).^2.*alpha.*(1-alpha);
a5 = -alpha.*l0.*K1.^2.*K2.^2;


c = [a0 a1 a2 a3 a4 a5];