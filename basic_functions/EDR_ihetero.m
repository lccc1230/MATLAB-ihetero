% Effective dose range
%
%  Author: Chang Lu (lu-c12@mails.tsinghua.edu.cn)
%  Paper: Quantitative analysis of ligand induced Hetero-dimerization

function [EC50, IC50, EDR] = EDR_ihetero (P1LP2max, K1, K2, a, P10, P20)

c=zeros(1,3);
p50=P1LP2max./2;

c(1) = 1;
c(2) = K1+K2-(P10-p50).*(P20-p50)./a./p50;
c(3) = K1.*K2;

r=roots(c);
if r(1)<r(2)
    L1=r(1);
    L2=r(2);
else
    L1=r(2);
    L2=r(1);
end

EC50 = L1+P20.*L1./(L1+K2)+P10.*L1./(L1+K1)+p50.*(K1.*K2-L1.^2)./(L1+K1)./(L1+K2);

IC50 = L2+P20.*L2./(L2+K2)+P10.*L2./(L2+K1)+p50.*(K1.*K2-L2.^2)./(L2+K1)./(L2+K2);

EDR = IC50./EC50;

