% Generate experimental data using coefficients and total concentrations of
% ligand and two proteins
% 
% --- for different total ligand concentrations and fixed protein
% concentration
%
%
%  Author: Chang Lu (lu-c12@mails.tsinghua.edu.cn)
%  Paper: Quantitative analysis of ligand induced Hetero-dimerization

function p1lp2 = generate_data(K1,K2,K4,l0,p10,p20)

alpha=K4./K2;

count=zeros(length(l0),1);
p1lp2=zeros(length(l0),1);

%tic;
for i=1:length(l0)

    c=coeff_ihetero(alpha, K1, K2, l0(i), p10, p20);
    r=roots(c);
    for j=1:length(r)
        temp=p1lp2_ihetero(K1, K2, r(j), l0(i), p10, p20);
        if ((r(j)>0) && (r(j)<l0(i)) && (temp<p10) && (temp<p20) && (temp>0)&& (temp<l0(i)))
%             p1lp2(i)=temp;
            p1lp2(i)=temp*(1+0.1*rand(1)); % add error
            count(i)=count(i)+1;
        end
    end
    if count(i)==0
        fprintf('no real root \n');
    elseif count(i)>1
        fprintf('more than one root \n');
    end
end
% toc;

