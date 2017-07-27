% Generate experimental data using coefficients and total concentrations of
% ligand and two proteins
% 
% Generated from paper by Douglass et. al. (JACS, 2013)

function p1lp2 = generate_data_Douglass(K1,K2,K4,l0,p10,p20)

alpha=K2./K4;

count=zeros(length(l0),1);
p1lp2=zeros(length(l0),1);

%tic;
for i=1:length(l0)

    c=coeff_ihetero_Douglass(alpha, K1, K2, l0(i), p10, p20);
    r=roots(c);
    for j=1:length(r)
        if ((r(j)>0) && (r(j)<l0(i)) && (r(j)<p10) && (r(j)<p20) )
            p1lp2(i)=r(j);
%             p1lp2(i)=r(j)*(1+0.1*rand(1)); % add error
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

