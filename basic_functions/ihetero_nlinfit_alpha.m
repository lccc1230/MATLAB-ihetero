%     K_2=1; K_1/K_2=n; 
%     c=[K1,K2,alpha]
%     c=lsqcurvefit('ihetero',c0,l0,p1lp2);
%     K1, K2, K4 always >0 0<A<10
%   used functions: 
%       coeff_ihetero(alpha, K1, K2, l0(i), p10(i), p20(i));
%       p1lp2_ihetero(K1, K2, r(j), l0(i), p10(i), p20(i));
%
%
%  Author: Chang Lu (lu-c12@mails.tsinghua.edu.cn)
%  Paper: Quantitative analysis of ligand induced Hetero-dimerization

function [y,l_free] = ihetero_nlinfit_alpha (varargin)

switch length(varargin)
    case 2
        c = varargin{1}; 
        s0= varargin{2};
        p10=s0(:,1);    p20=s0(:,2);    l0=s0(:,3);
    case 4
        l0 =varargin{4};
        p10=ones(size(l0)); p20=ones(size(l0));
        p10=p10.*varargin{1};
        p20=p20.*varargin{2};
        c  =varargin{3};
end

    K1=c(1);
    K2=c(2);
    alpha=c(3);
    A=c(4);
% alpha = K4./K2;

y=zeros(size(l0));
l_free=zeros(size(l0));
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
        if ((r(j)>0) && (r(j)<l0(i)))
            temp=p1lp2_ihetero(K1, K2, r(j), l0(i), p10(i), p20(i));
            if ( (temp<p10(i)) && (temp<p20(i)) && (temp>0))
                y(i)=temp;
                l_free(i)=r(j);
                l=l+1;
            end
        end
    end
    if l==0
        fprintf('no real root \n');
    elseif l>1
        fprintf('more than one root \n');
    end
end

y=A*y;
end