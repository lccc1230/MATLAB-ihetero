
%%  notes
%   k1, k2, k3 always >0 0<A<10 ???
%   used functions: 
%       coeff_ihetero(alpha, k1, k2, l0(i), p10(i), p20(i)); (generate fifth order equation parameters)
%       p1lp2_ihetero(k1, k2, r(j), l0(i), p10(i), p20(i)); (calculate p1lp2 given free ligand and coefficients)
%  possible errors:
%  ROOT ERROR!!! no reasonable root
%  ROOT ERROR!!! more than one reasonable root
%
%
%  Author: Chang Lu (lu-c12@mails.tsinghua.edu.cn)
%  Paper: Quantitative analysis of ligand induced Hetero-dimerization

%%
function [y,l_free] = ihetero_equil_conc (varargin)

switch length(varargin)
    case 2
        % use data points {l0, p1lp2}
        % to fit coefficients c using initial parameters c0 
        % used with function lsqcurvefit: 
        %   c=[k1,k2,alpha]
        %   c=lsqcurvefit('ihetero',c0,l0,p1lp2);
        c = varargin{1}; 
        s0= varargin{2};
        p10=s0(:,1);    p20=s0(:,2);    l0=s0(:,3);
    case 4
        % to calculate equilibrium concentrations of each species 
        % with kinetic constants c and total conc. l0, p10, p20
        l0 =varargin{4};
        p10=ones(size(l0)); p20=ones(size(l0));
        p10=p10.*varargin{1};
        p20=p20.*varargin{2};
        c  =varargin{3};
end

    k1=c(1);
    k2=c(2);
    alpha=c(3);
    A=c(4);
% alpha = k3./k2;

y=zeros(size(l0));
l_free=zeros(size(l0));
for i=1:length(l0)

    c=coeff_ihetero(alpha, k1, k2, l0(i), p10(i), p20(i));
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
            temp=p1lp2_ihetero(k1, k2, r(j), l0(i), p10(i), p20(i));
            if ( (temp<p10(i)) && (temp<p20(i)) && (temp>0))
                y(i)=temp;
                l_free(i)=r(j);
                l=l+1;
            end
        end
    end
    if l==0
        fprintf('ROOT ERROR!!! no reasonable root \n');
    elseif l>1
        fprintf('ROOT ERROR!!! more than one reasonable root \n');
    end
end

y=A*y;
end