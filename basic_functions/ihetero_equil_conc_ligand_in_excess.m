
%%  notes
%   K1, K2, K4 always >0 0<A<10 ???
%   used functions: 
%       coeff_ihetero(alpha, K1, K2, l0(i), p10(i), p20(i)); (generate fifth order equation parameters)
%       p1lp2_ihetero(K1, K2, r(j), l0(i), p10(i), p20(i)); (calculate p1lp2 given free ligand and coefficients)
%  possible errors:
%  ROOT ERROR!!! no reasonable root
%  ROOT ERROR!!! more than one reasonable root
%
%
%  Author: Chang Lu (lu-c12@mails.tsinghua.edu.cn)
%  Paper: Quantitative analysis of ligand induced Hetero-dimerization

%%
function [y] = ihetero_equil_conc_ligand_in_excess (varargin)

switch length(varargin)
%     case 2
%         % use data points {l0, p1lp2}
%         % to fit coefficients c using initial parameters c0 
%         % used with function lsqcurvefit: 
%         %   c=[K1,K2,alpha]
%         %   c=lsqcurvefit('ihetero',c0,l0,p1lp2);
%         c = varargin{1}; 
%         s0= varargin{2};
%         p10=s0(:,1);    p20=s0(:,2);    l0=s0(:,3);
    case 4
        % to calculate equilibrium concentrations of each species 
        % with kinetic constants c and total conc. l0, p10, p20
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

for i=1:length(l0)

    c1=1; c3=p10(i).*p20(i);
    c2=p10(i)+p20(i)+(l0(i)+K1).*(l0(i)+K2).*alpha./l0(i);
    
    r=roots([c1 -c2 c3]);
    if r(1)>0 && r(1)<p10(i) && r(1)<p20(i)
        y(i)=r(1);
    else
        y(i)=r(2);
    end

end

y=A*y;
end