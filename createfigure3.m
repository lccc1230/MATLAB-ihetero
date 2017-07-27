% createfigure
%
%  Author: Chang Lu (lu-c12@mails.tsinghua.edu.cn)
%  Paper: Quantitative analysis of ligand induced Hetero-dimerization

function createfigure3(X1, YMatrix1, X2, Y2)
%CREATEFIGURE2(X1,YMATRIX1,X2,YMATRIX2)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data
%  X2:  vector of x data
%  YMATRIX2:  matrix of y data

%  Auto-generated by MATLAB on 25-Nov-2013 15:58:52

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'XScale','log','XMinorTick','on',...
    'Position',[0.13 0.141730103806228 0.775 0.815],...
    'LineWidth',2,...
    'FontSize',16,'FontName','Arial');
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[-1 1000000]);
% ylim(axes1,[0 1.2]);
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to semilogx
semilogx1 = semilogx(X1,YMatrix1,'LineWidth',3,'Parent',axes1);
set(semilogx1(1),'Color',[0 0 0]);
set(semilogx1(2),'LineStyle','--','Color',[0 1 0]);
set(semilogx1(3),'LineStyle','--','Color',[1 0 0]);
set(semilogx1(4),'LineStyle',':','Color',[0 1 0]);
set(semilogx1(5),'LineStyle',':','Color',[1 0 0]);

% Create xlabel
xlabel('Inducer (nM)','Interpreter','latex','FontSize',24);

% Create ylabel
ylabel('Hetero-dimer','Interpreter','latex','FontSize',24);

% Create semilogx
semilogx(X2,Y2,'MarkerSize',15,'Marker','x','LineWidth',6,...
    'LineStyle','none',...
    'Color',[1 0 0]);
% grid on;
end