function createfigure1(data_length, ~, X2, YMatrix2, X3, Y3)
%CREATEFIGURE2(X1,YMATRIX1,X2,YMATRIX2)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data
%  X2:  vector of x data
%  YMATRIX2:  matrix of y data
%
%
%  Author: Chang Lu (lu-c12@mails.tsinghua.edu.cn)
%  Paper: Quantitative analysis of ligand induced Hetero-dimerization

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'XScale','log','XMinorTick','on',...
    'Position',[0.13 0.141730103806228 0.775 0.815],...
    'LineWidth',2,...
    'FontSize',16,'FontName','Arial');
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[-0.5 100000]);
% ylim(axes1,[0 1.2]);
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to semilogx
% count=0;
colo=hsv(data_length);
% colo=zeros(data_length, 3);
% colo(1,:)=[0 0 0]; colo(2,:)=[0 1 0]; colo(3,:)=[1 0 0];
% colo(4,:)=[0 1 0]; colo(5,:)=[1 0 0];
for i=1:data_length
%     if ~isempty(data)
%         X1=data(count+1:count+data_length(i),1);
%         Y1=data(count+1:count+data_length(i),2);
%         count=count+data_length(i);
%         semilogx(X1,Y1,'MarkerSize',8,'Marker','x','LineWidth',3,...
%         'LineStyle','none',...
%         'Color',colo(i,:),...
%         'Parent',axes1);
%     end
    % Create multiple lines using matrix input to semilogx
    Y2=YMatrix2(:,i);
    semilogx(X2,Y2,'LineWidth',3,'Color',colo(i,:),...
    'Parent',axes1);
end
% set(semilogx1(1),'MarkerSize',15);

% Create xlabel
xlabel('Inducer (nM)','Interpreter','latex','FontSize',24);

% Create ylabel
ylabel('Hetero-dimer','Interpreter','latex','FontSize',24);

% Create semilogx
semilogx(X3,Y3,'MarkerSize',15,'Marker','x','LineWidth',6,...
    'LineStyle','none',...
    'Color',[1 0 0]);
% grid on;
end
