function plotData_MC(xvector, data, monteCarlo, plotTitle, plotXLabel)
%  plotData_MC(xvector, data, monteCarlo, plotTitle, plotXLabel) 
%   Author: Ozgur Altinok
%   Date: 2013-05-30  

%%

scrsz = get(0,'ScreenSize');

% Create figure
figure1 = figure('Position',[0 0 scrsz(3)/2 scrsz(4)]);

% Create axes
axes1 = axes('Parent',figure1,'YScale','linear','YMinorTick','on',...
    'PlotBoxAspectRatio',[1 1 1],...
    'FontSize',24);
box(axes1,'on');
hold(axes1,'all');
 
xlim([-5 185]);
% ylim([0 ylimit]);

nbins = length(xvector);


% Create bar for MC
bar(xvector,monteCarlo,...
    'FaceColor',[0.75 0.75 0],...
    'EdgeColor',[0.75 0.75 0],...
    'BarWidth',0.8',...
    'BaseValue',0,...
    'LineStyle','-',...
    'Parent',axes1,...
    'DisplayName','Monte Carlo');

% Create bar for Data
bar(xvector,data,...
    'FaceColor','none',...
    'EdgeColor','k',...
    'BarWidth',0.8',...
    'BaseValue',0,...
    'LineStyle','-',...
    'LineWidth',3.0,...
    'Parent',axes1,...
    'DisplayName','Data');

statError = sqrt(data + 0.25)+0.5;
errorbar(xvector,data,statError,...
    'x','LineWidth',2.0,...
    'DisplayName','Statistical Error');


% Create xlabel
xlabel(plotXLabel,'FontSize',24);

% Create ylabel
ylabel('N(Events)','FontSize',24);

% Create title
title(plotTitle,'FontSize',24);

% Create legend
legend(axes1,'show');

end


