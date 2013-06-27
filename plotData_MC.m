function plotData_MC(bins, data, monteCarlo, plotTitle, plotXLabel)
%  plotData_MC(bins, data, monteCarlo, plotTitle, plotXLabel) 
%   Author: Ozgur Altinok
%   Date: 2013-06-27  

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
 
% X-Axis Limits
bin_width = bins(2)-bins(1);
x_min = min(bins) - bin_width/2;
x_max = max(bins) + bin_width/2;

xlim([x_min x_max]);

% Y-Axis Limits are Automatic
% ylim([0 ylimit]);

% Calculate histograms
data = hist(data,bins);
monteCarlo = hist(monteCarlo,bins);

% Create bar for MC
bar(bins,monteCarlo,...
    'FaceColor',[0.75 0.75 0],...
    'EdgeColor',[0.75 0.75 0],...
    'BarWidth',0.8',...
    'BaseValue',0,...
    'LineStyle','-',...
    'Parent',axes1,...
    'DisplayName','Monte Carlo');

% Create bar for Data
bar(bins,data,...
    'FaceColor','none',...
    'EdgeColor','k',...
    'BarWidth',0.8',...
    'BaseValue',0,...
    'LineStyle','-',...
    'LineWidth',3.0,...
    'Parent',axes1,...
    'DisplayName','Data');

% Calculate the statistical error on Data
statError = sqrt(data + 0.25)+0.5;
errorbar(bins,data,statError,...
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


