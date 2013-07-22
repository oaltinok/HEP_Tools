function plotVS(xData, yData, plotTitle, plotXLabel, plotYLabel)
%%  plotVS(xData, yData, plotTitle, plotXLabel, plotYLabel)
%       Author: Ozgur Altinok
%       Date: 2013-05-30 

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

% X-Axis Limits are Automatic
% xlim([x_min x_max]);


% Y-Axis Limits are Automatic
% ylim([0 ylimit]);

plot(xData,yData,'.');

% Create xlabel
xlabel(plotXLabel,'FontSize',24);

% Create ylabel
ylabel(plotYLabel,'FontSize',24);

% Create title
title(plotTitle,'FontSize',24);


end


