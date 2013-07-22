function plot1DHist(bin, data, plotTitle, plotXLabel, dataLabel)
%%  plot1DHist(bin, data, plotTitle, plotXLabel, dataLabel)
%       Author: Ozgur Altinok
%       Date: 2013-05-30    Modified: 2013-07-02


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
bin_width = bin(2)-bin(1);
x_min = min(bin) - bin_width/2;
x_max = max(bin) + bin_width/2;

xlim([x_min x_max]);


% Y-Axis Limits are Automatic
% ylim([0 ylimit]);

% Calculate histogram
data = hist(data,bin);

% Create bar for data
bar(bin,data,...
    'FaceColor',[0.75 0.75 0],...
    'EdgeColor',[0.75 0.75 0],...
    'BarWidth',0.8',...
    'BaseValue',0,...
    'LineStyle','-',...
    'Parent',axes1,...
    'DisplayName',dataLabel);

% Create xlabel
xlabel(plotXLabel,'FontSize',24);

% Create ylabel
ylabel('N(Events)','FontSize',24);

% Create title
title(plotTitle,'FontSize',24);

% Create legend
legend(axes1,'show');

end


