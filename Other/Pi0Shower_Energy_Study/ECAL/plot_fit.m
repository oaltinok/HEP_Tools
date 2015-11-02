function [] = plot_fit()

[count_array,evis_array,mean_array,min_array,max_array,weights] = GetValues()

%% Plot Only Fit Region
% Create figure
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',24,'FontWeight','bold');
ylim([1.0 6.0]);
hold on;

% Plot Low Energy Points
% errorbar(evis_array,mean_array,(mean_array-min_array),(max_array-mean_array),...
%     'MarkerSize',10,'MarkerFaceColor',[0 0 0],...
%     'MarkerEdgeColor',[0 0 0],...
%     'Marker','o',...
%     'LineStyle','none',...
%     'LineWidth',1,...
%     'DisplayName','MC Points with HWHM',...
%     'Parent',axes1,...
%     'Color',[0 0 0]);

errorbar(evis_array,mean_array,mean_array./sqrt(count_array),...
    'MarkerSize',10,'MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor',[0 0 0],...
    'Marker','o',...
    'LineStyle','none',...
    'LineWidth',1,...
    'DisplayName','MC Points with Statistical Error',...
    'Parent',axes1,...
    'Color',[0 0 0]);

% Plot Low Energy Fit
x = linspace(0,300,1000);
y = linspace(3.10,3.10,1000);
plot(x,y,'LineWidth',2,'Color','r','DisplayName','Jaewons kE = 3.10');

% Plot Low Energy Fit
x = linspace(0,300,1000);
y = linspace(3.2367,3.2367,1000);
plot(x,y,'LineWidth',2,'Color','g','DisplayName','Weighted Mean kE = 3.2367');

% Create xlabel
xlabel('Visible Energy [MeV]','FontSize',24,'FontWeight','bold');
% Create ylabel
ylabel('E_{True}/E_{Visible} (kE)','FontSize',24,'FontWeight','bold');
% Create legend
legend(axes1,'show');

hold off;

end

