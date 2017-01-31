function [] = Plot_Significance( )

MaRESVector = load('MaRESVector.txt');
ChiSqVector = load('ChiSqVector.txt');

MaRES = MaRESVector(:,2);
ChiSq = ChiSqVector(:,2);
DeltaChiSq = ChiSq - min(ChiSq);
significance = sqrt(DeltaChiSq);
plotCustom(MaRES, significance, 'Plots/Significance.png', '$$\sigma = \sqrt{\Delta \chi^2}$$');

end

function [] = plotCustom(x, y, plot_name, y_label)

% Create figure
figure1 = figure('visible','off', 'position', [0, 0, 1280,800]);

% Create axes
axes1 = axes('Parent',figure1,'FontSize',16);
box(axes1,'on');
hold(axes1,'all');

% Plot Data Points
plot (x,y, 'o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'MarkerSize',3,...
    'DisplayName','$$\sigma = \sqrt{\Delta \chi^2}$$');

xlim([min(x) max(x)]);
ylim([0.0 10]);

% Draw Sigma Lines
x_sigma = linspace(min(x), max(x), 1000);
y_3sigma = linspace(3, 3, 1000);
y_5sigma = linspace(5, 5, 1000);

hold on;
plot (x_sigma,y_5sigma, 'r-',...
    'LineWidth',2,...
    'DisplayName', '$$5\sigma$$');

plot (x_sigma,y_3sigma, 'b-',...
    'LineWidth',2,...
    'DisplayName', '$$3\sigma$$');


% Draw Nominal MaRES Value
x_MaRES = linspace(1.12, 1.12, 1000);
y_MaRES = linspace(0.0, max(y)*1.25, 1000);

plot (x_MaRES,y_MaRES, 'b--',...
    'LineWidth',1,...
    'DisplayName', 'MaRES = 1.12 GeV');


% Create xlabel
xlabel('MaRES [GeV]','FontWeight','bold','FontSize', 20);

% Create ylabel
ylabel(y_label,'Interpreter','latex','FontWeight','bold','FontSize', 20);

% Create legend
legend1 = legend(axes1,'show','Location','northwest');

set(legend1 ,'Interpreter','Latex');

print(figure1, plot_name,'-dpng')

end


