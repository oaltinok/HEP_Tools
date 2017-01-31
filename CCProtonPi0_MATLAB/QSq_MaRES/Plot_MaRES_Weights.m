function [] = Plot_MaRES_Weights( )

MaRESVector = load('MaRESVector.txt');
WeightsVector = load('Weights.txt');

MaRES = MaRESVector(:,2);
Weights= WeightsVector(:,2);

plotCustom(MaRES, Weights);

end

function [] = plotCustom(x, y)

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
    'MarkerSize',3, ...
    'DisplayName', 'Weights');

xlim([min(x) max(x)]);
ylim([0.6 1.4]);

% Draw 1 Line
x_1 = linspace(min(x), max(x), 1000);
y_1 = linspace(1.0, 1.0, 1000);

hold on;
plot (x_1 ,y_1, 'b--',...
    'LineWidth',1,...
    'DisplayName', '$$w = 1.0$$');

% Create xlabel
xlabel('MaRES [GeV]','FontWeight','bold','FontSize', 20);

% Create ylabel
ylabel('Weights','FontWeight','bold','FontSize', 20);

% Create legend
legend1 = legend(axes1,'show','Location','northeast');

set(legend1 ,'Interpreter','Latex');

print(figure1, 'Plots/MaRES_Weights.png','-dpng')

end


