function [ ] = PlotFactor()

A = 1.0;
Q0 = 0.114;

plotCustom(A, Q0);

end

function y = CalcFactor(x, A, Q0)

y = A ./ (1 + exp(1 - sqrt(x)/Q0));

end

function [] = plotCustom(A, Q0)

% Create figure
figure1 = figure('visible','off', 'position', [0, 0, 1280,800]);

% Create axes
axes1 = axes('Parent',figure1,'FontSize',16);
box(axes1,'on');
hold(axes1,'all');

QSq = linspace(0, 2, 1000);
Factor_MINOS = CalcFactor(QSq, 1.010, 0.156);

% Plot MINOS Factor
plot (QSq,Factor_MINOS, 'o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'MarkerSize',3,...
    'DisplayName','$$\Delta$$ Factor: MINOS');

xlim([min(QSq) max(QSq)]);
ylim([0.0 1.2]);

hold on;

Factor_Test = CalcFactor(QSq, A, Q0);

plot (QSq,Factor_Test , 'o',...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r',...
    'MarkerSize',3,...
    'DisplayName','$$\Delta$$ Factor: Best');

% Create xlabel
xlabel('Q^2 [GeV^2]','FontWeight','bold','FontSize', 20);

% Create ylabel
ylabel('\Delta Suppression Factor','FontWeight','bold','FontSize', 20);

% Create legend
legend1 = legend(axes1,'show','Location','northeast');

set(legend1 ,'Interpreter','Latex');

print(figure1, 'Suppressian_Factor.png','-dpng')

end
