function [MaRES, ChiSq] = Plot_MaRES_ChiSq( )

MaRESVector = load('MaRESVector.txt');
ChiSqVector = load('ChiSqVector.txt');

MaRES = MaRESVector(:,2);
ChiSq = ChiSqVector(:,2);

plotCustom(MaRES, ChiSq);

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
    'DisplayName', '$$\chi^2$$');

xlim([min(x) max(x)]);
ylim([0.0 max(y)*1.5]);

% Draw Nominal MaRES Value
x_MaRES = linspace(1.12, 1.12, 1000);
y_MaRES = linspace(0.0, max(y), 1000);

hold on;
plot (x_MaRES,y_MaRES, 'b-',...
    'LineWidth',1,...
    'DisplayName', 'MaRES = 1.12 GeV');


% Create xlabel
xlabel('MaRES [GeV]','FontWeight','bold','FontSize', 20);

% Create ylabel
ylabel('$$\chi^2$$','Interpreter','latex','FontWeight','bold','FontSize', 20);

% Create legend
legend1 = legend(axes1,'show','Location','northeast');

set(legend1 ,'Interpreter','Latex');

print(figure1, 'Plots/MaRES_ChiSq.png','-dpng')

end

function y = GetFitResult(x, pars)
y = pars(1).*x.^3 + pars(2).*x.^2 + pars(3).*x + pars(1);
end

