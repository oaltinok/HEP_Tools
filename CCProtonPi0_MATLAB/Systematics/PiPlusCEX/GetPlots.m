function [ ] = GetPlots( xsec )
mb_unit = 1.0*10^27;
N_C = N_CarbonAtoms();
N_C = N_C / mb_unit;


length = linspace(0,100,1000);
Prob = 1-exp(-length*xsec*N_C);
Prob_up = 1-exp(-length*(xsec*1.5)*N_C);
Prob_down = 1-exp(-length*(xsec*0.5)*N_C);

probs = [Prob;Prob_up;Prob_down];
wgts = [Prob_up./Prob;Prob_down./Prob];

PlotProbabilities(length, probs);
PlotWeights(length, wgts);


end

function PlotProbabilities(X1, YMatrix1)
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',2,'Parent',axes1);
set(plot1(1),'DisplayName','Central Value','Color',[0 0 0]);
set(plot1(2),'DisplayName','50% Increased XSec','Color',[0 0 1]);
set(plot1(3),'DisplayName','50% Decreased XSec','Color',[1 0 0]);

% Create xlabel
xlabel('Track Length [cm]');

% Create ylabel
ylabel('Probability');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',20,'FontWeight','bold');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','northwest');

end


function PlotWeights(X1, YMatrix1)
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',2,'Parent',axes1);
set(plot1(1),'DisplayName','50% Increased XSec','Color',[0 0 1]);
set(plot1(2),'DisplayName','50% Decreased XSec','Color',[1 0 0]);

ylim([0 2]);
% Create xlabel
xlabel('Track Length [cm]');

% Create ylabel
ylabel('Weight');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',20,'FontWeight','bold');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','northwest');

end


