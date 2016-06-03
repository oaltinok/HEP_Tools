function [] = ComparisonPlots( )

data = load('EventKinematics.txt');

[Enu_reco, Enu_truth, Enu_genie] = GetVariable(data, 1, 2, 3);
[QSq_reco, QSq_truth, QSq_genie] = GetVariable(data, 4, 5, 6);
[W_reco, W_truth, W_genie] = GetVariable(data, 7, 8, 9);

PlotVariable(Enu_truth, Enu_genie, 20, 'E_{\nu} from Truth [GeV]','E_{\nu} from GENIE [GeV]','Enu_truth_genie.png');
PlotVariable(Enu_reco, Enu_truth, 20, 'E_{\nu} from Reco [GeV]','E_{\nu} from Truth [GeV]','Enu_reco_truth.png');
PlotVariable(Enu_reco, Enu_genie, 20, 'E_{\nu} from Reco [GeV]','E_{\nu} from GENIE [GeV]','Enu_reco_genie.png');

PlotVariable(QSq_truth, QSq_genie, 2, 'Q^{2} from Truth [GeV^{2}]','Q^{2} from GENIE [GeV^{2}]','QSq_truth_genie.png');
PlotVariable(QSq_reco, QSq_truth, 2, 'Q^{2} from Reco [GeV^{2}]','Q^{2} from Truth [GeV^{2}]','QSq_reco_truth.png');
PlotVariable(QSq_reco, QSq_genie, 2, 'Q^{2} from Reco [GeV^{2}]','Q^{2} from GENIE [GeV^{2}]','QSq_reco_genie.png');

PlotVariable(W_truth, W_genie, 3, 'W from Truth [GeV]','W from GENIE [GeV]','W_truth_genie.png');
PlotVariable(W_reco, W_truth, 3, 'W from Reco [GeV]','W from Truth [GeV]','W_reco_truth.png');
PlotVariable(W_reco, W_genie, 3, 'W from Reco [GeV]','W from GENIE [GeV]','W_reco_genie.png');

end

function [] = PlotVariable(x, y, max_x, x_label, y_label, fig_name)
% Create figure
figure1 = figure('visible','off');

% Create axes
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontWeight','bold','FontSize',20);
box(axes1,'on');
hold(axes1,'all');


% Create plot for evis vs true
plot(x,y,...
        'Marker','*','LineStyle','none',...
        'DisplayName','MC Points');

xlim([0 max_x]);
ylim([0 max_x]);

hold on;

% x = y Line
m_x = linspace(0,max_x,1000);
m_y = m_x;
plot(m_x,m_y,...
    'k','LineWidth',2,...
    'DisplayName','y = x line');

% Create xlabel
xlabel(x_label,'FontWeight','bold','FontSize', 20);

% Create ylabel
ylabel(y_label,'FontWeight','bold','FontSize', 20);

% Create legend
legend(axes1,'show','Location','northwest');

print(figure1,fig_name,'-dpng')

end



function [reco, truth, genie] = GetVariable(data, reco_ind, truth_ind, genie_ind)

reco = data(:,reco_ind);
truth = data(:,truth_ind);
genie = data(:,genie_ind);

end

