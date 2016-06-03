function [ ] = QSquarePlots(  )

data = load('QSquare.txt');

res_events = data(:,1) == 2;
dis_events = data(:,1) == 3;

QSq_res_truth = data(res_events, 2);
QSq_res_genie = data(res_events, 3);

QSq_dis_truth = data(dis_events, 2);
QSq_dis_genie = data(dis_events, 3);

PlotVariable(QSq_res_truth, QSq_res_genie, QSq_dis_truth, QSq_dis_genie,...
    'Q^{2} from Truth [GeV^{2}]','Q^{2} from GENIE [GeV^{2}]', 'QSquare_Type.png')

end

function [] = PlotVariable(res_x, res_y, dis_x, dis_y, x_label, y_label, fig_name)
% Create figure
figure1 = figure('visible','off');

% Create axes
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontWeight','bold','FontSize',20);
box(axes1,'on');
hold(axes1,'all');


% Create plot for evis vs true
plot(res_x,res_y,...
        'Marker','*',...
        'Color','b',...
        'LineStyle','none',...
        'DisplayName','Resonance');
    
xlim([0 2]);
ylim([0 2]);
    
hold on;

plot(dis_x,dis_y,...
        'Marker','*',...
        'Color','r',...
        'LineStyle','none',...
        'DisplayName','DIS');

% x = y Line
m_x = linspace(0,2,1000);
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

