
function [evis,true,n_events] = GetTrueRecoFit_2Track()

% Read Data
data = load('extra_energy_reco_true_2Track.txt');

evis = data(:,1);
true = data(:,2);
n_events = data(:,3);

% createfigure(evis,true,n_events);
    
end

function createfigure(evis, true, n_events)

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontWeight','bold','FontSize',24);
box(axes1,'on');
hold(axes1,'all');


ylim([0 500]);
xlim([0 500]);


% Create plot for evis vs true
plot(evis,true,...
        'Marker','*','LineStyle','none',...
        'DisplayName','MC Points');

hold on;

% x = y Line
x = linspace(0,500,1000);
y = x;
plot(x,y,...
    'k','LineWidth',2,...
    'DisplayName','y = x line');

% Fit Line
p1 = 1.01;
p2 = 80.9;
x_fit = linspace(0,500,1000);
y_fit = p1*x + p2;
fit_legend = sprintf('y = %3.2fx + %3.2 line',p1,p2);
plot(x_fit,y_fit,...
    'r','LineWidth',2,...
    'DisplayName',fit_legend);

% Create xlabel
xlabel('Vertex Visible Energy [MeV]','FontWeight','bold','FontSize',24);

% Create ylabel
ylabel('E_{\nu} - (E_{\mu} + E_{\pi^{0}})','FontWeight','bold','FontSize',24);

% Create legend
legend(axes1,'show');


end

