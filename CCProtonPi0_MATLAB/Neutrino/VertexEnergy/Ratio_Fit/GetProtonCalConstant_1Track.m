
function [] = GetProtonCalConstant_1Track()

% Read Data
data = load('extra_energy_reco_ratio_1Track.txt');

evis = data(:,1);
avg = data(:,2);
n_events = data(:,3);

createfigure(evis,avg,n_events);
    
end

function createfigure(x, y, n_events)

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontWeight','bold','FontSize',24);
box(axes1,'on');
hold(axes1,'all');


ylim([0 3]);
xlim([0 200]);

% Create plot for x vs y
errorbar(x,y,y./sqrt(n_events),...
        'Marker','*','LineStyle','none',...
        'DisplayName','MC with Statistical Error');

hold on;

% Plot Fit
m = -0.001442;
c = 1.964;
x_fit = linspace(0,200,1000);
y_fit = m*x_fit + c;
fit_legend = sprintf('y = -(1.4E^{-3})x + %3.2f',c);
plot(x_fit,y_fit,...
    'k','LineWidth',2,...
    'DisplayName',fit_legend);


% Create xlabel
xlabel('Vertex Visible Energy [MeV]','FontWeight','bold','FontSize',24);

% Create ylabel
ylabel('T_{True}/E_{Visible}','FontWeight','bold','FontSize',24);

% Create legend
legend(axes1,'show');


end

