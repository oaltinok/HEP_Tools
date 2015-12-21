function [] = GetProtonCalConstant_2Track()

% Read Data
data = load('extra_energy_reco_ratio_2Track.txt');

evis = data(:,1);
avg = data(:,2);
n_events = data(:,3);

weighted_mean = CalcWeightedMean(avg,n_events);

createfigure(evis,avg,n_events,weighted_mean);
    
end

function [weighted_mean] = CalcWeightedMean(avg, n_events)

total_n_events = sum(n_events);
weighted_avgs = avg.*n_events;


weighted_mean = sum(weighted_avgs) ./ total_n_events;

end

function createfigure(x, y, n_events, weighted_mean)

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontWeight','bold','FontSize',24);
box(axes1,'on');
hold(axes1,'all');


ylim([0 3]);
xlim([0 100]);

% Create plot for x vs y
errorbar(x,y,y./sqrt(n_events),...
        'Marker','*','LineStyle','none',...
        'DisplayName','MC with Statistical Error');

hold on;

% Plot Mean
x_mean = linspace(0,100,1000);
y_mean = linspace(weighted_mean, weighted_mean,1000);
plot_legend = sprintf('%s%3.2f','Weighted Mean = ',weighted_mean);

plot(x_mean,y_mean,...
    'r','LineWidth',2,...
    'DisplayName',plot_legend );


% Create xlabel
xlabel('Vertex Visible Energy [MeV]','FontWeight','bold','FontSize',24);

% Create ylabel
ylabel('T_{True}/E_{Visible}','FontWeight','bold','FontSize',24);

% Create legend
legend(axes1,'show');


end

