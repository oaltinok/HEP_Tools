function [  ] = GetProtonCalConstant( )

%  GetCalibrationConstant('PC_theta_00.txt');
%  GetCalibrationConstant('PC_theta_60.txt');
%  GetCalibrationConstant('PC_short_protons_theta_all.txt');
  GetCalibrationConstant('CCProtonPi0_protons.txt');

end

function [] = GetCalibrationConstant(input_file)

% Read Data
data = load(input_file);

evis = data(:,1);
avg = data(:,2);
n_events = data(:,3);

% % Remove First Bin -- Too much uncertainity
% evis = evis(2:end);
% avg = avg(2:end);
% n_events = n_events(2:end);

% Select Only Events having more than 50 entries
n_events_high = n_events>50;

evis = evis(n_events_high);
avg = avg(n_events_high);
n_events = n_events(n_events_high);

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


ylim([0 2]);
xlim([0 450]);

% Create plot for x vs y
errorbar(x,y,y./sqrt(n_events),...
        'Marker','*','LineStyle','none',...
        'DisplayName','MC with Statistical Error');


hold on;


x_mean = linspace(0,450,1000);
y_mean = linspace(weighted_mean, weighted_mean,1000);
plot_legend = sprintf('%s%3.2f','Weighted Mean = ',weighted_mean);

plot(x_mean,y_mean,...
    'r','LineWidth',2,...
    'DisplayName',plot_legend );

% Plot Fit if there is
m = -0.003856;
c = 1.415;
x_fit = linspace(0,200,1000);
y_fit = m*x_fit + c;
fit_legend = sprintf('y = -(3.9E^{-3})x + %3.2f',c);
plot(x_fit,y_fit,...
    'k','LineWidth',2,...
    'DisplayName',fit_legend);


% Create xlabel
xlabel('Visible Energy [MeV]','FontWeight','bold','FontSize',24);

% Create ylabel
ylabel('T_{True}/E_{Visible}','FontWeight','bold','FontSize',24);

% Create legend
legend(axes1,'show');


end

