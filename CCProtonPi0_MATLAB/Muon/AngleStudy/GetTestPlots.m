function [ ] = GetTestPlots()

%% Read File and Get theta vectors
all_angles = load('MuonTheta_Test.txt');

reco_theta = all_angles(:,1);
true_theta = all_angles(:,2);
calc_theta = all_angles(:,3);


%% Plot

% Create figure
figure1 = figure;

subplot1 = subplot(2,2,1,'Parent',figure1,'FontSize',24);
box(subplot1,'on');
hold(subplot1,'all');
plot_comparison(subplot1, calc_theta, true_theta, 'Calculated \theta', 'minervaCoordSysTool \theta', 'True vs Calculated \theta wrt Beam');

subplot2 = subplot(2,2,2,'Parent',figure1,'FontSize',24);
box(subplot2,'on');
hold(subplot2,'all');
plot_comparison(subplot2, reco_theta, calc_theta, 'Reconstructed \theta', 'Calculated \theta', 'Calculated vs Reco \theta wrt Beam');

subplot3 = subplot(2,2,3,'Parent',figure1,'FontSize',24);
box(subplot3,'on');
hold(subplot3,'all');
plot_error(subplot3, calc_theta, true_theta, 'Calculated \theta', '(Calc-minervaCoordSysTool)/minervaCoordSysTool', '\theta wrt Beam Error');

subplot4 = subplot(2,2,4,'Parent',figure1,'FontSize',24);
box(subplot4,'on');
hold(subplot4,'all');
plot_error(subplot4, reco_theta, calc_theta, 'Reconstructed \theta', '(Reco-Calc)/Calc', '\theta wrt Beam Error');

end



function [] = plot_comparison(subplot, reco, truth, x_title, y_title, plot_title)
%% Plot Comparison

plot(reco,truth,'Parent',subplot,...
    'Marker','*',...
    'LineStyle','none');

axis('square');
xlim([0 25]);
ylim([0 25]);

title(plot_title,'FontWeight','bold','FontSize',16);
xlabel(x_title,'FontWeight','bold','FontSize',16);
ylabel(y_title,'FontWeight','bold','FontSize',16);

hold on;

x_line = linspace(0,25,1000);
y_line = x_line;
plot(x_line,y_line,'-r','LineWidth',2);

hold off;

end

function [] = plot_error(subplot, reco, truth, x_title, y_title, plot_title)
%% Plot Error

difference = (reco-truth)./truth;

plot(reco,difference,'Parent',subplot,...
    'Marker','*',...
    'LineStyle','none');

axis('square');
xlim([0 25]);
ylim([-1 1]);

title(plot_title,'FontWeight','bold','FontSize',16);
xlabel(x_title,'FontWeight','bold','FontSize',16);
ylabel(y_title,'FontWeight','bold','FontSize',16);

hold on;

x_line = linspace(0,25,1000);
y_line = linspace(0,0,1000);
plot(x_line,y_line,'-r','LineWidth',2);

hold off;

end



