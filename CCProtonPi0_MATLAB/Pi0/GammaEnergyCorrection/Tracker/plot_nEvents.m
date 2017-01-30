function [ ] = plot_nEvents( )

nEvents = [2606,2388,925,334,132,58];


% Create figure
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',24,'FontWeight','bold');
hold on;

xlim([0.5,6.5]);

bar(nEvents,...
    'DisplayName','N(Events)',...
    'Parent',axes1);

% % Plot  Fit
% [x,y] = getFitXY();
% plot(x,y,'LineWidth',2,...
%     'Color','r','DisplayName','Exponential Fit');

% Create xlabel
xlabel('Energy Regions','FontSize',24,'FontWeight','bold');
% Create ylabel
ylabel('N(Events)','FontSize',24,'FontWeight','bold');
% Create legend
legend(axes1,'show');

end

function [x,y] = getFitXY()
a = 4850;
b = -0.5211;

x = linspace(1,6,1000);
y = a*exp(b*x);
end

