function [ p ] = linefit( source_file )

[x, y] = get_values(source_file);

p = polyfit(x,y,1);

create_plot(x,y,p);

end

function create_plot(reco,true,p)

figure1 = figure('Position', [100, 100, 1049, 895]);

axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],...
    'FontWeight','bold',...
    'FontSize',24);

xlim(axes1,[0 1]);
ylim(axes1,[0 1]);
box(axes1,'on');
hold(axes1,'all');

% Add true vs reco Points
plot (reco,true, 'k.','MarkerSize',10);

% Form y = x 
x = linspace(0,1,1000);
y = x;
% Plot y = x Line
plot (x,y,'b-','LineWidth',1.5);

% Add true vs reco Fit Line
yfit = polyval(p,x);
plot (x,yfit, 'r-','LineWidth',1.5);

% Form y = 1.2x -- Best Fit for Gamma Energy
y_best = 1.2*x;
% Plot y = x Line
plot (x,y_best,'g-','LineWidth',1.5);


xlabel('Reconstructed P_{\gamma} [GeV]','FontWeight','bold','FontSize',24);
ylabel('True P_{\gamma} [GeV]','FontWeight','bold','FontSize',24);

% Add Legend'
legend_points = sprintf('abs(y-x/y) < %2.1f',0.3);
legend_y_x = sprintf('y = x');
legend_fit = sprintf('y_{fit} = %4.3fx + %4.3f', p(1), p(2));
legend_best = sprintf('y_{best} = %2.1fx',1.2);
legend(legend_points, legend_y_x, legend_fit, legend_best, 'Location','northwest');

end

function [reco, true] = get_values(source_file)
% Load File
All = load(source_file);

reco = All(:,1);
true = All(:,2);

% Calculate Error
error = (reco - true) ./ true;

% Get only the good ones (-0.1<x<0.1)
good_ind = abs(error) < 0.3;

% Use only good match
reco = reco(good_ind);
true = true(good_ind);
end


