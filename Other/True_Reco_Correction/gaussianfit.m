function [ f ] = gaussianfit( source_file )

[x,y] = get_values(source_file);

% Fit Range to get a meaningful Result 
%   Fitting around pi0 rest mass Only
fit_range_min = 8;
fit_range_max = 20;

f = fit(x(fit_range_min:fit_range_max),y(fit_range_min:fit_range_max),'gauss1')

create_plot(x,y,f);

end

function create_plot(x,y,f)

figure1 = figure('Position', [100, 100, 1049, 895]);

axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],...
    'FontWeight','bold',...
    'FontSize',24);

xlim(axes1,[0 max(x)]);
ylim(axes1,[0 1.1]);
box(axes1,'on');
hold(axes1,'all');

% Use x_smooth for Line Plots
x_smooth = linspace(0,max(x),1000);

% Plot Points
plot (x,y, 'k.','MarkerSize',10);

% Plot Correct Gaussian Peak
center = 134.98; % Pi0 Rest Mass [MeV]
sigma = 30; %[MeV]
y_expected = gaussmf(x_smooth,[sigma center]);
plot (x_smooth,y_expected,'b-','LineWidth',1.5);

% Plot Fit Result
fit_coeffs = coeffvalues(f);
y_fit = gaussmf(x_smooth,[fit_coeffs(3) fit_coeffs(2)]);
plot (x_smooth,y_fit,'r-','LineWidth',1.5);

% Add Labels
xlabel('Reconstructed Pi0 Invariant Mass [GeV]','FontWeight','bold','FontSize',24);
ylabel('Normalized N(Pions)','FontWeight','bold','FontSize',24);

% Add Legend
legend_points = sprintf('Data');
legend_y_expected = sprintf('mean = %3.2f, \\sigma = %3.2f',center,sigma);
legend_y_fit = sprintf('mean = %3.2f, \\sigma = %3.2f', fit_coeffs(2),fit_coeffs(3));
legend(legend_points, legend_y_expected, legend_y_fit, 'Location','northeast');

end

function [x,y] = get_values(source_file)
% Load File
All = load(source_file);

x = All(:,1);
y = All(:,2);

% Normalize to 1
x = x;
y = y/max(y);

end


