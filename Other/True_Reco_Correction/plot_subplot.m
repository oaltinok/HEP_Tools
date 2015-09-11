function plot_subplot(subplot,x,y,f,title_text)

% Use x_smooth for Line Plots
x_smooth = linspace(0,max(x),1000);

% Plot Points
% plot (x,y, 'k.','MarkerSize',10,'Parent',subplot);
bar (x,y,'Parent',subplot,...
    'FaceColor',[0.0431372560560703 0.517647087574005 0.780392169952393],...
    'BarWidth',1);

ylim([0 1.2]);

% Plot Correct Gaussian Peak
center = 134.98; % Pi0 Rest Mass [MeV]
sigma = 30; %[MeV]
y_expected = gaussmf(x_smooth,[sigma center]);
plot (x_smooth,y_expected,'g-','LineWidth',1.5,'Parent',subplot);

% Plot Fit Result
fit_coeffs = coeffvalues(f);
y_fit = gaussmf(x_smooth,[fit_coeffs(3) fit_coeffs(2)]);
plot (x_smooth,y_fit,'r-','LineWidth',1.5,'Parent',subplot);

% Add Labels
title(title_text,'FontWeight','bold','FontSize',24);
xlabel('Reconstructed Pi0 Invariant Mass [MeV]','FontSize',24);
ylabel('Normalized N(Pions)','FontSize',24);

% Add Legend
legend_points = sprintf('Normalized Value');
legend_y_expected = sprintf('mean = %3.2f, \\sigma = %3.2f',center,sigma);
legend_y_fit = sprintf('mean = %3.2f, \\sigma = %3.2f', fit_coeffs(2),fit_coeffs(3));
legend(subplot,'show',legend_points, legend_y_expected, legend_y_fit, 'Location','northeast');

end
