function create_plot(reco,true,p)

figure1 = figure('Position', [100, 100, 1049, 895]);

axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],...
    'FontWeight','bold',...
    'FontSize',24);

xlim(axes1,[0 1]);
ylim(axes1,[0 1]);
box(axes1,'on');
hold(axes1,'all');

% Form x = y 
x = linspace(0,1,1000);
y = x;

% Plot x = y Line
plot (x,y,'b-','LineWidth',1.5);

% Add true vs reco Fit Line
yfit = polyval(p,x);
plot (x,yfit, 'r-','LineWidth',1.5);

% Add true vs reco Points
plot (reco,true, 'k.','MarkerSize',10);

xlabel('Reconstructed P_{\gamma} [GeV]','FontWeight','bold','FontSize',24);
ylabel('True P_{\gamma} [GeV]','FontWeight','bold','FontSize',24);

% Add Legend
legend_text = sprintf('y = %fx + %f', p(1), p(2));
legend('y = x',legend_text,'Location','northwest');

end
