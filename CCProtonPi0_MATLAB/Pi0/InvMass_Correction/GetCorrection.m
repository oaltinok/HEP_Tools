function [ ] = GetCorrection()

data_file = load('data_invMass.txt');
mc_file = load('mc_invMass.txt');

data_x = data_file(:,1);
data_y = data_file(:,2);

mc_x = mc_file(:,1);
mc_y = mc_file(:,2);

plot_fit(data_x, data_y, 0);
plot_fit(mc_x, mc_y, 1);

end

function [] = plot_fit(x,y,isMC)

figure1 = figure;
axes1 = axes('Parent',figure1,...
    'PlotBoxAspectRatio',[1 1 1],...
    'FontSize',24,...
    'FontWeight','bold');

hold on;

xlim([0,500]);
ylim([0,max(y)*1.2]);


% Plot Data Points
if isMC == 1
    legend_text = 'MC Points';
else
    legend_text = 'Data Points';
end
p1 = errorbar(x,y,sqrt(y),...
    'MarkerFaceColor',[0 0 0],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 0],...
    'DisplayName', legend_text);


% Plot Pi0 Inv Mass Line
pi0_line_x = linspace(134.98,134.98, 1000);
pi0_line_y = linspace(min(y),max(y)*1.1, 1000);
p2 = plot(pi0_line_x, pi0_line_y,'b-',...
    'LineWidth',2,...
    'DisplayName', '\pi^{0} Invariant Mass');

% Plot Signal Region
signal_line_x = linspace(60,60, 1000);
signal_line_y = linspace(min(y),max(y)*1.1, 1000);
p3 = plot(signal_line_x, signal_line_y,'b--',...
    'LineWidth',2,...
    'DisplayName', 'Signal Selection');

signal_line_x = linspace(200,200, 1000);
signal_line_y = linspace(min(y),max(y)*1.1, 1000);
plot(signal_line_x, signal_line_y,'b--',...
    'LineWidth',2);


% Get Fit Parameters and Plot
[fit_x, fit_y, NormChiSq, max_y] = GetFitParam(x,y);
p4 = plot (fit_x,fit_y,'r-',...
    'LineWidth',2,...
    'DisplayName', 'Double Gaussian Fit');

% Write NormChiSq & MaxY on Plot
ChiSq_text = sprintf('Signal Region Norm(\\chi^{2}) = %3.2f',NormChiSq);
text(max(x)*0.50,max(y)*0.85,ChiSq_text,...
    'FontSize',24,...
    'FontWeight','bold');

maxY_text = sprintf('Fit Peaks at %3.2f',max_y);
text(max(x)*0.50, max(y)*0.78, maxY_text,...
    'FontSize',24,...
    'FontWeight','bold');

xlabel('Invariant Mass [MeV]','FontSize',24,'FontWeight','bold');
ylabel('N(Entries)','FontSize',24,'FontWeight','bold');

legend([p1 p2 p3 p4]);

end







