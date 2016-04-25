function [ data_x, data_y, mc_x, mc_y ] = GetCorrection()

data_file = load('data_invMass.txt');
mc_file = load('mc_invMass.txt');

data_x = data_file(:,1);
data_y = data_file(:,2);

mc_x = mc_file(:,1);
mc_y = mc_file(:,2);

plot_fit(data_x, data_y, 0);
plot_fit(mc_x, mc_y, 1);

data_x = data_x(1:28);
data_y = data_y(1:28);

mc_x = mc_x(1:28);
mc_y = mc_y(1:28);

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


% Calculate ChiSquare and Print
if isMC == 1
    [y_fit, max_y] = GetFitY_MC(x);
else
    [y_fit, max_y] = GetFitY_Data(x);
end
ChiSq = CalcChiSq(y(7:21),y_fit(7:21));
NormChiSq = ChiSq/(21-7+1);
ChiSq_text = sprintf('Signal Region Norm(\\chi^{2}) = %3.2f',NormChiSq);
text(max(x)*0.50,max(y)*0.85,ChiSq_text,...
    'FontSize',24,...
    'FontWeight','bold');

% Plot Fit Result smoother
x_fit = linspace(60,200,1000);
if isMC == 1
    [y_fit, max_y] = GetFitY_MC(x_fit);
else
    [y_fit, max_y] = GetFitY_Data(x_fit);
end
p4 = plot (x_fit,y_fit,'r-',...
    'LineWidth',2,...
    'DisplayName', 'Gaussian Fit');

MaxY_text = sprintf('Fit Peaks at %3.2f',max_y);
text(max(x)*0.50,max(y)*0.78,MaxY_text,...
    'FontSize',24,...
    'FontWeight','bold');


xlabel('Invariant Mass [MeV]','FontSize',24,'FontWeight','bold');
ylabel('N(Entries)','FontSize',24,'FontWeight','bold');

legend([p1 p2 p3 p4]);


end


function [fit_y,max_y] = GetFitY_Data(x)

a1 =       138.4;
b1 =         140;
c1 =       27.47;
a2 =       216.8;
b2 =       125.3;
b2 = b1;
c2 =       80.88;
fit_y =  a1.*exp(-((x-b1)./c1).^2) + a2.*exp(-((x-b2)./c2).^2);


       a1 =       323.5;
       b1 =       131.2;
       c1 =       59.58;
fit_y =  a1*exp(-((x-b1)./c1).^2);

max_y = x(fit_y == max(fit_y));

end

function [fit_y,max_y] = GetFitY_MC(x)

a1 =       815.9;
b1 =       133.1;
c1 =       24.33;
a2 =        1669;
b2 =       117.7;
b2 = b1;
c2 =       72.87;

fit_y =  a1.*exp(-((x-b1)./c1).^2) + a2.*exp(-((x-b2)./c2).^2);

       a1 =        2273;
       b1 =       123.4;
       c1 =       56.91;
fit_y =  a1*exp(-((x-b1)./c1).^2);



max_y = x(fit_y == max(fit_y));
end




