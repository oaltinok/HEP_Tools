function [ ] = GetFit( )

[data_x, data_y] = ReadFile('W_data.txt');
[MC_x, MC_y] = ReadFile('W_MC.txt');

GeV_to_MeV = 1000;
POT_Ratio = 0.150394;

data_x = data_x * GeV_to_MeV;
MC_x = MC_x * GeV_to_MeV;
MC_y = MC_y * POT_Ratio;

CreateFigure(data_x,data_y, MC_x, MC_y);


end

function [] = CreateFigure(data_x, data_y, mc_x, mc_y)
% Create figure
% figure1 = figure('visible','off');
figure1 = figure('Position', [100, 100, 1024, 860]);

% Create axes
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontWeight','bold','FontSize',16);
box(axes1,'on');
hold(axes1,'all');

% Create plot for Data
errorbar(data_x,data_y,sqrt(data_y),...
        'LineStyle','none',...
        'Marker','o',...
        'MarkerFaceColor','k',...
        'MarkerEdgeColor','k',...
        'Color','k',...
        'DisplayName','Data');

xlim([0 3000]);
ylim([0 8500]);

hold on;

% Create plot for MC
stairs(mc_x,mc_y,...
        'MarkerEdgeColor','r','LineWidth',2,'Color','r',...
        'DisplayName','MC');
    
curve_x = linspace(0,3000,3000);
    
% Plot Total Fit Curve
[fit_x, fit_y] = Fit_Curve(curve_x);
plot(fit_x, fit_y,...
    'Color','b',...
    'LineWidth',2,...
    'DisplayName','Fit Result')

% Plot Breit Wigner Part
a1 = 6332;
gamma = 257.1;
mass = 1213;
[Breit_Wigner_x, Breit_Wigner_y] = Breit_Wigner_Curve(curve_x, a1, gamma, mass);
plot(Breit_Wigner_x, Breit_Wigner_y,...
    'Color','m',...
    'LineStyle','--',... 
    'LineWidth',2,...
    'DisplayName','Fit: Breit-Wigner')

% Plot Gauss Part
a2 = 2467;
b2 = 1578;
c2 = 488;
[Gauss_x, Gauss_y] = Gauss_Curve(curve_x, a2, b2, c2);
plot(Gauss_x, Gauss_y,...
    'Color','c',...
    'LineStyle','--',... 
    'LineWidth',2,...
    'DisplayName','Fit: Gaussian')

% Calculate ChiSquare
[fit_x, fit_y] = Fit_Curve(data_x);
[ChiSq, NormChiSq] = CalcChiSq(data_y,fit_y);
disp(ChiSq);
disp(NormChiSq);

    
% Create xlabel
xlabel('W [MeV]','FontWeight','bold','FontSize', 16);

% Create ylabel
ylabel('N(Events)','FontWeight','bold','FontSize', 16);

% Create legend
legend(axes1,'show');

print(figure1,'W_Fit.png','-dpng')

end

function [x,y] = Breit_Wigner_Curve(x, amp, gamma, mass)

y = amp * ((gamma.^2)/4) ./ ((x-mass).^2 + ((gamma.^2)/4));

end

function [x,y] = Gauss_Curve(x, a, b, c)

y = a * exp(-((x-b)/c).^2);

end

function [x,y] = Fit_Curve(x)

% Breit Wigner Parameters
a1 = 6332;
gamma = 257.1;
mass = 1213;

% Gaussian Parameters
a2 = 2467;
b2 = 1578;
c2 = 488;

[dummy_x, breit_wigner_part] = Breit_Wigner_Curve(x, a1,gamma,mass);
[dummy_x, gauss_part] =  Gauss_Curve(x, a2, b2, c2);

y =  breit_wigner_part + gauss_part;

end



