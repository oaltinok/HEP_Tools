function [] = BackgroundSubtraction()

% Read File
mc = load('../Source/Pi0InvMass_MC.txt');

ratio = 3.24E20/2.30E21;

% Get Arrays
bins = mc(:,1);
mc_1Track_Background = mc(:,3) * ratio;
mc_2Track_Background = mc(:,5) * ratio;

% Create figure
figure1 = figure;

% Create & Plot subplots
subplot1 = subplot(2,2,1,'Parent',figure1,'FontSize',24);
box(subplot1,'on');
hold(subplot1,'all');
plot_subplot(subplot1, bins, mc_1Track_Background, 95.33, 0.0002472,0.04135,-0.02634, 2.17e-05, 'MC 1Track');

subplot2 = subplot(2,2,3,'Parent',figure1,'FontSize',24);
box(subplot2,'on');
hold(subplot2,'all');
plot_subplot(subplot2, bins, mc_2Track_Background, 231.4, 0.0001333,0.1125,-0.02645, 2.311e-05, 'MC 2Track');

end


function [y] = fit_background(a,b)
%% Background Function has Exponential Form
% Functional form was found fitting the background curve to y = a*exp(-bx)

x = linspace(0,600,2000);
y= a*exp(-b*x);

end

function plot_subplot(subplot,x,y,A,B,C,a,b,title_text)

% Use x_smooth for Line Plots
x_smooth = linspace(0,max(x),1000);

% Plot Points
bar (x,y,...
    'Parent',subplot,...
    'BarWidth',1,...
    'FaceColor',[0.749019622802734 0 0.749019622802734],...
    'DisplayName','Background');


% Plot Fit
y_fit = (A+B*x_smooth+C*x_smooth.^2).*exp(a*x_smooth+b*x_smooth.^2);
legend_y_fit = sprintf('(A + Bx + Cx^{2})exp(ax+bx^{2})');
plot (x_smooth,y_fit,'g-',...
    'LineWidth',1.5,...
    'Parent',subplot,...
    'DisplayName',legend_y_fit);

% Add Labels
title(title_text,'FontWeight','bold','FontSize',24);
xlabel('Reconstructed Pi0 Invariant Mass [MeV]','FontSize',24);
ylabel('Normalized N(Pions)','FontSize',24);

% Add Legend
legend1 = legend(subplot,'show');
set(legend1,'FontSize',16);

end



