function [  ] = GetFit_Subtraction( )

% Read File
mc = load('../Source/Pi0InvMass_MC.txt');
data = load('../Source/Pi0InvMass_Data.txt');

ratio = 3.24E20/2.30E21;

% Get Arrays
bins = mc(:,1);

mc_1Track_Signal = mc(:,2);
mc_2Track_Signal = mc(:,4);
data_1Track = data(:,2);
data_2Track = data(:,3);

% POT Normalize
mc_1Track_Signal = mc_1Track_Signal * ratio;
mc_2Track_Signal = mc_2Track_Signal * ratio;

% Subtract Background from Data
data_1Track = subtract_background(95.33, 0.0002472,0.04135,-0.02634, 2.17e-05, data_1Track, bins);
data_2Track = subtract_background(231.4, 0.0001333,0.1125,-0.02645, 2.311e-05, data_2Track, bins);

% Get Fit Results
f_mc_1Track = gaussianfit(bins,mc_1Track_Signal);
f_mc_2Track = gaussianfit(bins,mc_2Track_Signal);
f_data_1Track = gaussianfit(bins,data_1Track);
f_data_2Track = gaussianfit(bins,data_2Track);

% Create figure
figure1 = figure;

% Create & Plot subplots
subplot1 = subplot(2,2,1,'Parent',figure1,'FontSize',24);
box(subplot1,'on');
hold(subplot1,'all');
plot_subplot_mc(subplot1, bins, mc_1Track_Signal, f_mc_1Track, 'MC 1Track');

subplot2 = subplot(2,2,3,'Parent',figure1,'FontSize',24);
box(subplot2,'on');
hold(subplot2,'all');
plot_subplot_mc(subplot2, bins, mc_2Track_Signal, f_mc_2Track, 'MC 2Track');

subplot3 = subplot(2,2,2,'Parent',figure1,'FontSize',24);
box(subplot3,'on');
hold(subplot3,'all');
plot_subplot_data(subplot3, bins, data_1Track, f_data_1Track, 'Data 1Track');

subplot4 = subplot(2,2,4,'Parent',figure1,'FontSize',24);
box(subplot4,'on');
hold(subplot4,'all');
plot_subplot_data(subplot4, bins, data_2Track, f_data_2Track, 'Data 2Track');


end


function [data] = subtract_background(A,B,C,a,b,data,bins)
%% Background Function has Exponential Form
% Functional form was found fitting the background curve to y = a*exp(-bx)

bckg = (A + B*bins + C*bins.^2).*exp(a*bins + b*bins.^2);

data = data - bckg;

% Zero negative values
data(data<0) = 0;

end

function plot_subplot_mc(subplot,x,signal,f,title_text)

% Use x_smooth for Line Plots
x_smooth = linspace(0,max(x),1000);

% Plot Points
% plot (x,y, 'k.','MarkerSize',10,'Parent',subplot);
bar (x,signal,...
    'Parent',subplot,...
    'FaceColor',[0 0.498039215803146 0],...
    'BarWidth',1,...
    'DisplayName','Signal');

% Plot Correct Gaussian Peak
center = 134.98; % Pi0 Rest Mass [MeV]
sigma = 30; %[MeV]
y_expected = gaussmf(x_smooth,[sigma center]);
legend_y_expected = sprintf('Expected: mean = %3.2f, \\sigma = %3.2f',center,sigma);
plot (x_smooth,y_expected*max(signal),'g-',...
    'LineWidth',1.5,...
    'Parent',subplot,...
    'DisplayName',legend_y_expected);

% Plot Fit Result
fit_coeffs = coeffvalues(f);
y_fit = gaussmf(x_smooth,[fit_coeffs(3) fit_coeffs(2)]);
legend_y_fit = sprintf('Fitted: mean = %3.2f, \\sigma = %3.2f', fit_coeffs(2),fit_coeffs(3));
plot (x_smooth,y_fit*max(signal),'r-',...
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

function plot_subplot_data(subplot,x,y,f,title_text)

% Use x_smooth for Line Plots
x_smooth = linspace(0,max(x),1000);

% Plot Points
% plot (x,y, 'k.','MarkerSize',10,'Parent',subplot);
bar (x,y,...
    'Parent',subplot,...
    'BarWidth',1,...
    'DisplayName','Background Subtracted Data');

% Plot Correct Gaussian Peak
center = 134.98; % Pi0 Rest Mass [MeV]
sigma = 30; %[MeV]
y_expected = gaussmf(x_smooth,[sigma center]);
legend_y_expected = sprintf('Expected: mean = %3.2f, \\sigma = %3.2f',center,sigma);
plot (x_smooth,y_expected*max(y),'g-',...
    'LineWidth',1.5,...
    'Parent',subplot,...
    'DisplayName',legend_y_expected);

% Plot Fit Result
fit_coeffs = coeffvalues(f);
y_fit = gaussmf(x_smooth,[fit_coeffs(3) fit_coeffs(2)]);
legend_y_fit = sprintf('Fitted: mean = %3.2f, \\sigma = %3.2f', fit_coeffs(2),fit_coeffs(3));
plot (x_smooth,y_fit*max(y),'r-',...
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

