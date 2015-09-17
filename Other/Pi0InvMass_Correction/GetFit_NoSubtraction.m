function [] = GetFit_NoSubtraction()

% Read File
mc = load('../Source/Pi0InvMass_MC.txt');
data = load('../Source/Pi0InvMass_Data.txt');

ratio = 3.24E20/2.30E21;

% Get Arrays
bins = mc(:,1);

mc_1Track_Signal = mc(:,2);
mc_1Track_Background = mc(:,3);

mc_2Track_Signal = mc(:,4);
mc_2Track_Background = mc(:,5);

data_1Track = data(:,2);
data_2Track = data(:,3);

% POT Normalize
mc_1Track_Signal = mc_1Track_Signal * ratio;
mc_2Track_Signal = mc_2Track_Signal * ratio;
mc_1Track_Background = mc_1Track_Background * ratio;
mc_2Track_Background = mc_2Track_Background * ratio;
mc_1Track = mc_1Track_Signal + mc_1Track_Background;
mc_2Track = mc_2Track_Signal + mc_2Track_Background;

% Get Fit Results
f_mc_1Track = gaussianfit(bins,mc_1Track);
f_mc_2Track = gaussianfit(bins,mc_2Track);
f_data_1Track = gaussianfit(bins,data_1Track);
f_data_2Track = gaussianfit(bins,data_2Track);

% Create figure
figure1 = figure;

% Create & Plot subplots
subplot1 = subplot(2,2,1,'Parent',figure1,'FontSize',24);
box(subplot1,'on');
hold(subplot1,'all');
plot_subplot_mc(subplot1, bins, mc_1Track_Signal, mc_1Track_Background, f_mc_1Track, 'MC 1Track');

subplot2 = subplot(2,2,3,'Parent',figure1,'FontSize',24);
box(subplot2,'on');
hold(subplot2,'all');
plot_subplot_mc(subplot2, bins, mc_2Track_Signal, mc_2Track_Background, f_mc_2Track, 'MC 2Track');

subplot3 = subplot(2,2,2,'Parent',figure1,'FontSize',24);
box(subplot3,'on');
hold(subplot3,'all');
plot_subplot_data(subplot3, bins, data_1Track, f_data_1Track, 'Data 1Track');

subplot4 = subplot(2,2,4,'Parent',figure1,'FontSize',24);
box(subplot4,'on');
hold(subplot4,'all');
plot_subplot_data(subplot4, bins, data_2Track, f_data_2Track, 'Data 2Track');

end

function plot_subplot_data(subplot,x,y,f,title_text)

% Use x_smooth for Line Plots
x_smooth = linspace(0,max(x),1000);

% Plot Points
% plot (x,y, 'k.','MarkerSize',10,'Parent',subplot);
bar (x,y,...
    'Parent',subplot,...
    'BarWidth',1,...
    'DisplayName','All Data Points');

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

function plot_subplot_mc(subplot,x,signal,background,f,title_text)

% Use x_smooth for Line Plots
x_smooth = linspace(0,max(x),1000);

% Plot Points
% plot (x,y, 'k.','MarkerSize',10,'Parent',subplot);
bar1 = bar (x,[background signal],'stack',...
    'Parent',subplot,...
    'BarWidth',1);
set(bar1(1),'FaceColor',[0.749019622802734 0 0.749019622802734],...
    'DisplayName','Background');
set(bar1(2),'FaceColor',[0 0.498039215803146 0],...
    'DisplayName','Signal');

% Plot Correct Gaussian Peak
center = 134.98; % Pi0 Rest Mass [MeV]
sigma = 30; %[MeV]
y_expected = gaussmf(x_smooth,[sigma center]);
legend_y_expected = sprintf('Expected: mean = %3.2f, \\sigma = %3.2f',center,sigma);
plot (x_smooth,y_expected*max(background+signal),'g-',...
    'LineWidth',1.5,...
    'Parent',subplot,...
    'DisplayName',legend_y_expected);

% Plot Fit Result
fit_coeffs = coeffvalues(f);
y_fit = gaussmf(x_smooth,[fit_coeffs(3) fit_coeffs(2)]);
legend_y_fit = sprintf('Fitted: mean = %3.2f, \\sigma = %3.2f', fit_coeffs(2),fit_coeffs(3));
plot (x_smooth,y_fit*max(background+signal),'r-',...
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

