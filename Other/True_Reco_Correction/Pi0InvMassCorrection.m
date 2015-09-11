function [  ] = Pi0InvMassCorrection(  )

% Read File
mc = load('../Source/Pi0InvMass_MC.txt');
data = load('../Source/Pi0InvMass_Data.txt');

% Get Arrays
bins = mc(:,1);
mc_1Track = mc(:,2);
mc_2Track = mc(:,3);
data_1Track = data(:,2);
data_2Track = data(:,3);

% Normalize to 1
mc_1Track = mc_1Track / max(mc_1Track);
mc_2Track = mc_2Track / max(mc_2Track);
data_1Track = data_1Track / max(data_1Track);
data_2Track = data_2Track / max(data_2Track);

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
plot_subplot(subplot1, bins, mc_1Track, f_mc_1Track, 'MC 1Track');

subplot2 = subplot(2,2,3,'Parent',figure1,'FontSize',24);
box(subplot2,'on');
hold(subplot2,'all');
plot_subplot(subplot2, bins, mc_2Track, f_mc_2Track, 'MC 2Track');

subplot3 = subplot(2,2,2,'Parent',figure1,'FontSize',24);
box(subplot3,'on');
hold(subplot3,'all');
plot_subplot(subplot3, bins, data_1Track, f_data_1Track, 'Data 1Track');

subplot4 = subplot(2,2,4,'Parent',figure1,'FontSize',24);
box(subplot4,'on');
hold(subplot4,'all');
plot_subplot(subplot4, bins, data_2Track, f_data_2Track, 'Data 2Track');



end

