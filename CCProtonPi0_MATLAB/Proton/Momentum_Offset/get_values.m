% Get Values
P_All = load('FailChecks_2Prong.txt');
P_true = P_All(:,1);
P_reco = P_All(:,2);

% Calculate Error
P_err = (P_reco - P_true) ./ P_true;

% Get only the good ones (-0.1<x<0.1)
good_ind = abs(P_err) < 0.1;

% Use only good match 
P_true = P_true(good_ind);
P_reco = P_reco(good_ind);

% Fit Results --> true = a * reco + b
a = 1.097;
b = -45.43;
fit_line = a * P_reco + b;

% x=y Line
x_y_line = P_reco;

plot(P_reco,P_true, 'k.');

axis('square');
xlim([0 2000]);
ylim([0 2000]);
xlabel('Proton Momentum Reconstructed [MeV]');
ylabel('Proton Momentum True [MeV]');

hold on;

plot(P_reco,fit_line,'b-');
plot(P_reco,x_y_line,'g-');

