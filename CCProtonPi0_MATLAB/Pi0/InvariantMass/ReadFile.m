% Read File
mc = load('../Source/Pi0InvMass_MC.txt');
data = load('../Source/Pi0InvMass_Data.txt');

ratio = 3.24E20/2.30E21;

% Get Arrays
bins = mc(:,1);

bins = bins(2:end);

mc_1Track_Background = mc(:,3) * ratio;
mc_2Track_Background = mc(:,5) * ratio;

mc_1Track_Background = mc_1Track_Background(2:end);
mc_2Track_Background = mc_2Track_Background(2:end);





