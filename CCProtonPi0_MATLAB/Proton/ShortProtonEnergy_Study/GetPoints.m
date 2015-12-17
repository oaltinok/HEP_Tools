function [ evis,avg,n_events ] = GetPoints(  )

% Read Data
data = load('CCProtonPi0_protons.txt');

evis = data(:,1);
avg = data(:,2);
n_events = data(:,3);

evis = evis(evis<200);
avg = avg(evis<200);
n_events = n_events(evis<200);

end

