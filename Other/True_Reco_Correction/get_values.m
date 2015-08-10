function [reco, true] = get_values(source_file)
% Load File
All = load(source_file);

reco = All(:,1);
true = All(:,2);

% Calculate Error
error = (reco - true) ./ true;

% Get only the good ones (-0.1<x<0.1)
good_ind = abs(error) < 0.3;

% Use only good match
reco = reco(good_ind);
true = true(good_ind);
end
