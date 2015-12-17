function [count_array,evis_array,mean_array,min_array,max_array,weights] = GetValues()

count_array = [389,438,193,81,45];
evis_array = [25,75,125,175,250];
mean_array = [3.341,3.137,3.254,3.223,3.298];
min_array = [2,2.5,2.5,2.5,2.5];
max_array = [4.5,4,4,4,4];
rms_array = max_array - min_array;
weights = 3 - rms_array;

% Calculate Weighted Mean
weighted_mean_array = mean_array .* weights;
weighted_mean = sum(weighted_mean_array)/sum(weights);
disp(weighted_mean);

% errorbar(evis_array,mean_array,(mean_array-min_array),(max_array-mean_array),'kx');


end
