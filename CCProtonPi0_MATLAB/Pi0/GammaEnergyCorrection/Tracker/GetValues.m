function [count_array,evis_array,mean_array,min_array,max_array,weights] = GetValues()

% See Excel Sheet HWHM_Table for values
% count_array = [2606,2388,925,332,182,62];
% evis_array = [25,75,125,175,250,500];
% mean_array = [1.48,1.478,1.46,1.443,1.433,1.312];
% min_array = [1.1,1.2,1.2,1.2,1.2,1.1];
% max_array = [1.85,1.65,1.6,1.6,1.55,1.5];
% rms_array = max_array - min_array;
% weights = 1 - rms_array;

count_array = [2606,2388,925,332,182];
evis_array = [25,75,125,175,250];
mean_array = [1.48,1.478,1.46,1.443,1.433];
min_array = [1.1,1.2,1.2,1.2,1.2];
max_array = [1.85,1.65,1.6,1.6,1.55];
rms_array = max_array - min_array;
weights = 1 - rms_array;

% errorbar(evis_array,mean_array,(mean_array-min_array),(max_array-mean_array),'kx');


end
