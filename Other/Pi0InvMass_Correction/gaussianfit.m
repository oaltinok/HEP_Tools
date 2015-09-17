function [ f ] = gaussianfit(x,y)

% Fit Range to get a meaningful Result 
%   Fitting around pi0 rest mass (135 MeV)
%   Range: -60, +60 MeV
fit_range_min = 8;
fit_range_max = 20;

f = fit(x(fit_range_min:fit_range_max),y(fit_range_min:fit_range_max),'gauss1')

end




