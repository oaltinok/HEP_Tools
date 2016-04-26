function [fit_x, fit_y, NormChiSq, max_y, curvefit, output] = GetFitParam(x,y)
%% Function to Calculate Fit Parameters
% Fit following range to a Double Gaussian and 
fit_range_min = 1;
fit_range_max = 28;

[curvefit, gof, output] = fit(x(fit_range_min:fit_range_max),y(fit_range_min:fit_range_max),'gauss2');
fit_coeffs = coeffvalues(curvefit);


% Get Function Coeffs
a1 = fit_coeffs(1);
b1 = fit_coeffs(2);
c1 = fit_coeffs(3);
a2 = fit_coeffs(4);
b2 = fit_coeffs(5);
c2 = fit_coeffs(6);

% Get ChiSq
fit_y = a1.*exp(-((x-b1)./c1).^2) + a2.*exp(-((x-b2)./c2).^2);
ChiSq = CalcChiSq(y(7:21),fit_y(7:21));
NormChiSq = ChiSq/(21-7+1);

% Get Fit Curve smoother (recalculate fit_y)
fit_x = linspace(0,275,1000);
fit_y = a1.*exp(-((fit_x-b1)./c1).^2) + a2.*exp(-((fit_x-b2)./c2).^2);

% Find x Value maximizes the function
max_y = fit_x(fit_y == max(fit_y));



end
