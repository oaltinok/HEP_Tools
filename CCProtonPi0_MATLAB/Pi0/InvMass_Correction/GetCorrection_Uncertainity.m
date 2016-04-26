function [] = GetCorrection_Uncertainity()

data_file = load('data_invMass.txt');
mc_file = load('mc_invMass.txt');

data_x = data_file(:,1);
data_y = data_file(:,2);

mc_x = mc_file(:,1);
mc_y = mc_file(:,2);

[fit_x, fit_y, NormChiSq, max_y, curvefit, fit_output] = GetFitParam(data_x,data_y);
plot_variance(curvefit);

[fit_x, fit_y, NormChiSq, max_y, curvefit, fit_output] = GetFitParam(mc_x,mc_y);
plot_variance(curvefit);


end

function [] = plot_variance(curvefit)

figure1 = figure;
axes1 = axes('Parent',figure1,...
    'FontSize',24,...
    'FontWeight','bold');

xlim([0,500]);

hold on;

% Plot Default Curve
Zero_Correction = zeros(6);
[x, y, max_y] = GetCurve_1Sigma(curvefit, 0);
maxY_text = sprintf('Default Peak = %3.2f',max_y);
p1 = plot (x,y,'-',...
    'LineWidth', 2,...
    'color','k',...
    'DisplayName', maxY_text);

[x, y, max_y] = GetCurve_1Sigma(curvefit,1);
maxY_text = sprintf('-1\\sigma Peak = %3.2f',max_y);
p2 = plot (x,y,'--',...
    'LineWidth', 2,...
    'color','b',...
    'DisplayName', maxY_text);

[x, y, max_y] = GetCurve_1Sigma(curvefit,2);
maxY_text = sprintf('+1\\sigma Peak = %3.2f',max_y);
p3 = plot (x,y,'--',...
    'LineWidth', 2,...
    'color','r',...
    'DisplayName', maxY_text);

% Plot Pi0 Inv Mass Line
pi0_line_x = linspace(134.98,134.98, 1000);
pi0_line_y = linspace(min(y),max(y)*1.1, 1000);
p4 = plot(pi0_line_x, pi0_line_y,'b-',...
    'LineWidth',2,...
    'DisplayName', '\pi^{0} Invariant Mass');

xlabel('Invariant Mass [MeV]','FontSize',24,'FontWeight','bold');
ylabel('N(Entries)','FontSize',24,'FontWeight','bold');

legend('show');

end

function [x, y, max_y] = GetCurve_1Sigma(curvefit, ind)

if ind == 0
    fit_coeffs = coeffvalues(curvefit);
    
    a1 = fit_coeffs(1);
    b1 = fit_coeffs(2);
    c1 = fit_coeffs(3);
    a2 = fit_coeffs(4);
    b2 = fit_coeffs(5);
    c2 = fit_coeffs(6);
    
else
    one_sigma_var = confint(curvefit,0.68)
    
    a1 = one_sigma_var(ind,1);
    b1 = one_sigma_var(ind,2);
    c1 = one_sigma_var(ind,3);
    a2 = one_sigma_var(ind,4);
    b2 = one_sigma_var(ind,5);
    c2 = one_sigma_var(ind,6);
    
end

x = linspace(0,275,10000);
y = a1.*exp(-((x-b1)./c1).^2) + a2.*exp(-((x-b2)./c2).^2);
% Find x Value maximizes the function
max_y = x(y == max(y));

end



