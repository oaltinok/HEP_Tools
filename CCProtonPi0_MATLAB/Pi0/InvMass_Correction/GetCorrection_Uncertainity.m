function [data_CovB,mc_CovB] = GetCorrection_Uncertainity()

data_file = load('data_invMass.txt');
mc_file = load('mc_invMass.txt');

data_x = data_file(:,1);
data_y = data_file(:,2);

mc_x = mc_file(:,1);
mc_y = mc_file(:,2);

[fit_x, fit_y, NormChiSq, max_y, data_curvefit, fit_output, data_CovB] = GetFitParam(data_x,data_y);
[fit_x, fit_y, NormChiSq, max_y, mc_curvefit, fit_output, mc_CovB] = GetFitParam(mc_x,mc_y);


% plot_Exact_var(data_curvefit,data_CovB);
% plot_Exact_var(mc_curvefit,mc_CovB);

plot_1sigma_var(data_curvefit);
plot_1sigma_var(mc_curvefit);

end

function [] = plot_1sigma_var(curvefit)

figure1 = figure;
axes1 = axes('Parent',figure1,...
    'FontSize',24,...
    'FontWeight','bold');

xlim([0,500]);

hold on;

% Plot Default Curve
Zero_Correction = zeros(6);
[x, y, max_y] = GetCurve_1Sigma(curvefit, 0);
maxY_text = sprintf('Default Peak = %3.1f',max_y);
p1 = plot (x,y,'-',...
    'LineWidth', 2,...
    'color','k',...
    'DisplayName', maxY_text);

[x, y, max_y] = GetCurve_1Sigma(curvefit,1);
maxY_text = sprintf('-1\\sigma Peak = %3.1f',max_y);
p2 = plot (x,y,'--',...
    'LineWidth', 2,...
    'color','b',...
    'DisplayName', maxY_text);

[x, y, max_y] = GetCurve_1Sigma(curvefit,2);
maxY_text = sprintf('+1\\sigma Peak = %3.1f',max_y);
p3 = plot (x,y,'--',...
    'LineWidth', 2,...
    'color','r',...
    'DisplayName', maxY_text);

% Plot Pi0 Inv Mass Line
yl = ylim;
pi0_line_x = linspace(134.98,134.98, 1000);
pi0_line_y = linspace(min(yl),max(yl)*0.9, 1000);
p4 = plot(pi0_line_x, pi0_line_y,'b-',...
    'LineWidth',2,...
    'DisplayName', '\pi^{0} Invariant Mass');

xlabel('Invariant Mass [MeV]','FontSize',24,'FontWeight','bold');
ylabel('N(Entries)','FontSize',24,'FontWeight','bold');

legend('show');

end

function [] = plot_Exact_var(curvefit, CovB)

figure1 = figure;
axes1 = axes('Parent',figure1,...
    'FontSize',24,...
    'FontWeight','bold');

xlim([0,500]);

hold on;

% Plot Default Curve
[x, y, max_y] = GetCurve_Exact(curvefit, 0, CovB);
maxY_text = sprintf('Default Peak = %3.1f',max_y);
plot (x,y,'-',...
    'LineWidth', 2,...
    'color','k',...
    'DisplayName', maxY_text);


C = {'b','r','g','m','y','c'};

for ii = 1:6
    [x, y, max_y] = GetCurve_Exact(curvefit,ii, CovB);
    maxY_text = sprintf('\\lambda_{%d} Peak = %3.1f',ii,max_y);
    plot (x,y,'--',...
        'LineWidth', 2,...
        'color',C{ii},...
        'DisplayName', maxY_text);
end

% Plot Pi0 Inv Mass Line
yl = ylim;
pi0_line_x = linspace(134.98,134.98, 1000);
pi0_line_y = linspace(min(yl),max(yl)*0.8, 1000);
plot(pi0_line_x, pi0_line_y,'b-',...
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
    one_sigma_var = confint(curvefit,0.682689492137086);
    
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

function [x, y, max_y] = GetCurve_Exact(curvefit, ind, CovB)

fit_coeffs = coeffvalues(curvefit);

if ind == 0
    a1 = fit_coeffs(1);
    b1 = fit_coeffs(2);
    c1 = fit_coeffs(3);
    a2 = fit_coeffs(4);
    b2 = fit_coeffs(5);
    c2 = fit_coeffs(6);
else
    [V,D] = eig(CovB);
    var = V(:,ind).*D(ind,ind);
    
    a1 = fit_coeffs(1) + var(1);
    b1 = fit_coeffs(2) + var(2);
    c1 = fit_coeffs(3) + var(3);
    a2 = fit_coeffs(4) + var(4);
    b2 = fit_coeffs(5) + var(5);
    c2 = fit_coeffs(6) + var(6);
end

x = linspace(0,275,10000);
y = a1.*exp(-((x-b1)./c1).^2) + a2.*exp(-((x-b2)./c2).^2);
% Find x Value maximizes the function
max_y = x(y == max(y));

end






