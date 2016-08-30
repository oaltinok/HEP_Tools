function [ ] = GetPlots( var_name )

ReadTables(var_name);

truth = MergeColumns(1);
cv = MergeColumns(2);
err_stat = MergeColumns(3);
err_syst = MergeColumns(4);

area_ratio = sum(cv) / sum(truth);
truth = truth .* area_ratio;

difference = cv - truth;
residual = (difference ./ truth);
residual(isnan(residual)) = 0.0; % In some cases truth is Zero leading a NaN in residual
pull = difference./err_stat;
pull(isnan(pull)) = 0;


%% Plots
PlotVal(var_name, 'Fractional Residual (Unfolded-Truth)/Truth', residual , strcat('Plots/',var_name,'_diff.png'),1,[-0.5,0.5]);
PlotVal(var_name, 'Pull', pull , strcat('Plots/',var_name,'_pull.png'),1,[-2e4, 2e4]);
PlotVal(var_name, 'Fractional Stat Error', err_stat, strcat('Plots/',var_name,'_err_stat.png'),0,[0.0,0.5]);
PlotVal(var_name, 'Fractional Syst Error', err_syst, strcat('Plots/',var_name,'_err_syst.png'),0,[0.0,0.5]);

PlotSums(var_name, 'Fractional Residual (Unfolded-Truth)/Truth', residual, strcat('Plots/',var_name,'_diff_sum.png'),[0, 3]);
PlotSums(var_name, 'Pull (Unfolded-Truth)/Stat_Err', pull, strcat('Plots/',var_name,'_pull_sum.png'),[1e4, 10e4]);
PlotSums(var_name, 'Fractional Stat Error', err_stat, strcat('Plots/',var_name,'_err_stat_sum.png'),[0, 1]);
PlotSums(var_name, 'Fractional Syst Error', err_syst, strcat('Plots/',var_name,'_err_syst_sum.png'),[1, 3]);

end

function [] = PlotSums(var_name,  y_label, YMatrix, fig_name, y_limits)

figure1 = figure('visible','off');

YMatrix = sum(abs(YMatrix));

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Axis Limits
x = 0:1:10;
xlim([0, 10]);

plot(x,YMatrix,'-b', 'LineWidth', 2);

% N(Iterations) = 4 Line
% y_limits = [min(YMatrix), max(YMatrix)*1.01];
ylim(y_limits);
x_4 = linspace(4,4,1000);
y_4 = linspace(y_limits(1), y_limits(2), 1000);
plot(x_4,y_4,'--r', 'LineWidth', 2);

% Create xlabel
xlabel('Iteration','FontWeight','bold');
ax = gca;
ax.XTick = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

% Create ylabel
y_label = ['Total ', y_label];
ylabel(y_label, 'FontWeight','bold');

% Create title
title_text = [var_name, ' ', y_label];
title(title_text,'Interpreter','none');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',16,'FontWeight','bold');

% Save Figure
saveas(figure1,fig_name);

end

function [] = PlotVal(var_name, y_label, YMatrix1, fig_name, plot_zero_line, y_limits)

figure1 = figure('visible','off');

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Axis Limits
matrix_size = size(YMatrix1);
xlim([1, matrix_size(1)]);

% Create multiple lines using matrix input to plot
plot1 = plot(YMatrix1,'Parent',axes1, 'LineWidth', 1);
set(plot1(1),'DisplayName','No Iteration', 'LineWidth', 2, 'Color', 'r');
set(plot1(2),'DisplayName','Iteration 1');
set(plot1(3),'DisplayName','Iteration 2');
set(plot1(4),'DisplayName','Iteration 3');
set(plot1(5),'DisplayName','Iteration 4');
set(plot1(6),'DisplayName','Iteration 5');
set(plot1(7),'DisplayName','Iteration 6');
set(plot1(8),'DisplayName','Iteration 7');
set(plot1(9),'DisplayName','Iteration 8');
set(plot1(10),'DisplayName','Iteration 9');
set(plot1(11),'DisplayName','Iteration 10');


if plot_zero_line
    zero_line_x = linspace(1,matrix_size(1),1000);
    zero_line_y = linspace(0,0,1000);
    
    plot2 = plot(zero_line_x, zero_line_y, '--k', 'LineWidth', 2);
    set(get(get(plot2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

ylim(y_limits);

% Create xlabel
xlabel('Bin Number','FontWeight','bold');

% Create ylabel
ylabel(y_label,'FontWeight','bold');

% Create title
title_text = [var_name, ' ', y_label];
title(title_text,'Interpreter','none');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',16,'FontWeight','bold');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','northeastoutside');

% Save Figure
saveas(figure1,fig_name);
end

function [] = ReadTables(var_name)
global table_0; 
global table_1; 
global table_2; 
global table_3;
global table_4;
global table_5;
global table_6; 
global table_7; 
global table_8;
global table_9;
global table_10;

table_0 = load(strcat('Input/Table_',var_name,'_0.txt'));
table_1 = load(strcat('Input/Table_',var_name,'_1.txt'));
table_2 = load(strcat('Input/Table_',var_name,'_2.txt'));
table_3 = load(strcat('Input/Table_',var_name,'_3.txt'));
table_4 = load(strcat('Input/Table_',var_name,'_4.txt'));
table_5 = load(strcat('Input/Table_',var_name,'_5.txt'));
table_6 = load(strcat('Input/Table_',var_name,'_6.txt'));
table_7 = load(strcat('Input/Table_',var_name,'_7.txt'));
table_8 = load(strcat('Input/Table_',var_name,'_8.txt'));
table_9 = load(strcat('Input/Table_',var_name,'_9.txt'));
table_10 = load(strcat('Input/Table_',var_name,'_10.txt'));
end

function [var] = MergeColumns(col_ind)
global table_0;
global table_1;
global table_2;
global table_3;
global table_4;
global table_5;
global table_6;
global table_7;
global table_8;
global table_9;
global table_10;

var = [ table_0(:,col_ind), table_1(:,col_ind), table_2(:,col_ind),...
        table_3(:,col_ind), table_4(:,col_ind), table_5(:,col_ind),...
        table_6(:,col_ind), table_7(:,col_ind), table_8(:,col_ind),...
        table_9(:,col_ind), table_10(:,col_ind)];

end

