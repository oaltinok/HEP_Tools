function [  ] = GetPlots( var_name )

global pot_ratio;
pot_ratio = 0.150295;

ReadTables(var_name);

truth = MergeColumns(1) .* pot_ratio;
cv = MergeColumns(2);
difference = cv - truth;
percent_difference = (difference ./ truth) .* 100;
err_stat = MergeColumns(3);
err_syst = MergeColumns(4);


PlotVal(var_name, 'Percent (Reco-Truth)/Truth', percent_difference, strcat('Plots/',var_name,'_diff.png'),1);
PlotVal(var_name, 'Fractional Stat Error', err_stat, strcat('Plots/',var_name,'_err_stat.png'),0);
PlotVal(var_name, 'Fractional Syst Error', err_syst, strcat('Plots/',var_name,'_err_syst.png'),0);

end

function [] = PlotVal(var_name, y_label, YMatrix1, fig_name, plot_zero_line)

figure1 = figure('visible','off');

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Axis Limits
matrix_size = size(YMatrix1);
xlim([1, matrix_size(1)]);

% Create multiple lines using matrix input to plot
plot1 = plot(YMatrix1,'Parent',axes1, 'LineWidth', 2);
set(plot1(1),'DisplayName','Iteration 1');
set(plot1(2),'DisplayName','Iteration 2');
set(plot1(3),'DisplayName','Iteration 3');
set(plot1(4),'DisplayName','Iteration 4');
set(plot1(5),'DisplayName','Iteration 5');
set(plot1(6),'DisplayName','Iteration 6');
set(plot1(7),'DisplayName','Iteration 7');
set(plot1(8),'DisplayName','Iteration 8');
set(plot1(9),'DisplayName','Iteration 9');
set(plot1(10),'DisplayName','Iteration 10');

if plot_zero_line
    zero_line_x = linspace(1,matrix_size(1),1000);
    zero_line_y = linspace(0,0,1000);
    
    plot2 = plot(zero_line_x, zero_line_y, '--k', 'LineWidth', 2);
    set(get(get(plot2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    ylim([-50,50]);
else
    ylim([0,0.50]);
end


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

var = [ table_1(:,col_ind), table_2(:,col_ind), table_3(:,col_ind),...
        table_4(:,col_ind), table_5(:,col_ind), table_6(:,col_ind),...
        table_7(:,col_ind), table_8(:,col_ind), table_9(:,col_ind),...
        table_10(:,col_ind)];

end

