function [  ] = CompareSCALModels( )

model_1 = load('model1.dat'); 
model_2 = load('model2.dat'); 

error = Inf;

plot_comparison(model_1,model_2,error);

disp('Old Model');
get_table(model_1,error);
disp(' ');
disp('New Model');
get_table(model_2,error);


end

function [] = get_table(model,max_error)

[trkr_t, trkr_r, scal_t, scal_r] = getCols(model);

% Replace True == 0 Values
trkr_t(trkr_t == 0) = 0.001;
scal_t(scal_t == 0) = 0.001;

trkr_error = getError(trkr_t,trkr_r);
scal_error = getError(scal_t,scal_r);

trkr = (trkr_error < max_error);
scal = (scal_error < max_error);

trkr_text = sprintf('Tracker N(Events) = %.0f, max_error = %3.2f', sum(trkr), max_error);
scal_text = sprintf('Side ECAL N(Events) = %.0f, max_error = %3.2f', sum(scal), max_error);

disp(trkr_text);
disp(scal_text);

end


function [error] = getError(true, reco)

error = abs((reco-true)./true);

end

function [c1,c2,c3,c4] = getCols(m)

c1 = m(:,1);
c2 = m(:,2);
c3 = m(:,3);
c4 = m(:,4);

end

function [] = plot_reco_true(plot_id, true,reco,axis_lim,label,error)

box(plot_id,'on');
hold(plot_id,'all');

% Plot Points
plot(reco,true,'Parent',plot_id,'Marker','*','LineStyle','none');

% Plot x = y Line
x = linspace(0,1500,1000);
y = x;
plot(x,y,'Parent',plot_id,'LineWidth',2,'Color',[1 0 0]);

xlim([0 axis_lim]);
ylim([0 axis_lim]);

% Plot x = y_err
y_err = get_y_err(x,error);
plot(x,y_err,'Parent',plot_id,'LineWidth',2,'Color',[1 0 0],'LineStyle','--');

% Plot x = -y_err
y_err = get_y_err(x,-error);
plot(x,y_err,'Parent',plot_id,'LineWidth',2,'Color',[1 0 0],'LineStyle','--');

% Create xlabel
xlabel(sprintf('%s%s',label,' Reco N(Hits)'),'FontWeight','bold','FontSize',16);

% Create ylabel
ylabel(sprintf('%s%s',label,' True N(Hits)'),'FontWeight','bold','FontSize',16);

end

function [] = plot_comparison(model_1,model_2,error)

[m1_trkr_t, m1_trkr_r, m1_scal_t, m1_scal_r] = getCols(model_1);
[m2_trkr_t, m2_trkr_r, m2_scal_t, m2_scal_r] = getCols(model_2);

% Create figure
figure1 = figure;

subplot1 = subplot(2,2,1,'Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',16,'FontWeight','bold');
plot_reco_true(subplot1, m1_trkr_t,m1_trkr_r, 1200,'Old Model Tracker',error);

subplot2 = subplot(2,2,2,'Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',16,'FontWeight','bold');
plot_reco_true(subplot2, m2_trkr_t,m2_trkr_r, 1200,'New Model Tracker',error);

subplot3 = subplot(2,2,3,'Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',16,'FontWeight','bold');
plot_reco_true(subplot3, m1_scal_t,m1_scal_r, 500,'Old Model Side ECAL',error);

subplot4 = subplot(2,2,4,'Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',16,'FontWeight','bold');
plot_reco_true(subplot4, m2_scal_t,m2_scal_r, 500,'New Model Side ECAL',error);

% Create figure
figure2 = figure;

subplot1 = subplot(2,2,1,'Parent',figure2,'PlotBoxAspectRatio',[1 1 1],'FontSize',16,'FontWeight','bold');
plot_reco_true(subplot1, m1_trkr_t,m1_trkr_r, 200,'Old Model Tracker',error);

subplot2 = subplot(2,2,2,'Parent',figure2,'PlotBoxAspectRatio',[1 1 1],'FontSize',16,'FontWeight','bold');
plot_reco_true(subplot2, m2_trkr_t,m2_trkr_r, 200,'New Model Tracker',error);

subplot3 = subplot(2,2,3,'Parent',figure2,'PlotBoxAspectRatio',[1 1 1],'FontSize',16,'FontWeight','bold');
plot_reco_true(subplot3, m1_scal_t,m1_scal_r, 100,'Old Model Side ECAL',error);

subplot4 = subplot(2,2,4,'Parent',figure2,'PlotBoxAspectRatio',[1 1 1],'FontSize',16,'FontWeight','bold');
plot_reco_true(subplot4, m2_scal_t,m2_scal_r, 100,'New Model Side ECAL',error);
end


function [y_err] = get_y_err(x,error)

y_err = x.*(1+error);

end

function [y_dev] = get_y_dev(x,deviation)

y_dev = x + deviation;

end