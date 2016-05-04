function [] = Get3DPlots()

%% Read Data
data = load('WSq.txt');

WSq = data(:,1);
Enu = data(:,2);
QSq = data(:,3);

muon_E = data(:,4);
muon_P = data(:,5);

muon_theta = data(:,6);
muon_theta_truth = data(:,7);

%% Get Values for Negative WSq
neg_Enu = GetNegativeWSq(Enu, WSq);
neg_QSq = GetNegativeWSq(QSq, WSq);
neg_muon_E = GetNegativeWSq(muon_E, WSq);
neg_muon_P = GetNegativeWSq(muon_P, WSq);
neg_muon_theta = GetNegativeWSq(muon_theta, WSq);
neg_muon_theta_truth = GetNegativeWSq(muon_theta_truth, WSq);

%% Get Values for Positive WSq
pos_Enu = GetPositiveWSq(Enu, WSq);
pos_QSq = GetPositiveWSq(QSq, WSq);
pos_muon_E = GetPositiveWSq(muon_E, WSq);
pos_muon_P = GetPositiveWSq(muon_P, WSq);
pos_muon_theta = GetPositiveWSq(muon_theta, WSq);
pos_muon_theta_truth = GetPositiveWSq(muon_theta_truth, WSq);


% createfigure(pos_Enu, pos_muon_E, neg_Enu, neg_muon_E, 'Neutrino Energy [GeV]', 'Muon Energy [GeV]');
% createfigure(pos_Enu, pos_QSq, neg_Enu, neg_QSq, 'Neutrino Energy [GeV]', 'Q^{2} [GeV^{2}]');
% createfigure(pos_muon_E, pos_muon_P, neg_muon_E, neg_muon_P, 'Muon Energy [GeV]', 'Muon Momentum [GeV]');
createfigure(pos_muon_E, pos_muon_theta, neg_muon_E, neg_muon_theta, 'Muon Energy [GeV]', 'Muon Theta wrt Beam [degree]');
% createfigure(pos_muon_theta, pos_muon_theta_truth, neg_muon_theta, neg_muon_theta_truth, 'Muon Reconstructed Theta [degree]', 'Muon Truth Theta [degree]');

%% Compare QSq Parameters
% Muon Energy vs Muon Momentum
% plot(pos_muon_E, pos_muon_P , 'k.');
% hold on;
% plot(neg_muon_E, neg_muon_P , 'MarkerFaceColor',[1 0 0],'Marker','o','LineStyle','none',...
%     'Color',[1 0 0]);
% xlabel('Muon Energy [GeV]','FontSize',24);
% ylabel('Muon Momentum [GeV]','FontSize',24);

% Muon Energy vs Muon Angle
% plot(pos_muon_E, pos_muon_theta , 'k.');
% hold on;
% plot(neg_muon_E, neg_muon_theta , 'MarkerFaceColor',[1 0 0],'Marker','o','LineStyle','none',...
%     'Color',[1 0 0]);
% xlabel('Muon Energy [GeV]','FontSize',24);
% ylabel('Muon Theta [degree]','FontSize',24);

% Muon Angle vs Truth
% plot(pos_muon_theta, pos_muon_theta_truth , 'k.');
% hold on;
% plot(neg_muon_theta, neg_muon_theta_truth , 'MarkerFaceColor',[1 0 0],'Marker','o','LineStyle','none',...
%     'Color',[1 0 0]);
% xlabel('Muon Reconstructed Theta [degree]','FontSize',24);
% ylabel('Muon Truth Theta [degree]','FontSize',24);
% ylim([0,180]);
% xlim([0,180]);

end

function [neg] = GetNegativeWSq(var,WSq)
neg = var(WSq<0);
end

function [pos] = GetPositiveWSq(var,WSq)
pos = var(WSq>=0);
end



function createfigure(pos_x,pos_y, neg_x, neg_y, label_x, label_y )
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
axis('square');
hold(axes1,'on');

% Create plot
plot(pos_x,pos_y,'DisplayName','W^{2} \geq 0','Marker','.','LineStyle','none',...
    'Color',[0 0 0]);

% Create plot
plot(neg_x,neg_y,'DisplayName','W^{2} < 0','MarkerFaceColor',[1 0 0],'Marker','o',...
    'LineStyle','none',...
    'Color',[1 0 0]);

% Create xlabel
xlabel(label_x,'FontWeight','bold');

% Create ylabel
ylabel(label_y,'FontWeight','bold');

% xlim([0 40])
ylim([0 40])

box(axes1,'on');

% Set the remaining axes properties
set(axes1,'FontSize',24,'FontWeight','bold');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','northwest');
end

