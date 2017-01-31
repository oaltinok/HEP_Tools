Data_W_All = load('Files/Data_W_All.txt');
MC_W_All = load('Files/MC_W_All.txt');
MC_W_DeltaRES = load('Files/MC_W_DeltaRES.txt');
MC_W_OtherRES = load('Files/MC_W_OtherRES.txt');
MC_W_NonRES = load('Files/MC_W_NonRES.txt');

data_x_All = Data_W_All(:,1);
data_y_All = Data_W_All(:,2);

mc_x_All = MC_W_All(:,1);
mc_y_All = MC_W_All(:,2);

mc_x_DeltaRES = MC_W_DeltaRES(:,1);
mc_y_DeltaRES = MC_W_DeltaRES(:,2);

mc_x_NonRES = MC_W_NonRES(:,1);
mc_y_NonRES = MC_W_NonRES(:,2);

mc_x_OtherRES = MC_W_OtherRES(:,1);
mc_y_OtherRES = MC_W_OtherRES(:,2);
