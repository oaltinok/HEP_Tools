function [] = dEdX_Profile( )
%% Actual Calculations

profiles = FormProfilesCellArray();

normalized = NormalizeAllProfiles(profiles);

% Test1(profiles,normalized,5);

% PlotAllProfiles(normalized);

PlotSurface(normalized);


end

function [] = PlotSurface(all_profiles)

% % Create figure
figure1 = figure('Position', [100, 100, 1049, 895]);

% Create axes
axes1 = axes('Parent',figure1,'FontSize',24);
% view(axes1,[-37.5 30]);
grid(axes1,'on');
hold(axes1,'all');

% Create Title
title('dEdX Profile for Normalized Shower Energy','FontWeight','bold','FontSize',24);
% Create ylabel
xlabel('Planes','FontWeight','bold','FontSize',24);
% Create xlabel
ylabel('Shower Length (nPlanes)','FontWeight','bold','FontSize',24);
% Create ylabel
zlabel('Normalized Plane Energy','FontWeight','bold','FontSize',24);

Z = GetSurfaceMatrix(all_profiles);

surf(Z','Parent',axes1);

xlim([1 60]);
ylim([1 60]);

% Create colorbar
colorbar('peer',axes1,'FontSize',24);

end

function [Z] = GetSurfaceMatrix(all_profiles)
%% SurfaceMatrix = (nPlanes,nGroups)

Z = zeros(100,100);

for shower_size = 1:100
    % Get Group
    profile_group = all_profiles{shower_size};
    % Get Number of Profiles in the Group
    nProfiles = length(profile_group);
    
    % If There are Profiles in the group: 
    if nProfiles > 0
        % Find the Most Energetic Profile
        best_profile = GetMostEnergeticProfile(profile_group);
    
        % Save best profile's normalized Energy to Z
        for plane = 1:length(best_profile) 
            Z(plane,shower_size) = best_profile(plane);
        end
    end
end

end

function [best_profile] = GetMostEnergeticProfile(profile_group)
%% Finds the Most Energetic Profile in a given Profile Group

% Empty Array to Hold Maximum Values
max_array = zeros(1, length(profile_group));

% Loop over all the Profiles in the Group and save their Maximum Value
for profile = 1:length(profile_group)
    profile_max = max(profile_group{profile});
    max_array(profile) = profile_max; 
end

% Get Maximum Indice and return the profile with that Indice
[M,I] = max (max_array);
best_profile = profile_group{I};

end


function [] = PlotAllProfiles(all_profiles)
%% Plots All Profiles

for ii = 1:length(all_profiles)
    nProfiles = length(all_profiles{ii});
    if nProfiles > 0
        disp(ii)
        PlotProfileSet(all_profiles{ii},ii);
    end
end

end

function [] = PlotProfileSet(pro_set,nPlanes)
%% Plot a Profile Set
% Profile Set - Every member has same number of planes

% Create figure
figure1 = figure('Position', [100, 100, 1049, 895],'Visible', 'off');

% Create axes
axes1 = axes('Parent',figure1,'FontSize',24);
box(axes1,'on');
hold(axes1,'all');

% Create Title
title_text = sprintf('dEdX Profile for %d Planes',nPlanes);
title(title_text,'FontWeight','bold','FontSize',24);
% Create xlabel
xlabel('N(Planes)','FontWeight','bold','FontSize',24);
% Create ylabel
ylabel('Normalized Plane Energy','FontWeight','bold','FontSize',24);

% Color Map
cc=hsv();

xlim([1 nPlanes]);

% Plot
for ii = 1:length(pro_set);
    plot(pro_set{ii},'-', ...
        'LineWidth',2,...
        'Parent',axes1,...
        'color',cc(ii,:));
    hold on;
end


% % Save Plot
fname = sprintf('/Users/oaltinok/Desktop/dEdXProfile/%d_Planes_Profile.png', nPlanes);
saveas(figure1,fname);

end


function [] = Test1(profiles,normalized,testID)
%% Test1 outputs a selected profile 

p = profiles{testID};
pn = normalized{testID};

max_value = GetMaxValue(p);

fprintf('\n\n--- Test Begin ---\n\n');
n = length(p);
disp('Actual Values');
for ii = 1:n;
    disp(p{ii})
end
fprintf('\nMax Value = %f\n\n',max_value);
disp('Normalized Values');
for ii = 1:n;
    disp(pn{ii})
end
fprintf('\n\n--- Test End ---\n\n');

end

function [profiles] = NormalizeAllProfiles(profiles)
%% Normalize Profile to 1 -- Each Profile Set Normalized uniquely

nProfiles = length(profiles);

% Loop Over all Profiles - Each Profile has its own set of profiles
for ii = 1:nProfiles
    max_value = GetMaxValue(profiles{ii});
    profiles{ii} = NormalizeProfileSet(profiles{ii}, max_value);
end

end


function [pro_set] = NormalizeProfileSet(pro_set,max_value)
%% Each Profile in the Set is normalized to the MAX value of the SET

n = length(pro_set);

% Loop Over all members of that profile set
for ii = 1:n
    pro_set{ii} = NormalizeSingleProfile(pro_set{ii},max_value);
end

end


function [single_profile] = NormalizeSingleProfile(single_profile,max_value)
%% Normalize Each Energy Value to the MAX Energy in their Profile Set

n = length(single_profile);

% Loop Over all energy values in that profile
for ii = 1:n
    single_profile(ii) = single_profile(ii) ./ max_value;
end

end


function [actual_max] = GetMaxValue(pro_set)
%% Find Maximum Value in a Set of Profiles

max_values = [];        % Empty Array to hold max values in the profile set
n = length(pro_set);    % Number of elements in that Profile Set

% Loop Over all members of that profile set
for ii = 1:n
    temp_max = max(pro_set{ii});
    max_values = [max_values, temp_max];    % Save all max for each profile member
end

actual_max = max(max_values);               % Find Max from all Max Values

end


function [profiles] = FormProfilesCellArray()
%% Forms profiles array from input file

% Pre-Allocate Memory for the structure
nMin = 1;
nMax = 100;
profiles = cell(nMin,nMax);

% Read Input File get All Profiles
input_profiles = ReadFile();
nInput = length(input_profiles);

% According to nPlanes fill profiles Cell Array
for ii = 1:nInput
    nPlanes = length(input_profiles{ii});
    tempCell{1} = input_profiles{ii};
    profiles{nPlanes} = [profiles{nPlanes}, tempCell];
end

end


function [ all_profiles ] = ReadFile()
%% Reads File and returns a cell array

% Define File and an empty Cell Array
file_source = '../Source/dEdX_Profile.txt';
fid = fopen(file_source);
all_profiles = {};

% Read Line by Line until end of line
tline = fgetl(fid);
while ischar(tline)
    x = str2num(tline);
    all_profiles = [all_profiles, x]; % Append all Lines to the Cell Array
    tline = fgetl(fid);
end

fclose(fid);

end

