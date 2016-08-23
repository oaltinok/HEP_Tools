function [ n_carbons ] = N_CarbonAtoms()
%% Calculates Number of Carbon Atoms in Target

N_avogadro = 6.022*10^23;   % 1/mol
density = 1.04;             % g/cm3
atomic_wgt = 12.01;         % g/mol

% Number of Carbon Atoms in 1 cm3
n_cm3 = density * N_avogadro / atomic_wgt;

% Number of Carbon Atoms in Target
target_V = 0.1; % cm3
n_carbons = n_cm3 * target_V;

end