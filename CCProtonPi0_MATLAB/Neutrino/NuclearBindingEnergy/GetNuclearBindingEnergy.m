function [ Eb_nucleon ] = GetNuclearBindingEnergy( n_protons, n_neutrons, atomic_mass )
%% Calculates Nuclear Binding Energy
%   Input Arguments: (n_protons, n_neutrons, atomic_mass)
%       Note: Use Carbon Atomic Mass = 12.0u 
%   Returns Nuclear Binding Energy per Nucleon in MeV
%   Uses Equation from
%   K. Krane. (1996). Modern Physics Second Edition(pp.380-381), USA, John Wiley & Sons, Inc. 



% Constants
neutron_mass = 1.008665;        % [u]
hydrogen_mass = 1.007825;       % [u]
c_sq = 931.5;                   % [MeV/u] -- Convertes Mass Units to Energy

% Nuclear Binding Energy Eq. 12.4
Eb = (n_neutrons.*neutron_mass + n_protons.*hydrogen_mass - atomic_mass) .* c_sq;

% Binding Energy per Nucleon
Eb_nucleon = Eb ./ (n_protons + n_neutrons);

end

