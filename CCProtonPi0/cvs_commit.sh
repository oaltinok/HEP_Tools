cvs commit -m "v1_2 2014_07_09
TG4Trajectory Map Removed
Muon Reconstruction
	Using Global Variables
		m_MuonTrack
		m_MuonProng
		m_MuonParticle
	SetMuonParticleData() Improved
	
Proton Reconstruction
	Fixed a bug causing empty NTuple Branches even with a successful proton reconstruction
	Using Global Variables
		m_ProtonProngs
		m_ProtonParticles
	Stop Algorithm if CreateTrackedParticles() fails
		New Variable Cut_Particle_None
	Removed Proton Score Cut
	getProtonProng() Improved
	SetProtonParticleData() Improved
	Branches initialized to SENTINEL = -9.9
	
Global Variable: m_PrimaryVertex
	Still need to pass the variable to CCPi0 Functions - Some of the functions outside the Global Scope
	
Found a new bug in correctProtonProngEnergy() 
	If the function changes proton energy, NTuple Branches for Momentum and Energy does not filled correctly
	Will fix this with v1_3" .