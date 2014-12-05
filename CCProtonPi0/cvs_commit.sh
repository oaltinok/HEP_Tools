cvs commit -m "v1_07 - 2014_12_05

General Updates:
    truthIsPlausible() Updated
        using muonIsPlausible()
    Removed reco_* variables (unused variables)
    Event Kinematics Calculations are Implemented
    Event Kinematics Need Update: See @todo in function header

New Functions:
    findLeadingProton()
    setLeadingProton_4P()


tagTruth() Updated
    Using GENIE:Event Record
        Creating Tables for Final State Particles
        Saving Information for Pi0, Pi0 Mother and Pi0 GrandMother
    setSignalKinematics() Updated
        Using FS Particle Trajectories
        Saving only 4-Momentum and Theta wrt Beam
        Removed all others

Vertex Reconstruction Updates:
    smear true vertex no longer used

Muon Reconstruction Updates:
    Muon MINOS Match Information now saved under CCProtonPi0_muon_ variables
        muon_hasMinosMatchTrack
        muon_hasMinosMatchStub

Proton Reconstruction Updates:
    Proton reconstruction is controlled with a new parameter: 
        m_reconstruct_NoProtonEvents
    When it is TRUE, Reconstruction Continues even there is no proton candidate" .