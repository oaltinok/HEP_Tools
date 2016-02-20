/*
   See CCProtonPi0.h header for Class Information
   */
#ifndef CCProtonPi0_Truth_cpp 
#define CCProtonPi0_Truth_cpp 1

#include "CCProtonPi0.h"

bool CCProtonPi0::isTrueVertexFiducial(Minerva::GenMinInteraction* truthEvent) const
{
    bool isFidVol = FiducialPointTool->isFiducial(truthEvent, m_fidHexApothem, m_fidUpStreamZ, m_fidDownStreamZ );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isFidVol", isFidVol );
    return isFidVol;
}

//------------------------------------------------------------------------------
// tagSignal() 
//      Tags Event as Signal or NOT
//      Signal: CC-Neutrino Interaction inside Fiducial Volume, 
//              FS Particles: Single Pi0, No Other (Proton and Neutrons are OK)
//      returns TRUE if it is a Signal Event    
//------------------------------------------------------------------------------
bool CCProtonPi0::tagSignal(Minerva::GenMinInteraction* truthEvent) const
{
    //--------------------------------------------------------------------------
    // get TG4Trajectories
    //--------------------------------------------------------------------------
    const SmartRefVector<Minerva::TG4Trajectory> pri_trajectories = truthEvent->trajectories();
    SmartRefVector<Minerva::TG4Trajectory>::const_iterator it_mcpart;

    if (pri_trajectories.size() != truthEvent->fSpdg().size() && pri_trajectories.size() + 1 != truthEvent->fSpdg().size() ) {
        warning()<<"Number of GENIE FS particles "<<truthEvent->fSpdg().size()<<
            " does not match number of Geant4 primary trajectories "<<pri_trajectories.size()<<"!"<<endmsg;

        for ( SmartRefVector<Minerva::TG4Trajectory>::const_iterator itTraj = pri_trajectories.begin(); itTraj != pri_trajectories.end(); ++itTraj ) {
            warning()<<"    G4 particle: "<<(*itTraj)->GetPDGCode()<<endmsg;
        }
        for ( std::vector<int>::const_iterator itGENIE = truthEvent->fSpdg().begin(); itGENIE != truthEvent->fSpdg().end(); ++itGENIE ) {
            warning()<<"    GENIE part:  "<<(*itGENIE)<<endmsg;
        }
    }

    //--------------------------------------------------------------------------
    // Count the Number of FS Particles
    //--------------------------------------------------------------------------
    int particle_PDG;
    int N_FSParticles = 0;
    int N_proton    = 0;
    int N_neutron   = 0;
    int N_pi0       = 0;
    int N_muon      = 0;
    int N_other     = 0;

    // Get Number of Final State Particles
    N_FSParticles = pri_trajectories.size();

    for (it_mcpart = pri_trajectories.begin(); it_mcpart != pri_trajectories.end(); ++it_mcpart) {
        particle_PDG =  (*it_mcpart)->GetPDGCode();
        if (  particle_PDG == PDG::proton ) N_proton++;
        else if ( particle_PDG == PDG::neutron ) N_neutron++;
        else if ( particle_PDG == PDG::pi0 ) N_pi0++;
        else if ( particle_PDG == PDG::muon ) N_muon++;
        else{
            N_other++;
            debug()<<"tagSignal Other PDG = "<<particle_PDG<<endmsg;
        }
    }
    debug()<<"N(Final State Particles) = "<<N_FSParticles<<endmsg;
    debug()<<"N(Proton) = "<<N_proton<<" N(Neutron) = "<<N_neutron<<" N(Pi0) = "<<N_pi0<<" N(Other) = "<<N_other<<endmsg;

    truthEvent->setIntData("N_FSParticles",  N_FSParticles);
    truthEvent->setIntData("N_proton",  N_proton);
    truthEvent->setIntData("N_pi0",  N_pi0);
    truthEvent->setIntData("N_other",  N_other);

    //--------------------------------------------------------------------------
    // Find Signal --  Inside Fiducial Volume
    //                 CC Neutrino Interaction
    //                 FS Particles: muon, pi0, X (No Meson)
    //--------------------------------------------------------------------------
    int t_current = truthEvent->current();
    int t_neutrinoPDG = truthEvent->incoming();
    double Enu = truthEvent->incomingEnergy();
    const double max_Enu = 20000; 

    bool isFidVol = truthEvent->filtertaglist()->isFilterTagTrue("isFidVol");
    bool isCCNeutrino = (t_current == 1 && t_neutrinoPDG == PDG::nu_mu);
    bool isFSGood = (N_pi0 == 1 && N_other == 0);
    bool isNeutrinoEnergyLow = Enu <= max_Enu; 

    bool isSignal = false;
    if(isFidVol && isCCNeutrino && isFSGood && isNeutrinoEnergyLow){
        isSignal = true;
        debug()<<"Found a Signal Event!"<<endmsg;
    }

    truthEvent->filtertaglist()->setOrAddFilterTag( "isSignal", isSignal );

    return isSignal;
}

void CCProtonPi0::tagBackgroundWithPi0(Minerva::GenMinInteraction* truthEvent) const
{
    Minerva::TG4Trajectories::iterator it_mcpart;
    Minerva::TG4Trajectories* trajects = NULL;
    if( exist<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default ) ) {
        trajects = get<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default );
    } else {
        warning() << "Could not get TG4Trajectories from " << Minerva::TG4TrajectoryLocation::Default << endmsg;
        return;
    }

    // just to make sure
    if( !trajects ) return;

    std::vector<int> primary_Pi0_ID;
    std::vector<int> primary_particle_ID;
    int particle_PDG = 0;
    int particle_ID = 0;
    int mother_ID = 0;

    int nPi0_Primary = 0;
    int nPi0_Secondary = 0;

    // Loop over Primary Trajectories to Count Primary Pi0
    for (it_mcpart = trajects->begin(); it_mcpart != trajects->end(); ++it_mcpart) {

        // Skip NonPrimary Trajectories
        if ( (*it_mcpart)->GetProcessName() != "Primary") continue;

        particle_PDG = (*it_mcpart)->GetPDGCode();
        particle_ID = (*it_mcpart)->GetTrackId();
        primary_particle_ID.push_back(particle_ID);

        if (particle_PDG == PDG::pi0){
            nPi0_Primary++; 
            primary_Pi0_ID.push_back(particle_ID);
        }
    }

    // Loop over Secondary Trajectories to Check Secondary Pi0
    for (it_mcpart = trajects->begin(); it_mcpart != trajects->end(); ++it_mcpart) {

        // Skip Primary Trajectories
        if ( (*it_mcpart)->GetProcessName() == "Primary") continue;

        particle_PDG = (*it_mcpart)->GetPDGCode();
        mother_ID = (*it_mcpart)->GetParentId();

        // Loop over only Secondary Particles
        if ( isMotherPrimary(primary_particle_ID, mother_ID) ){
            // Check Trajectory Mother is NOT a Primary Pi0 - Scattering
            if ( particle_PDG == PDG::pi0 && !isMotherPrimary(primary_Pi0_ID, mother_ID) ){
                nPi0_Secondary++;
            }
        }
    }

    // Get Total Number of Pi0 in Final State
    int nPi0_Total = nPi0_Primary + nPi0_Secondary;

    // Set Background Type
    bool NoPi0 = false;
    bool SinglePi0 = false;
    bool MultiPi0 = false;
    if (nPi0_Total == 0) NoPi0 = true;
    if (nPi0_Total == 1) SinglePi0 = true;
    if (nPi0_Total > 1) MultiPi0 = true;

    truthEvent->setIntData("Bckg_nPi0_Primary",nPi0_Primary);
    truthEvent->setIntData("Bckg_nPi0_Secondary",nPi0_Secondary);
    truthEvent->setIntData("Bckg_nPi0_Total",nPi0_Total);
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_NoPi0", NoPi0 );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_SinglePi0", SinglePi0 );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_MultiPi0", MultiPi0 );
}

//------------------------------------------------------------------------------
// tagBackground() 
//      Tags Background Event
//      Signal: CC-Neutrino Interaction inside Fiducial Volume, 
//              FS Particles: Single Pi0 + X (No Meson)
//------------------------------------------------------------------------------
void CCProtonPi0::tagBackground(Minerva::GenMinInteraction* truthEvent) const
{
    Minerva::TG4Trajectories::iterator it_mcpart;
    Minerva::TG4Trajectories* trajects = NULL;
    if( exist<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default ) ) {
        trajects = get<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default );
    } else {
        warning() << "Could not get TG4Trajectories from " << Minerva::TG4TrajectoryLocation::Default << endmsg;
        return;
    }

    // just to make sure
    if( !trajects ) return;

    int particle_PDG = -1;
    int nAntiMuon = 0;
    int nMuon = 0;
    int nNucleon = 0;
    int nPiCharged = 0;
    int nPiCharged_ChargeExchanged = 0;
    int nPiZero = 0;
    int nOther = 0;
    std::vector<int> charged_pion_ID;

    // Count FS Primary Particles 
    for (it_mcpart = trajects->begin(); it_mcpart != trajects->end(); ++it_mcpart) {

        // Skip Secondary Particles
        if ( (*it_mcpart)->GetProcessName() != "Primary") continue;

        particle_PDG = (*it_mcpart)->GetPDGCode();

        // Count Primary Particles
        if ( particle_PDG == -(PDG::muon) ) nAntiMuon++;
        else if (particle_PDG == PDG::muon) nMuon++;
        else if ( (particle_PDG == PDG::proton) || (particle_PDG == PDG::neutron) ) nNucleon++;
        else if ( std::abs(particle_PDG) == PDG::pi ) nPiCharged++;
        else if ( particle_PDG == PDG::pi0 ) nPiZero++;
        else{ 
            nOther++;
            debug()<<"Other PDG =  "<<particle_PDG<<endmsg;
        }
        // Save Charged Pion IDs for Charge Exchange
        if (std::abs(particle_PDG) == PDG::pi){
            charged_pion_ID.push_back((*it_mcpart)->GetTrackId());
        }
    } 

    int mother_ID = 0;
    // Loop Again to Find Charged Exchanged Pions
    if (nPiCharged > 0){
        for (it_mcpart = trajects->begin(); it_mcpart != trajects->end(); ++it_mcpart) {

            // Skip Primary Particles
            if ( (*it_mcpart)->GetProcessName() == "Primary") continue;

            particle_PDG = (*it_mcpart)->GetPDGCode();
            mother_ID = (*it_mcpart)->GetParentId();

            // Check Pi0 Mothers 
            if ( particle_PDG == PDG::pi0 ){
                if (isMotherPrimary(charged_pion_ID,mother_ID)){
                    nPiCharged_ChargeExchanged++;
                    nPiCharged--;
                    nPiZero++;
                }
            }

            if (nPiCharged == 0) break;
        } 
    }

    // -------------------------------------------------------------------------
    // Background Types 
    // -------------------------------------------------------------------------
    bool isNC = isInteractionNC(truthEvent);
    bool isAntiNeutrino = false;
    bool isQELike = false;
    bool isSingleChargedPion = false;
    bool isSingleChargedPion_ChargeExchanged = false;
    bool isDoublePionWithPi0 = false;
    bool isDoublePionWithoutPi0 = false;
    bool isMultiPionWithPi0 = false;
    bool isMultiPionWithoutPi0 = false;
    bool isOther = false;

    if (isNC) ; // Do Nothing
    else if (nAntiMuon == 1) isAntiNeutrino = true;
    else if (nOther > 0) isOther = true;
    else if (nPiZero == 0 && nPiCharged == 0) isQELike = true;
    else if (nPiZero == 0 && nPiCharged == 1) isSingleChargedPion = true;
    else if (nPiZero == 1 && nPiCharged == 0 && nPiCharged_ChargeExchanged == 1) isSingleChargedPion_ChargeExchanged = true;
    else if (nPiZero == 1 && nPiCharged == 1) isDoublePionWithPi0 = true;
    else if (nPiZero == 2 && nPiCharged == 0) isDoublePionWithPi0 = true;
    else if (nPiZero == 0 && nPiCharged == 2) isDoublePionWithoutPi0 = true;
    else if (nPiZero > 0 && (nPiCharged + nPiZero) > 2) isMultiPionWithPi0 = true;
    else if (nPiZero == 0 && nPiCharged > 2) isMultiPionWithoutPi0 = true;
    else{ 
        isOther = true;
        debug()<<"No Background Type Matched!"<<endmsg;
    }
    // -------------------------------------------------------------------------
    // Check for Background Branching: hasMichel
    // -------------------------------------------------------------------------
    bool hasMichel = false;
    hasMichel = checkMichel(truthEvent);

    // -------------------------------------------------------------------------
    // Debugging
    // -------------------------------------------------------------------------
    debug()<<"Background Type"<<endmsg;
    debug()<<"nAntiMuon = "<<nAntiMuon<<endmsg;
    debug()<<"nOther = "<<nOther<<endmsg;
    debug()<<"nPiZero = "<<nPiZero<<endmsg;
    debug()<<"nPiCharged = "<<nPiCharged<<endmsg;
    debug()<<"nPiCharged_ChargeExchanged = "<<nPiCharged_ChargeExchanged<<endmsg;
    debug()<<"---------"<<endmsg;
    debug()<<"isNC = "<<isNC<<endmsg;
    debug()<<"isAntiNeutrino = "<<isAntiNeutrino<<endmsg;
    debug()<<"isQELike = "<<isQELike<<endmsg;
    debug()<<"isSingleChargedPion = "<<isSingleChargedPion<<endmsg;
    debug()<<"isSingleChargedPion_ChargeExchanged = "<<isSingleChargedPion_ChargeExchanged<<endmsg;
    debug()<<"isDoublePionWithPi0 = "<<isDoublePionWithPi0<<endmsg;
    debug()<<"isDoublePionWithoutPi0 = "<<isDoublePionWithoutPi0<<endmsg;
    debug()<<"isMultiPionWithPi0 = "<<isMultiPionWithPi0<<endmsg;
    debug()<<"isMultiPionWithoutPi0 = "<<isMultiPionWithoutPi0<<endmsg;
    debug()<<"isOther = "<<isOther<<endmsg;
    debug()<<"---------"<<endmsg;
    // -------------------------------------------------------------------------

    truthEvent->setIntData("Bckg_nOther",nOther);
    truthEvent->setIntData("Bckg_nPiCharged",nPiCharged);
    truthEvent->setIntData("Bckg_nPiCharged_ChargeExchanged",nPiCharged_ChargeExchanged);
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_NC", isNC );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_AntiNeutrino", isAntiNeutrino );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_QELike", isQELike );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_SingleChargedPion", isSingleChargedPion );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_SingleChargedPion_ChargeExchanged", isSingleChargedPion_ChargeExchanged );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_DoublePionWithPi0", isDoublePionWithPi0 );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_DoublePionWithoutPi0", isDoublePionWithoutPi0 );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_MultiPionWithPi0", isMultiPionWithPi0 );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_MultiPionWithoutPi0", isMultiPionWithoutPi0 );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_Other", isOther );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_withMichel", hasMichel );
}

//------------------------------------------------------------------------------
// checkMichel() 
//      Check MC Trajectories for a PiPlus to MuPlus Decay
//      Returns TRUE if such a decay exists and fills Michel Truth Branches
//------------------------------------------------------------------------------
bool CCProtonPi0::checkMichel(Minerva::GenMinInteraction* truthEvent) const
{
    bool hasMichel = false;
    int N_trueMichels = 0;

    Minerva::TG4Trajectories::iterator it_mcpart;
    Minerva::TG4Trajectories* trajects = NULL;
    if( exist<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default ) ) {
        trajects = get<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default );
    } else {
        warning() << "Could not get TG4Trajectories from " << Minerva::TG4TrajectoryLocation::Default << endmsg;
        return false;
    }

    // just to make sure
    if( !trajects ) return false;

    // Save Vertex Information
    Gaudi::LorentzVector vtx = truthEvent->Vtx();

    int particle_PDG;
    int particle_ID;
    int mother_ID;
    std::vector<int> PiPlus_traj_ID;
    std::vector<double> PiPlus_traj_length;
    std::vector<double> PiPlus_dist_vtx;
    std::vector<double> PiPlus_momentum;

    // Save all Pi+ Info
    for (it_mcpart = trajects->begin(); it_mcpart != trajects->end(); ++it_mcpart) {

        particle_PDG =  (*it_mcpart)->GetPDGCode();
        particle_ID = (*it_mcpart)->GetTrackId();

        // Save Pi+ Information
        if (particle_PDG == PDG::pi ){
            PiPlus_traj_ID.push_back(particle_ID);
            Gaudi::LorentzVector piplus_momentum = (*it_mcpart)->GetInitialMomentum();
            Gaudi::LorentzVector piplus_start = (*it_mcpart)->GetInitialPosition();
            Gaudi::LorentzVector piplus_end = (*it_mcpart)->GetFinalPosition();

            double piplus_length = calcDistance(piplus_start.X(), piplus_start.Y(), piplus_start.Z(),
                    piplus_end.X(), piplus_end.Y(), piplus_end.Z());

            // PiPlus Initial Point Distance to Vertex
            double piplus_dist_vtx = calcDistance(  vtx.X(), vtx.Y(), vtx.Z(),
                    piplus_start.X(), piplus_start.Y(), piplus_start.Z());                                                

            PiPlus_traj_length.push_back(piplus_length);
            PiPlus_momentum.push_back(piplus_momentum.P());
            PiPlus_dist_vtx.push_back(piplus_dist_vtx);
        }
    }

    // Loop again to locate Mu-Plus as a daughter of Pi-Plus
    for (it_mcpart = trajects->begin(); it_mcpart != trajects->end(); ++it_mcpart) {

        particle_PDG =  (*it_mcpart)->GetPDGCode();
        mother_ID = (*it_mcpart)->GetParentId();


        // Check particle is a Mu-Plus
        if (particle_PDG == -(PDG::muon) ){

            // Check Mu-Plus Mother is one of the Pi-Plus
            int michelPion_ind = getMichelPion(PiPlus_traj_ID, mother_ID);
            if ( michelPion_ind != -1 ){

                Gaudi::LorentzVector michelMuon_end = (*it_mcpart)->GetFinalPosition();

                // Skip if the michel Muon end point is NOT in Reconstructable Volume
                if ( !FiducialPointTool->isFiducial(michelMuon_end.X(),michelMuon_end.Y(),michelMuon_end.Z(), m_recoHexApothem, m_recoUpStreamZ, m_recoDownStreamZ ) ) {
                    continue;
                }

                hasMichel = true;

                // Save Michel Information Only for the 1st Michel
                N_trueMichels++;
                if (N_trueMichels > 1) continue;

                // -------------------------------------------------------------
                // Save Michel Electron Information
                // -------------------------------------------------------------
                int muon_ID = (*it_mcpart)->GetTrackId();
                saveMichelElectron(truthEvent, muon_ID);

                // -------------------------------------------------------------
                // Get Michel Muon Information
                // -------------------------------------------------------------
                std::vector<double> michelMuonPoint;
                Gaudi::LorentzVector michelMuon_momentum = (*it_mcpart)->GetInitialMomentum();
                Gaudi::LorentzVector michelMuon_start = (*it_mcpart)->GetInitialPosition();

                // Anti-Muon End Point
                michelMuonPoint.push_back(michelMuon_end.X());
                michelMuonPoint.push_back(michelMuon_end.Y());
                michelMuonPoint.push_back(michelMuon_end.Z());

                // Anti-Muon Length
                double michel_length = calcDistance(michelMuon_start.X(), michelMuon_start.Y(), michelMuon_start.Z(),
                        michelMuon_end.X(), michelMuon_end.Y(), michelMuon_end.Z());

                // Distance AntiMuon End point to Vertex
                double vtx_michel_end = calcDistance(   vtx.X(), vtx.Y(), vtx.Z(),
                        michelMuon_end.X(), michelMuon_end.Y(), michelMuon_end.Z());

                // Fill NTuple
                truthEvent->setContainerDoubleData("michelMuon_endPoint",michelMuonPoint);
                truthEvent->setDoubleData("michelMuon_P", michelMuon_momentum.P());
                truthEvent->setDoubleData("michelMuon_length", michel_length);
                truthEvent->setDoubleData("michelMuon_end_dist_vtx", vtx_michel_end);
                truthEvent->setDoubleData("michelPion_P", PiPlus_momentum[michelPion_ind]);
                truthEvent->setDoubleData("michelPion_length", PiPlus_traj_length[michelPion_ind]);
                truthEvent->setDoubleData("michelPion_begin_dist_vtx", PiPlus_dist_vtx[michelPion_ind]);

            }
        } 
    }

    debug()<<"N_trueMichels = "<<N_trueMichels<<endmsg;
    truthEvent->setIntData("N_trueMichelElectrons", N_trueMichels);
    return hasMichel;
}

//------------------------------------------------------------------------------
// writeBackgroundType() 
//      Prints Background Type with its Branching Information
//          Uses info()
//------------------------------------------------------------------------------
void CCProtonPi0::writeBackgroundType(Minerva::GenMinInteraction* truthEvent) const
{
    // RETURN if it is a Signal Event
    if(truthEvent->filtertaglist()->isFilterTagTrue("isSignal")) return;

    // Print Background Type
    info()<<"Background With Pi0"<<endmsg;
    if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_NoPi0")) info() <<"\tBackground without Pi0"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_SinglePi0")) info() <<"\tBackground with Single Pi0"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_MultiPi0")) info() <<"\tBackground with Multi Pi0"<<endmsg;
    else warning()<<"No BackgroundWithPi0!"<<endmsg;

    info()<<"Background Type"<<endmsg;
    if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_NC")) info()<<"\tBackground: NC"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_AntiNeutrino")) info() <<"\tBackground: AntiNeutrino"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_QELike")) info() <<"\tBackground: QE Like (No Pion)"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_SingleChargedPion")) info() <<"\tBackground: Single Charged Pion"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_SingleChargedPion_ChargeExchanged")) info() <<"\tBackground: Single Charged Pion (Charge Exchanged)"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_DoublePionWithPi0")) info() <<"\tBackground: Double Pion With Pi0"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_DoublePionWithoutPi0")) info() <<"\tBackground: Double Pion Without Pi0"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_MultiPionWithPi0")) info() <<"\tBackground: Multiple Pion With Pi0"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_MultiPionWithoutPi0")) info() <<"\tBackground: Multiple Pion Without Pi0"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_Other")) info() <<"\tBackground: Other"<<endmsg;
    else warning()<< "No Background Type!"<<endmsg;

    // Print Branching Information
    bool hasMichel = truthEvent->filtertaglist()->isFilterTagTrue("isBckg_withMichel");

    info()<<"Background Branching: "<<endmsg;
    info()<<"\tMichel = "<<hasMichel<<endmsg;
}

bool CCProtonPi0::isInteractionNC(Minerva::GenMinInteraction* truthEvent) const
{
    bool isNC = false;
    int current = truthEvent->current();

    if ( current == 2 ) isNC = true;
    else isNC = false;

    truthEvent->filtertaglist()->setOrAddFilterTag( "isNC", isNC );

    return isNC;
}

void CCProtonPi0::setPi0GenieRecord(Minerva::GenMinInteraction* truthEvent) const
{
    //--------------------------------------------------------------------------
    // Get Event Record
    //--------------------------------------------------------------------------
    const Minerva::GenMinEventRecord* eventRecord = truthEvent->eventRecord();
    const std::vector<int> eRecPartID = eventRecord->eRecPartID();
    const std::vector<int> eRecPartStatus = eventRecord->eRecPartStatus();
    const std::vector<int> eRecMother = eventRecord->eRecMother();
    double t_pi0_status = -9;
    int t_pi0_Mother = -9;
    int t_pi0_MotherStatus= -9;
    int t_pi0_GrandMother= -9;
    int t_pi0_GrandMotherStatus= -9;
    int ind_Mother;
    int ind_GrandMother;

    for (unsigned int g_part = 0; g_part < eRecPartID.size(); g_part++){
        if (eRecPartID[g_part] == 111){
            // Set Pi0 Status - Initial Status
            t_pi0_status = eRecPartStatus[g_part];

            // Set Pi0 Mother
            ind_Mother = eRecMother[g_part];
            t_pi0_Mother = eRecPartID[ind_Mother];
            t_pi0_MotherStatus = eRecPartStatus[ind_Mother];

            // Set Pi0 GrandMother
            ind_GrandMother = eRecMother[ind_Mother];
            t_pi0_GrandMother = eRecPartID[ind_GrandMother];
            t_pi0_GrandMotherStatus = eRecPartStatus[ind_GrandMother];

            break;
        }
    }

    if (t_pi0_Mother == 2000000001){
        debug()<<"Setting pi0_Mother to = "<<-5<<endmsg; 
        t_pi0_Mother = -5;
    }

    truthEvent->setIntData("pi0_status", t_pi0_status);
    truthEvent->setIntData("pi0_Mother", t_pi0_Mother);
    truthEvent->setIntData("pi0_MotherStatus", t_pi0_MotherStatus);
    truthEvent->setIntData("pi0_GrandMother", t_pi0_GrandMother);
    truthEvent->setIntData("pi0_GrandMotherStatus", t_pi0_GrandMotherStatus);

}

void CCProtonPi0::setSignalKinematics(Minerva::GenMinInteraction* truthEvent) const
{
    int Pi0_ID = -1; // Pi0_ID will be updated while setting primary trajectory kinematics

    setSignal_PrimaryTrajectoryKinematics(truthEvent, Pi0_ID);
    setSignal_SecondaryTrajectoryKinematics(truthEvent,Pi0_ID);
}

void CCProtonPi0::setSignal_PrimaryTrajectoryKinematics(Minerva::GenMinInteraction *truthEvent, int &Pi0_ID) const
{
    const double SENTINEL = -9.9;

    // Get MC Trajectories
    Minerva::TG4Trajectories* trajects = NULL;
    Minerva::TG4Trajectories::iterator it_mcpart;
    if( exist<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default ) ) {
        trajects = get<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default );
    } else {
        warning() << "Could not get TG4Trajectories from " << Minerva::TG4TrajectoryLocation::Default << endmsg;
        return;
    }

    // just to make sure
    if( !trajects ) return;

    // Information to be Saved
    std::vector<double> muon_4P(4,SENTINEL);
    std::vector<double> pi0_4P(4,SENTINEL);
    std::vector<double> proton_4P(4,SENTINEL);
    
    double muon_P;
    double muon_theta;
    double pi0_P;
    double pi0_theta;
    double proton_P;
    double proton_theta;

    // Loop Vars
    Gaudi::LorentzVector temp_4P;
    int particle_PDG;
    int particle_ID;
    int max_protonE = 0.0;

    // Loop Over Primary Trajectories
    for (it_mcpart = trajects->begin(); it_mcpart != trajects->end(); ++it_mcpart) {

        if ( (*it_mcpart)->GetProcessName() != "Primary") continue;

        particle_ID = (*it_mcpart)->GetTrackId();
        particle_PDG = (*it_mcpart)->GetPDGCode();
        temp_4P = (*it_mcpart)->GetInitialMomentum();
        
        if (particle_PDG == PDG::muon){
            muon_P = temp_4P.P();
            muon_theta = temp_4P.theta();
            muon_4P[0] = temp_4P.px();  
            muon_4P[1] = temp_4P.py(); 
            muon_4P[2] = temp_4P.pz(); 
            muon_4P[3] = temp_4P.E(); 
            debug()<<"True Muon 4P = "<<temp_4P<<endmsg;
        }else if (particle_PDG == PDG::pi0){
            Pi0_ID = particle_ID; // Update Pi0_ID
            pi0_P = temp_4P.P();
            pi0_theta = temp_4P.theta();
            pi0_4P[0] = temp_4P.px();  
            pi0_4P[1] = temp_4P.py(); 
            pi0_4P[2] = temp_4P.pz(); 
            pi0_4P[3] = temp_4P.E();
            debug()<<"True Pi0 4P = "<<temp_4P<<endmsg;
        }else if ((particle_PDG == PDG::proton) && (max_protonE < temp_4P.E()) ){
            proton_P = temp_4P.P();
            proton_theta = temp_4P.theta();
            proton_4P[0] = temp_4P.px();  
            proton_4P[1] = temp_4P.py(); 
            proton_4P[2] = temp_4P.pz(); 
            proton_4P[3] = temp_4P.E();  
            // Update max_ProtonE
            max_protonE = temp_4P.E();
            debug()<<"True Proton 4P = "<<temp_4P<<endmsg;
        }
    }

    // Save to NTuples
    truthEvent->setContainerDoubleData("muon_4P", muon_4P);
    truthEvent->setContainerDoubleData("pi0_4P", pi0_4P);
    truthEvent->setContainerDoubleData("proton_4P", proton_4P);
    truthEvent->setDoubleData("muon_P", muon_P);
    truthEvent->setDoubleData("pi0_P", pi0_P);
    truthEvent->setDoubleData("proton_P", proton_P);
    truthEvent->setDoubleData("muon_theta", muon_theta);
    truthEvent->setDoubleData("pi0_theta", pi0_theta);
    truthEvent->setDoubleData("proton_theta", proton_theta);
}

void CCProtonPi0::setSignal_SecondaryTrajectoryKinematics(Minerva::GenMinInteraction *truthEvent, int &Pi0_ID) const
{
    const double SENTINEL = -9.9;

    // Get MC Trajectories
    Minerva::TG4Trajectories* trajects = NULL;
    Minerva::TG4Trajectories::iterator it_mcpart;
    if( exist<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default ) ) {
        trajects = get<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default );
    } else {
        warning() << "Could not get TG4Trajectories from " << Minerva::TG4TrajectoryLocation::Default << endmsg;
        return;
    }

    // just to make sure
    if( !trajects ) return;

    // Information to be Saved
    std::vector<double> gamma1_4P(4,SENTINEL);
    std::vector<double> gamma2_4P(4,SENTINEL);
    std::vector<double> gamma1_init_pos(3,SENTINEL);
    std::vector<double> gamma1_final_pos(3,SENTINEL);
    std::vector<double> gamma2_init_pos(3,SENTINEL);
    std::vector<double> gamma2_final_pos(3,SENTINEL);
    bool isGamma1_conv_inside = false;
    bool isGamma2_conv_inside = false;

    // Loop Vars
    Gaudi::LorentzVector temp_4P;
    Gaudi::LorentzVector temp_init_pos;
    Gaudi::LorentzVector temp_final_pos;
    int particle_PDG;
    int mother_ID;
    int nGammas = 0;

    for (it_mcpart = trajects->begin(); it_mcpart != trajects->end(); ++it_mcpart) {

        if ( nGammas == 2) break;
        if ( (*it_mcpart)->GetProcessName() == "Primary") continue;

        temp_4P = (*it_mcpart)->GetInitialMomentum();
        temp_init_pos = (*it_mcpart)->GetInitialPosition();
        temp_final_pos = (*it_mcpart)->GetFinalPosition();
        particle_PDG = (*it_mcpart)->GetPDGCode();
        mother_ID = (*it_mcpart)->GetParentId();

        if (mother_ID == Pi0_ID && particle_PDG == PDG::gamma) {
            if (nGammas == 0){
                gamma1_4P[0] = temp_4P.px();  
                gamma1_4P[1] = temp_4P.py(); 
                gamma1_4P[2] = temp_4P.pz(); 
                gamma1_4P[3] = temp_4P.E(); 

                gamma1_init_pos[0] = temp_init_pos.px();  
                gamma1_init_pos[1] = temp_init_pos.py(); 
                gamma1_init_pos[2] = temp_init_pos.pz(); 

                gamma1_final_pos[0] = temp_final_pos.px();  
                gamma1_final_pos[1] = temp_final_pos.py(); 
                gamma1_final_pos[2] = temp_final_pos.pz(); 
                isGamma1_conv_inside = FiducialPointTool->isFiducial(temp_final_pos.x(), temp_final_pos.y(), temp_final_pos.z(), m_recoHexApothem, m_recoUpStreamZ, m_recoDownStreamZ);

                nGammas++;
            }else if (nGammas == 1){
                gamma2_4P[0] = temp_4P.px();  
                gamma2_4P[1] = temp_4P.py(); 
                gamma2_4P[2] = temp_4P.pz(); 
                gamma2_4P[3] = temp_4P.E(); 

                gamma2_init_pos[0] = temp_init_pos.px();  
                gamma2_init_pos[1] = temp_init_pos.py(); 
                gamma2_init_pos[2] = temp_init_pos.pz(); 

                gamma2_final_pos[0] = temp_final_pos.px();  
                gamma2_final_pos[1] = temp_final_pos.py(); 
                gamma2_final_pos[2] = temp_final_pos.pz(); 
                isGamma2_conv_inside = FiducialPointTool->isFiducial(temp_final_pos.x(), temp_final_pos.y(), temp_final_pos.z(), m_recoHexApothem, m_recoUpStreamZ, m_recoDownStreamZ);

                nGammas++;
            }
        }   
    }

    // Make sure Gamma1 is the more energetic one
    if (gamma1_4P[3] < gamma2_4P[3]){
        gamma1_4P.swap(gamma2_4P);
        gamma1_init_pos.swap(gamma2_init_pos);
        gamma1_final_pos.swap(gamma2_final_pos);
    }

    // Fill NTuples
    truthEvent->setContainerDoubleData("gamma1_4P",  gamma1_4P);
    truthEvent->setContainerDoubleData("gamma2_4P",  gamma2_4P);
    truthEvent->setContainerDoubleData("gamma1_init_pos",  gamma1_init_pos);
    truthEvent->setContainerDoubleData("gamma2_init_pos",  gamma2_init_pos);
    truthEvent->setContainerDoubleData("gamma1_final_pos",  gamma1_final_pos);
    truthEvent->setContainerDoubleData("gamma2_final_pos",  gamma2_final_pos);
    truthEvent->filtertaglist()->setOrAddFilterTag("isGamma1_conv_inside", isGamma1_conv_inside );
    truthEvent->filtertaglist()->setOrAddFilterTag("isGamma2_conv_inside", isGamma2_conv_inside );
}


//------------------------------------------------------------------------------
// writeFSParticleTable() writes Final State Particle Table using truth information
//------------------------------------------------------------------------------
void CCProtonPi0::writeFSParticleTable(bool isSignal) const
{     
    Minerva::TG4Trajectories* trajects = NULL;
    if( exist<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default ) ) {
        trajects = get<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default );
    } else {
        warning() << "Could not get TG4Trajectories from " << Minerva::TG4TrajectoryLocation::Default << endmsg;
        return;
    }

    // just to make sure
    if( !trajects ) return;

    Gaudi::LorentzVector temp_4p;
    double gammaE_threshold = 20; //MeV
    int particle_PDG;
    int particle_ID;
    int mother_ID = -1;
    std::vector<int> primary_traj_ID;
    std::vector<int> pri_piplus_ID;
    std::vector<int> pri_piminus_ID;


    // Using Reverse Iterator to have the table in correct order
    Minerva::TG4Trajectories::reverse_iterator rit_mcpart;

    std::cout <<"==============================================================="<<std::endl;
    std::cout <<"Final State Particle Table"<<std::endl;
    if(isSignal){
        std::cout<<">> Marked as Signal! <<"<<std::endl;
    }
    std::cout<<std::left;
    std::cout.width(12); std::cout<<"ID";
    std::cout.width(24); std::cout<<"Process";
    std::cout.width(12); std::cout<<"PDG";
    std::cout.width(12); std::cout<<"Mother";
    std::cout.width(12); std::cout<<"Momentum";
    std::cout.width(12); std::cout<<"Energy";
    std::cout<<std::endl;

    // Print all Primary Particles
    for (rit_mcpart = trajects->rbegin(); rit_mcpart != trajects->rend(); ++rit_mcpart) {

        if ( (*rit_mcpart)->GetProcessName() == "Primary"){
            temp_4p = (*rit_mcpart)->GetInitialMomentum();
            particle_PDG = (*rit_mcpart)->GetPDGCode();
            particle_ID = (*rit_mcpart)->GetTrackId();

            primary_traj_ID.push_back(particle_ID);

            std::cout.width(12); std::cout<<particle_ID;    // ID
            std::cout.width(24); std::cout<<(*rit_mcpart)->GetProcessName();// Process
            std::cout.width(12); std::cout<<particle_PDG;                   // PDG
            std::cout.width(12); std::cout<<(*rit_mcpart)->GetParentId();   // Mother PDG
            std::cout.width(12); std::cout<<temp_4p.P();
            std::cout.width(12); std::cout<<temp_4p.E();
            std::cout<<std::endl;

            // Save Primary Trajectory Information
            primary_traj_ID.push_back(particle_ID);
        }

    }


    // Print Special 2nd Particles
    for (rit_mcpart = trajects->rbegin(); rit_mcpart != trajects->rend(); ++rit_mcpart) {

        // Skip Primary Trajectories
        if ( (*rit_mcpart)->GetProcessName() == "Primary") continue;

        temp_4p = (*rit_mcpart)->GetInitialMomentum();
        particle_PDG =  (*rit_mcpart)->GetPDGCode();
        mother_ID = (*rit_mcpart)->GetParentId();

        if ( isMotherPrimary(primary_traj_ID, mother_ID) ){
            if (particle_PDG == -(PDG::muon) ){
                std::cout<<"Michel!"<<std::endl;
                std::cout.width(12); std::cout<<(*rit_mcpart)->GetTrackId();    // ID
                std::cout.width(24); std::cout<<(*rit_mcpart)->GetProcessName();// Process
                std::cout.width(12); std::cout<<particle_PDG;                   // PDG
                std::cout.width(12); std::cout<<(*rit_mcpart)->GetParentId();   // Mother PDG
                std::cout.width(12); std::cout<<temp_4p.P();
                std::cout.width(12); std::cout<<temp_4p.E();
                std::cout<<std::endl;
            }else if(particle_PDG == PDG::pi0){
                std::cout<<"Secondary Pi0!"<<std::endl;
                std::cout.width(12); std::cout<<(*rit_mcpart)->GetTrackId();    // ID
                std::cout.width(24); std::cout<<(*rit_mcpart)->GetProcessName();// Process
                std::cout.width(12); std::cout<<particle_PDG;                   // PDG
                std::cout.width(12); std::cout<<(*rit_mcpart)->GetParentId();   // Mother PDG
                std::cout.width(12); std::cout<<temp_4p.P();
                std::cout.width(12); std::cout<<temp_4p.E();
                std::cout<<std::endl;
            }else if(particle_PDG == PDG::gamma && temp_4p.E() > gammaE_threshold  ){
                std::cout<<"Secondary Gamma with Energy > 20 MeV!"<<std::endl;
                std::cout.width(12); std::cout<<(*rit_mcpart)->GetTrackId();    // ID
                std::cout.width(24); std::cout<<(*rit_mcpart)->GetProcessName();// Process
                std::cout.width(12); std::cout<<particle_PDG;                   // PDG
                std::cout.width(12); std::cout<<(*rit_mcpart)->GetParentId();   // Mother PDG
                std::cout.width(12); std::cout<<temp_4p.P();
                std::cout.width(12); std::cout<<temp_4p.E();
                std::cout<<std::endl;
            }
        }
    }
    std::cout<<"--------------------------------------------------------------"<<std::endl;

}

void CCProtonPi0::writeEventRecord(Minerva::GenMinInteraction* truthEvent, bool isSignal) const
{
    //--------------------------------------------------------------------------
    // Get Event Record
    //--------------------------------------------------------------------------
    const Minerva::GenMinEventRecord* eventRecord = truthEvent->eventRecord();
    const std::vector<int> eRecPartID = eventRecord->eRecPartID();
    const std::vector<int> eRecPartStatus = eventRecord->eRecPartStatus();
    const std::vector<int> eRecMother = eventRecord->eRecMother();
    const std::vector<double> eRecPartPx = eventRecord->eRecMomentumPx();
    const std::vector<double> eRecPartPy = eventRecord->eRecMomentumPy();
    const std::vector<double> eRecPartPz = eventRecord->eRecMomentumPz();
    const std::vector<double> eRecPartE = eventRecord->eRecEnergy();
    double momentum = 0.0;

    std::cout <<"Event Record"<<std::endl;
    if(isSignal){
        std::cout<<">> Marked as Signal! <<"<<std::endl;
    }
    std::cout<<std::left;
    std::cout.width(4); std::cout<<"ID";
    std::cout.width(12); std::cout<<"PDG";
    std::cout.width(12); std::cout<<"Status";
    std::cout.width(12); std::cout<<"Mother";
    std::cout.width(12); std::cout<<"Momentum";
    std::cout.width(12); std::cout<<"Total E";
    std::cout<<std::endl;

    for(unsigned int g_part = 0; g_part < eRecPartID.size(); g_part++ ){
        // Get Momentum
        momentum = sqrt(pow(eRecPartPx[g_part],2) + 
                pow(eRecPartPy[g_part],2) + 
                pow(eRecPartPz[g_part],2) );

        // Write Table               
        std::cout.width(4); std::cout<<g_part;
        std::cout.width(12); std::cout<<eRecPartID[g_part];
        std::cout.width(12); std::cout<<eRecPartStatus[g_part];
        std::cout.width(12); std::cout<<eRecMother[g_part];
        std::cout.width(12); std::cout<<momentum;
        std::cout.width(12); std::cout<<eRecPartE[g_part];
        std::cout<<std::endl;
    }

    std::cout <<"==============================================================="<<std::endl;
}


//------------------------------------------------------------------------------
// setTargetMaterial() Stores Target Material using Truth information
//------------------------------------------------------------------------------
void CCProtonPi0::setTargetMaterial(Minerva::GenMinInteraction* truthEvent) const
{
    Gaudi::LorentzVector t_vtx = truthEvent->Vtx(); 
    int t_targetMaterial = -1; 

    G4ThreeVector G4vtx(t_vtx.x(), t_vtx.y(), t_vtx.z());
    G4ThreeVector G4vec(0, 0, 1); // The next function needs a direction vector as argument. Here's a dummy.
    G4Navigator* geantNav = new G4Navigator;
    geantNav->SetWorldVolume(m_gigaCnvSvc->world());
    G4VPhysicalVolume* physVol=geantNav->LocateGlobalPointAndSetup(G4vtx, &G4vec, false, true);
    G4Material* material=physVol->GetLogicalVolume()->GetMaterial(); 
    delete geantNav;

    if ( material->GetName() == "/dd/Materials/Minerva/PlasticScint" ) {
        t_targetMaterial = 0;
    }
    else if ( material->GetName() == "/dd/Materials/Minerva/PSTitaniumDioxide" ) {
        t_targetMaterial = 1;
    }
    else if ( material->GetName() == "/dd/Materials/Minerva/GreyEpoxy" ) {
        t_targetMaterial = 2;
    }
    else if ( material->GetName() == "/dd/Materials/Minerva/Lexan" ) {
        t_targetMaterial = 3;
    }
    else if ( material->GetName() == "/dd/Materials/Minerva/GreenFiber" ) {
        t_targetMaterial = 4;
    }
    else if ( material->GetName() == "/dd/Materials/Minerva/PureLead" ) {
        t_targetMaterial = 5;
    }
    else if ( material->GetName() == "/dd/Materials/Minerva/PureCarbon" ) {
        t_targetMaterial = 6;
    }  
    else if ( material->GetName() == "/dd/Materials/Minerva/PureAluminum" ) {
        t_targetMaterial = 7;
    }   
    else if ( material->GetName() == "/dd/Materials/Minerva/StainlessSteel" ) {
        t_targetMaterial = 8;
    }   
    else if ( material->GetName() == "/dd/Materials/Air" ) {
        t_targetMaterial = 9;
    }
    else t_targetMaterial = 10;    

    truthEvent->setIntData("target_material", t_targetMaterial);
}

//------------------------------------------------------------------------------
// isMotherPrimary() returns True if the particle mother is a Primary Particle
//------------------------------------------------------------------------------
bool CCProtonPi0::isMotherPrimary(std::vector<int>& motherList, int mother ) const
{
    for( unsigned int i = 0 ; i < motherList.size(); i++ ) {
        if( mother == motherList[i]) return true;
    }

    return false;
}

//------------------------------------------------------------------------------
// getMichelPion() returns the indice of Michel Pion or -1
//------------------------------------------------------------------------------
int CCProtonPi0::getMichelPion(std::vector<int>& piList, int ID ) const
{
    for( unsigned int i = 0 ; i < piList.size(); i++ ) {
        if( ID == piList[i]) return i;
    }

    return -1;
}

//------------------------------------------------------------------------------
// saveMichelElectron() - Saves Michel Electron Information
//------------------------------------------------------------------------------
void CCProtonPi0::saveMichelElectron(Minerva::GenMinInteraction* truthEvent, int muon_ID) const
{
    Minerva::TG4Trajectories::iterator it_mcpart;
    Minerva::TG4Trajectories* trajects = NULL;
    if( exist<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default ) ) {
        trajects = get<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default );
    } else {
        warning() << "Could not get TG4Trajectories from " << Minerva::TG4TrajectoryLocation::Default << endmsg;
        return;
    }

    // just to make sure
    if( !trajects ) return;

    // Loop Over all Trajectories to Find Michel Electrons
    for (it_mcpart = trajects->begin(); it_mcpart != trajects->end(); ++it_mcpart) {

        int particle_PDG =  (*it_mcpart)->GetPDGCode();
        int mother_ID = (*it_mcpart)->GetParentId();

        if (particle_PDG == -(PDG::electron) && mother_ID == muon_ID ){

            Gaudi::LorentzVector michelElectron_4P = (*it_mcpart)->GetInitialMomentum();
            truthEvent->setDoubleData("michelElectron_P", michelElectron_4P.P());
            truthEvent->setDoubleData("michelElectron_E", michelElectron_4P.E());

            debug()<<"Found Michel Electron with PDG = "<<particle_PDG<<" E = "<<michelElectron_4P.E()<<endmsg;

            break;
        }
    }
}

#endif

