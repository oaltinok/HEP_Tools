CCPROTONPI0_V="v2_97"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:

---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Important Bug Fixes
		FillLatErrorBand_ProtonEnergy_Birks_invMass()
			I was not filling 1 track events 

		FillLatErrorBand_EM_EnergyScale_invMass()
			I was using zero shift, it should propagate the actual shift on pi0 Inv Mass
			I was filling the histograms before invMass selection, it should be after

		BckgConstraint Err Band
			I was using a single error band with 2 universes for all three bckg constraints 
			Now I am using 3 error bands with 2 universes each. 
				Each error band associated with specific bckg constraint
				This method ensures that the errors are quadrature summed
		
		W Selection included in Lateral Error Band Fills
			W is changing whenever I change, muon and pi0 kinematics
			W < 1.8 must be satisfied in all universes.

		Muon Energy shifts calculated correctly
			I was shifting Muon Energy same amount with Muon Momentum
			Now I shift Muon Momentum by dP and recalculate Shifted Muon Energy

		FillLatErrorBand_SingleUniverse(MnvH2D*)
			No longer updating BinError manually
			Original MnvLatErrorBand2D::Fill() function do not update errors
			Similarly MnvVertErrorBand2D::Fill() function do not update errors

	Using FillHistogram_ByHand() on TruthAnalysis also
		Want to make sure I am doing exactly same thing on both cases
		According to my previous tests, FillHistogram_ByHand()  and FillHistogram() are same
			But I just wanted to be safe, than sorry
		Note: BckgConstraints are not applied on TruthAnalysis, because histograms filled only for signal events

	W Binning Modified
		W is becoming an interesting parameter in my results, started optimizing it

	Plotter Systematics Functions are Improved
		Cleaned up all Practice Functions
		Plot and Table functions are improved for more informative output

" .

cvs tag -F ${CCPROTONPI0_V} .

