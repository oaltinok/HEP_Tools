CCPROTONPI0_V="v2_90"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Implemented and Tested FillVertErrorBand() and FillLatErrorBand() for MnvH2D

	Revised EM Uncertainty Correction
		Corrected Energy is calculated from Corrected Momentum
		Corrected Kinetic Energy is calculated from Corrected Energy

	Fixed a Bug in DrawDataMC_WithRatio() for XSec Variables
		I added a line to Normalize Histograms for another task.
		Normalization is already done by CrossSection class
		Removed 2nd Normalization from Plot function

	CrossSection Calculations Tested
		Added functionality to do calculations without systematics
			Remove error bands
" .

cvs tag -F ${CCPROTONPI0_V} .

