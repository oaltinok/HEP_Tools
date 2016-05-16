CCPROTONPI0_V="v2_80"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Added Muon Theta from All Nodes
		Theta estimate on vertex has a bias
		Phil showed that if we use a more downstream node the estimate is better
		Now saving theta estimate in all nodes
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Base Class: NTupleAnalysis
		Implemented GetMnvH1D function under Base Class:
			Uses dynamic_cast instead of c style casting
			Uses h->SetDirectory(NULL); to make sure there is no connection between TFile and MnvH1D
			returns a <new> MnvH1D
		Rewritten all MnvH1D and MnvH2D functions as template<class MnvHistoType>
		
	SideBandTool
		Removed Test Functions and cout statements used for debugging
		Moved and Improved GetMnvH1D to NTupleAnalysis

	CrossSection
		Using GetMnvH1D to make sure the SetDirectory(NULL) always called 
		Modified Flux Integration code to use FluxReweighter
    
    CutList
        Reduce mode no longer overwrites the ReducedNTuple.root and CutHistograms.root
                Protection against trying to reduce by mistake
                Now user needs to keep track of ROOT file names
		
	Plotter
		Removed separate DrawDataMC_CrossSection() functions
			Using a variable to control CrossSection plots on DrawDataMC()
		Added new function DrawDataMC_WithRatio()
			Plots Data-MC comparison and attaches the ratio below it

" .

cvs tag -F ${CCPROTONPI0_V} .

