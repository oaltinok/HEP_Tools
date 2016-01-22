CCPROTONPI0_V="v2_52"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Shower Recovery Methods
		1Shower Recovery Methods
			See AngleScan Improvements for implementation details
			Use Smaller Cone Angle
			Use Different Search View
		3 Shower Recovery Methods
			Use Direction Cut 
				If 2 shower have GOOD direction, keep event

	AngleScan Improved
		New Mode:
			AllowSmallConeAngle
			If it is TRUE, bin size will be 2 degrees (default is 4 degrees)
		New Modes 
			AngleScan_U: Start search from U view
			AngleScan_V: Start search from V view
			Duplicated from original class
				I may integrate SearchView into original class, but for now they will stay different
				I do not have time to work on that for now
			Do not include Small Angle Mode
		Removed MaxDigit Energy Cleaning
			It is not effective, no need to do use it
		Removed zDistanceFromLessThan
			It is not used
				
-------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	No Major Changes
	Study for Shower Recovery Methods
" .

cvs tag -F ${CCPROTONPI0_V} .
