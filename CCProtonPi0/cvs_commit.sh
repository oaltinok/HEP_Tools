CCPROTONPI0_V="v2_84"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	CrossSection Improved
		New Background Estimate Method
			In new method, dont use TFractionFitter just use MC shape
			MC shape already fitted in SideBand Studies
			Background Estimate using TFractionFitter top of Side Band Constraints is optional
" .

cvs tag -F ${CCPROTONPI0_V} .

