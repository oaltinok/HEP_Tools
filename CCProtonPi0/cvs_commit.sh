CCPROTONPI0_V="v2_96"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:

---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Fixed ProtonEnergy Lateral Error Band Fill
		I was filling only for Events with a Proton
		Actually, I need to fill them for events without a Proton as well.

	New Class: BckgConstrainer
		Uses a map to hold Background Weights for all universes
" .

cvs tag -F ${CCPROTONPI0_V} .

