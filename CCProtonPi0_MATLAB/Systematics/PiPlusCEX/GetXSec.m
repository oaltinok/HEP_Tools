function [XSec] = GetXSec(n_cex)

mb_unit = 1.0*10^27;
flux = 1.0*10^7;  
N_Target = N_CarbonAtoms(); % per cm3

XSec = n_cex/(flux*N_Target)*mb_unit;

end



