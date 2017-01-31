%% Breit Wigner
amp .* (x/0.938).^3 .* ((mass.*gamma) ./ ((x.^2 - mass.^2).^2 + (mass.*gamma).^2))

%% Data Fit
deltaRES1 .* (x/0.938).^3 .* ((deltaRES2.*deltaRES3) ./ ((x.^2 - deltaRES2.^2).^2 + (deltaRES2.*deltaRES3).^2)) + otherRES1*exp(-((x-otherRES2)/otherRES3)^2) + nonRES1*exp(-((x-nonRES2)/nonRES3)^2) + nonRES4*exp(-((x-nonRES5)/nonRES6)^2) 

deltaRES1 .* (x/0.938).^3 .* ((1.215.*deltaRES3) ./ ((x.^2 - 1.215.^2).^2 + (1.215.*deltaRES3).^2)) + otherRES1*exp(-((x-1.473)/otherRES3)^2) + nonRESamp .* (28.66*exp(-((x-1.197)/0.1498)^2) + 55.57*exp(-((x-1.439)/0.3686)^2)) 
