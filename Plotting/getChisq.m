function chiSquare = getChisq(dataY,modelY,variance)
%  chiSquare = getChisq(dataY,modelY,variance)
%   Author: Ozgur Altinok
%   Date: 2013-06-07 

%%

chiSquareCont = ((dataY - modelY) ./ variance) .^ 2;
chiSquare = sum(chiSquareCont);

end