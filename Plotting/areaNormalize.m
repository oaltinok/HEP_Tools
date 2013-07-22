function [data, MC] = areaNormalize(data,MC)
%  data = areaNormalize(data,MC)
%   Author: Ozgur Altinok
%   Date: 2013-06-29
%   Area Normalizes the data and MC

%%
sumData = sum(data);
sumMC = sum(MC);

if (sumMC > sumData)
    ratio = sumMC / sumData;
    MC = MC;
    data = data .* ratio;
else
    ratio = sumData / sumMC;
    data  = data;
    MC = MC .* ratio;
end


end