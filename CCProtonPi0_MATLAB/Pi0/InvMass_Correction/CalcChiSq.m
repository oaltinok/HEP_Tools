function [ ChiSq ] = CalcChiSq( measured, expected)

ChiSq = sum(((measured-expected).^2)./measured);

end

