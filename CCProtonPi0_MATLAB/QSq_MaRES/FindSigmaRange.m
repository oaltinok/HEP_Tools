function [MaRES_Roots] = FindSigmaRange( sigma )

p1 = 65.8;
p2 = -188;
p3 = 128.1;
p4 = 8.735;

p = [p1 p2 p3 (p4-sigma.^2)];

MaRES_Roots = roots(p);

end

