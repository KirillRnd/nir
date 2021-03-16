function v = vFromV(V,r,mug)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
u = rToU(r);
L = L_KS(u); 
h = (norm(V)^2)/2-mug/norm(r);
v = L'*V/(2*sqrt(-2*h));
end

