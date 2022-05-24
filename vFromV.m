function w = vFromV(V, r, mug, phi)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
u = rToU(r, phi);
L = L_KS(u); 
h = (norm(V)^2)/2-mug/norm(r);
dtds=u'*u/sqrt(-2*h);
dxds=V*dtds;
w = L'*dxds/(2*u'*u);
end

