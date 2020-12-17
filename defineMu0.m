function mu0 = defineMu0(N,r0,V0,rf,tft0,mug)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
M0=getMeanAnomaly(r0,V0,mug);

theta0=getTrueAnomaly(r0,V0,mug);
thetaf=theta0+getAngleRoRf(r0,rf,V0);

e = getEccentricity(r0,V0,mug);
Ef=2*atan(sqrt((1-e)/(1+e))*tan(thetaf/2));
Mf=Ef-e*sin(Ef);
h=norm(cross(r0,V0));
a0=(h^2)/(mug*(1-e^2));
mu0=((Mf+2*pi*N+M0)/tft0)^2*a0^3;
end

