function E = getEccentricAnomaly(r,V,mug)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
theta = getTrueAnomaly(r,V,mug);
e = getEccentricity(r,V,mug);
E=2*atan(sqrt((1-e)/(1+e))*tan(theta/2));
if E < 0
    E = 2*pi+E;
end
end

