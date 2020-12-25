function M = getMeanAnomaly(r,V,mug)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
E = getEccentricAnomaly(r,V,mug);
e = getEccentricity(r,V,mug);
M=E-e*sin(E);
end

