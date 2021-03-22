function e = getEccentricity(r,V,mug)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
c=cross(r,V);
f=cross(V,c)-mug*r/norm(r);
e=norm(f)/mug;
end

