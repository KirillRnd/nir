function theta = getTrueAnomaly(r,V,mug)
%UNTITLED Summary of this function goes here
h=cross(r,V);
e=cross(V,h)/mug-r/norm(r);
theta=acos(e*r'/(norm(e)*norm(r)));
if r'*V < 0
    theta=2*pi-theta;
end
end

