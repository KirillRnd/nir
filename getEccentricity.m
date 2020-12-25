function e = getEccentricity(r,V,mug)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
h=cross(r,V);
e_vec=cross(V,h)/mug-r/norm(r);
e=norm(e_vec);
end

