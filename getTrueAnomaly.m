function theta = getTrueAnomaly(r,V,mug)
%UNTITLED Summary of this function goes here
h=cross(r,V);
e=cross(V,h)/mug-r/norm(r);
if norm(e) == 0
    e = [1 0 0];
end

theta_cos=e'*r/(norm(e)*norm(r));
theta_sin=norm(cross(e, r))/norm(e)/norm(r);
theta = atan2(theta_cos, theta_sin);
if theta < 0
    theta=2*pi+theta;
end
end

