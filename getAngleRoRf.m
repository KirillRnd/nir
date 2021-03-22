function phi = getAngleRoRf(r0,rf,V0)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
h=cross(r0,V0);
n=h/norm(h);
rf_1 = rf-n*(rf'*n);
phi = acos((r0'*rf_1)/(norm(rf_1)*norm(r0)));
if cross(r0,rf_1)'*n < 0
    phi = 2*pi -phi;
end
end

