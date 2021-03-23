function phi = getAngleRoRf(r0,rf,V0)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
c=cross(r0,V0);
n=c/norm(c);
rf_1 = rf-n*(rf'*n);
phi_cos = (r0'*rf_1)/(norm(rf_1)*norm(r0));
phi_sin = norm(cross(r0,rf_1))/(norm(rf_1)*norm(r0));
phi=atan2(phi_sin, phi_cos);
if cross(r0,rf_1)'*n < 0
    phi = 2*pi-phi;
end
end

