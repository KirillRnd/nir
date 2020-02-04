function f = firstStageRungeK(y, gamma, t)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
mug = 398600.4415*(10^(3*3));
m0 = 100000;

r = y(1:3);
v = y(4:6);
m_rash = y(7);
u_dv = 3000*9.81*v/norm(v);
f = [v;-mug/norm(r)^3*r+gamma*u_dv/(m0 - m_rash);gamma];
end

