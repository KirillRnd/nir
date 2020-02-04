function [v_opt,t_opt] = zad_exp(gamma)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here



mu_Earth = 6.672*10^13*6;
mu_sun = 2*10^19*6.672;
v_Earth = 29861.7;
r0 = [380000000;0;0];
v0 = [0;sqrt(mu_Earth/norm(r0));0];
pv_0 = [0;1;0];
pv_dot_0 = [0;0;0];
u_dv = 29430*v0/norm(v0); %скорость истечения

m0 = 100000;
t = 0;
ro_deistv = 924217000;
mug = mu_Earth;
y = [r0;v0;pv_0;pv_dot_0;0];
flag = 1;
for t = 0:1:86400*124
    r = y(1:3);
    v = y(4:6);
    pv = y(7:9);
    pv_dot = y(10:12);
    m_rash = y(13);
    if (flag)
        u_dv = 29430*v/norm(v);
    else
        u_dv = 29430*pv/norm(pv);
    end
    if (norm(r) >= ro_deistv && flag)
        mug = mu_sun;
        r = r+149597870700*r/norm(r);
        y(1:3) = r;
        v = v+v/norm(v)*v_Earth;        
        y(4:6) = v;
        flag = 0;
        v_opt = v;
        t_opt = t;
        break;
    end
    dy = [v;-mug/norm(r)^3*r+gamma*u_dv/(m0 - m_rash);...
        pv_dot;mug/norm(r)^3*(3*r'*pv*r/norm(r)^2-pv);gamma];
    y = y + dy*1;
end

end

