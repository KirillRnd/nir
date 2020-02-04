function [v_opt,t_opt] = zad_exp(gamma)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

mug = 398600.4415*(10^(3*3));
v_Earth = 29861.7;
r0 = [380*1e+06;0;0];
v0 = [0;sqrt(mug/norm(r0));0];
ro_deistv = 924217000;
y = [r0;v0;0];
h = 1;
for i = 0:1:86400*124
    
    k1 = firstStageRungeK(y,gamma,h*i);
    k2 = firstStageRungeK(y+h/2*k1, gamma, h*(i+1/2));
    k3 = firstStageRungeK(y+h/2*k2, gamma, h*(i+1/2));
    k4 = firstStageRungeK(y+h*k3, gamma, h*(i+1));

    Y=h/6*(k1+2*k2+2*k3+k4);
    y=y+Y;
    r = y(1:3);
    v = y(4:6);
    if (norm(r) >= ro_deistv)
        v = v+v/norm(v)*v_Earth;     
        v_opt = v;
        t_opt = i*h;
        break;
    end
end

end
