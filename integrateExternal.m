function res = integrateExternal(tau,z, b, case_traj, tspan,r0,V0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
y0 = cat(2,r0,V0,z',...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(eye(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(eye(3),[1,9])...
    )';

[t, y] = ode45(@(t,y)integrateTraectory(t,y),tspan, y0);

drdz=reshape(y(end,13:30),[3,6]);
ddrdzdt=reshape(y(end,49:66),[3,6]);

dpVdz=reshape(y(end,31:48),[3,6]);
if case_traj == 1
    f = cat(1,drdz, dpVdz);
elseif  case_traj == 2
    f = cat(1,drdz, ddrdzdt);
end
res=-f^(-1)*b';
end

