function res = integrateExternal(z, b, case_traj)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[t, y] = ode45(@(y)integrateTraectory(y), z);

drdz=reshape(y(end,13:30),[3,6]);
ddrdzdt=reshape(y(end,49:66),[3,6]);

dpVdz=reshape(y(end,31:48),[3,6]);
if case_traj == 1
    f = [drdz, ddrdzdt];
elseif  case_traj == 2
    f = [drdz, dpVdz];
end
res=-f^(-1)*b;
end

