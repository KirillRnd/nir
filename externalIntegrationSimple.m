function res = externalIntegrationSimple(tau,z,b,y0,tspan,case_traj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
y0_z=y0;
y0_z(7:12)=z;

options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);

[t,y] = ode113(@(t,y)internalIntegrationSimple(t,y),tspan,y0_z,options);
drdpr=reshape(y(end,13:21),[3,3]);
drdpv=reshape(y(end,22:30),[3,3]);
drdz=cat(2,drdpr,drdpv);

ddrdprdt=reshape(y(end,31:39),[3,3]);
ddrdpvdt=reshape(y(end,40:48),[3,3]);
ddrdzdt=cat(2,ddrdprdt,ddrdpvdt);

dpvdpr=reshape(y(49:57),[3,3]);
dpvdpv=reshape(y(58:66),[3,3]);
dpvdz=cat(2,dpvdpr,dpvdpv);

if case_traj == 1
    dfdz = cat(1,drdz,dpvdz);
elseif case_traj == 2
    dfdz = cat(1,drdz,ddrdzdt);
end
res=-dfdz\b;
end

