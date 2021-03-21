function res = externalIntegration(tau,z,b,dUdr,ddUdrdr,jac_ddUdrdr,y0,tspan,mu_tau,V0,Vf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
y0_z=y0;
y0_z(4:6)=V0*sqrt(mu_tau(tau)/mu_tau(1));
y0_z(88:90)=0.5*sqrt(mu_tau(1)/mu_tau(tau))*(1-mu_tau(0)/mu_tau(1))*V0;
y0_z(7:12)=z;

options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);

[t,y] = ode113(@(t,y) internalIntegration(t,y,dUdr,ddUdrdr,jac_ddUdrdr,mu_tau,tau),tspan,y0_z,options);
%plot(y(:,1),y(:,2));
%axis equal
drdpv=reshape(y(end,13:21),[3,3]);
dvdpv=reshape(y(end,22:30),[3,3]);

drdz=cat(2,drdpv,dvdpv);

ddrdpvdt=reshape(y(end,49:57),[3,3]);
ddVdpvdt=reshape(y(end,58:66),[3,3]);
ddrdzdt=cat(2,ddrdpvdt,ddVdpvdt);

dfdz = cat(1,drdz,ddrdzdt);

drdtau=y(end,85:87);
dvdtau=y(end,88:90);
dfdtau = cat(2,drdtau,dvdtau-0.5*sqrt(mu_tau(1)/mu_tau(tau))*(1-mu_tau(0)/mu_tau(1))*Vf);
res=-dfdz\(dfdtau+b)';
tau
if tau == 1.0
    cond(dfdz)
end

