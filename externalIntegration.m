function res = externalIntegration(tau, z, b,dUdr,ddUdrdr,jac_ddUdrdr,y0,tspan)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mug = 132712.43994*(10^6)*(10^(3*3));
ae = 149597870700;
%y0=[1 0 0 0 -mug^(1/2) 0 0 0 0 0 0 0]';
y0_z=y0;
y0_z(7:12)=z;
T=2*pi*sqrt((1*ae)^3/mug);

options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);

[t,y] = ode45(@(t,y) internalIntegration(t,y,dUdr,ddUdrdr,jac_ddUdrdr),tspan,y0_z,options);
%plot(y(:,1),y(:,2));
%axis equal
drdpv=reshape(y(end,13:21),[3,3]);
dvdpv=reshape(y(end,22:30),[3,3]);

drdz=cat(2,drdpv,dvdpv);

ddrdpvdt=reshape(y(end,49:57),[3,3]);
ddVdpvdt=reshape(y(end,58:66),[3,3]);
ddrdzdt=cat(2,ddrdpvdt,ddVdpvdt);

dfdz = cat(1,drdz,ddrdzdt);

res=-(dfdz^-1)*b';
tau
end

