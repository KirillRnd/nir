function res = internalIntegration(t,y,dUdr,ddUdrdr,jac_ddUdrdr,mu_tau,tau)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
mu0=mu_tau(0);
mu_t=mu_tau(tau);
mug=mu_tau(1);
res=zeros(84,1);

r=y(1:3);
V=y(4:6);
pv=y(7:9);
dpvdt=y(10:12);%dpv/dt

drdpr=reshape(y(13:21),[3,3]);
drdpv=reshape(y(22:30),[3,3]);

drdz=cat(2,drdpr,drdpv);
ddrdzdt=y(31:48);

dpvdpr=reshape(y(49:57),[3,3]);
dpvdpv=reshape(y(58:66),[3,3]);
dpvdz=cat(2,dpvdpr,dpvdpv);

ddpvdzdt=y(67:84);

res(1:3)=V;
res(4:6)=(mu_t/mug)*dUdr(r)+pv;
res(7:9)=dpvdt;
res(10:12)=(mu_t/mug)*ddUdrdr(r)*pv;
res(13:30)=ddrdzdt;
res(49:66)=ddpvdzdt;

dddrdzdtdt=(mu_t/mug)*ddUdrdr(r)*drdz+dpvdz;
dddrdprdtdt=dddrdzdtdt(1:3,1:3);
dddrdpvdtdt=dddrdzdtdt(1:3,4:6);
res(31:39)=reshape(dddrdprdtdt,[9,1]);
res(40:48)=reshape(dddrdpvdtdt,[9,1]);


dddpvdzdtdt=(mu_t/mug)*(jac_ddUdrdr(r,pv)*drdz+ddUdrdr(r)*dpvdz);
dddpvdprdtdt=dddpvdzdtdt(1:3,1:3);
dddpvdpvdtdt=dddpvdzdtdt(1:3,4:6);
res(67:75)=reshape(dddpvdprdtdt,[9,1]);
res(76:84)=reshape(dddpvdpvdtdt,[9,1]);

drdtau=y(85:87);
ddrdtaudt=y(88:90);
dpvdtau=y(91:93);
ddpvdtaudt=y(94:96);
dddrdtaudtdt=(1-mu0/mug)*dUdr(r)+(mu_t/mug)*ddUdrdr(r)*drdtau+dpvdtau;
dddpvdtaudtdt=(1-mu0/mug)*ddUdrdr(r)*pv+(mu_t/mug)*(jac_ddUdrdr(r,pv)*drdtau+ddUdrdr(r)*dpvdtau);
res(85:87)=ddrdtaudt;
res(88:90)=dddrdtaudtdt;
res(91:93)=ddpvdtaudt;
res(94:96)=dddpvdtaudtdt;
end

