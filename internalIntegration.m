function res = internalIntegration(t,y,dUdr,ddUdrdr,jac_ddUdrdr)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
res=zeros(84,1);

r=y(1:3);
V=y(4:6);
pv=y(7:9);
Pv=y(10:12);

drdpv=reshape(y(13:21),[3,3]);
dvdpv=reshape(y(22:30),[3,3]);

drdz=cat(2,drdpv,dvdpv);

dpvdpv=reshape(y(31:39),[3,3]);
dPvdpv=reshape(y(40:48),[3,3]);
dpvdz=cat(2,dpvdpv,dPvdpv);

drdzdt=y(49:66);
dpvdzdt=y(67:84);

res(1:3)=V;
res(4:6)=dUdr(r)+pv;
res(7:9)=Pv;
res(10:12)=ddUdrdr(r)*pv;
res(13:30)=drdzdt;
res(31:48)=dpvdzdt;

ddrdzdt=ddUdrdr(r)*drdz+dpvdz;
ddrdpvdt=ddrdzdt(1:3,1:3);
ddVdpvdt=ddrdzdt(1:3,4:6);
res(49:57)=reshape(ddrdpvdt,[9,1]);
res(58:66)=reshape(ddVdpvdt,[9,1]);


dPVdpvdt=jac_ddUdrdr(r,pv)*drdz+ddUdrdr(r)*dpvdz;
dpvdpvdt=dPVdpvdt(1:3,1:3);
dPvdpvdt=dPVdpvdt(1:3,4:6);
res(67:75)=reshape(dpvdpvdt,[9,1]);
res(76:84)=reshape(dPvdpvdt,[9,1]);
end

