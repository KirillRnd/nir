function res = internalIntegrationSimple(t,y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
mug=1;
res=zeros(84,1);

r=y(1:3);
V=y(4:6);
pv=y(7:9);
dpvdt=y(10:12);%dpv/dt

%pv=y(7:9);
%dpvdt=y(10:12);%dpv/dt

drdpr=reshape(y(13:21),[3,3]);
drdpv=reshape(y(22:30),[3,3]);

drdz=cat(2,drdpr,drdpv);
ddrdzdt=y(31:48);

dpvdpr=reshape(y(49:57),[3,3]);
dpvdpv=reshape(y(58:66),[3,3]);
dpvdz=cat(2,dpvdpr,dpvdpv);

ddpvdzdt=y(67:84);

res(1:3)=V;
res(4:6)=mug*sym_dUdr(r)+pv;%РИТЭГ
%res(4:6)=(mu_t/mug)*dUdr(r)+pv/(norm(r)^2);
res(7:9)=dpvdt;
res(10:12)=mug*sym_ddUdrdr(r)*pv;%РИТЭГ
%res(10:12)=(mu_t/mug)*ddUdrdr(r)*pv-r*norm(pv)^2/norm(r)^4;
res(13:30)=ddrdzdt;
res(49:66)=ddpvdzdt;

dddrdzdtdt=mug*sym_ddUdrdr(r)*drdz+dpvdz;
dddrdprdtdt=dddrdzdtdt(1:3,1:3);
dddrdpvdtdt=dddrdzdtdt(1:3,4:6);
res(31:39)=reshape(dddrdprdtdt,[9,1]);
res(40:48)=reshape(dddrdpvdtdt,[9,1]);


dddpvdzdtdt=mug*(jac_sym_ddUdrdr(r,pv)*drdz+sym_ddUdrdr(r)*dpvdz);
dddpvdprdtdt=dddpvdzdtdt(1:3,1:3);
dddpvdpvdtdt=dddpvdzdtdt(1:3,4:6);
res(67:75)=reshape(dddpvdprdtdt,[9,1]);
res(76:84)=reshape(dddpvdpvdtdt,[9,1]);

end

