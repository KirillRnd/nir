mug = 132712.43994*(10^6)*(10^(3*3));
%mug=1;
r = sym('r', [3 1],'real');
V = sym('V', [3 1],'real');

pV = sym('pV', [3 1],'real');
dpVdt = sym('dpVdt', [3 1],'real');

F=-mug*r/norm(r)^3;
dFdr=jacobian(F, r);

drdt=V;
dvdt=F+pV;

DpVdt=dpVdt;
ddpvdtdt=dFdr*pV;

drdz = sym('drdz', [3 6],'real');
ddrdzdt = sym('ddrdzdt', [3 6],'real');

dpVdz = sym('dpVdz', [3 6],'real');
ddpVdtdz = sym('ddpVdtdz', [3 6],'real');

Ddrdzdt=ddrdzdt;
dddrdzdtdt=dFdr*drdz+dpVdz;

DpVdz=ddpVdtdz;
dddpVdzdtdt=jacobian(dFdr*pV, r)*drdz+dFdr*dpVdz;


f =[drdt', dvdt', DpVdt', ddpvdtdt', reshape(Ddrdzdt,[18,1])', reshape(dddrdzdtdt,[18,1])', reshape(DpVdz,[18,1])', reshape(dddpVdzdtdt,[18,1])']';
symF = matlabFunction(f,'File','symF','Optimize',true, 'Vars', {r,V,pV,dpVdt,drdz,ddrdzdt,dpVdz,ddpVdtdz});
%symF = matlabFunction(f,'File','symF','Optimize',true, 'Vars', {u,v,h,pu,pv,ph,ptau});
%symJ = matlabFunction(J)