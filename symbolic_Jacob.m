%mug = 132712.43994*(10^6)*(10^(3*3));
mug=1;
u = sym('u', [1 4],'real')';
v = sym('v', [1 4],'real')';
h = sym('h','real');
tau = sym('tau','real');

pu = sym('pu', [1 4],'real')';
pv = sym('pv', [1 4],'real')';
ph = sym('ph','real');
ptau = sym('ptau','real');

u2 = u'*u;
L = [[u(1) -u(2) -u(3) u(4)];
    [u(2) u(1) -u(4) -u(3)];
    [u(3) u(4) u(1) u(2)];
    [u(4) -u(3) u(2) -u(1)]]; 
dLvdu = jacobian(L*v, u');

a=L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*(u2)*u/((-2*h)^(3/2)));

duds=v;
dhds=2*v'*L'*a;
dvds=-u/4-(u2)/(4*h)*L'*a-(dhds/(2*h))*v;
dtauds=(mug+4*(u'*v)*dhds+u2*u'*L'*a)/((-2*h)^(3/2));

H=-a'*a/2+pu'*duds+pv'*dvds+ph'*dhds+ptau'*dtauds;

dpuds=-gradient(H, u');
dpvds=-gradient(H, v');
dphds=-gradient(H, h);
dptauds=-gradient(H, tau);

%y = [u', v', h, tau, pu', pv', ph, ptau];

f = [duds', dvds', dhds, dtauds, dpuds', dpvds', dphds,  dptauds]';

%J = jacobian(f, y);

symF = matlabFunction(f,'File','symF','Optimize',true, 'Vars', {u,v,h,pu,pv,ph,ptau});
%symJ = matlabFunction(J)