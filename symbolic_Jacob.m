mug = 132712.43994*(10^6)*(10^(3*3));

amax = sym('amax');

u = sym('u', [1 4])';
v = sym('v', [1 4])';
h = sym('h');
tau = sym('tau');

pu = sym('pu', [1 4])';
pv = sym('pv', [1 4])';
ph = sym('ph');
ptau = sym('ptau');

u2 = u'*u;
L = [[u(1) -u(2) -u(3) u(4)];
    [u(2) u(1) -u(4) -u(3)];
    [u(3) u(4) u(1) u(2)];
    [u(4) -u(3) u(2) -u(1)]]; 
dLvdu = jacobian(L*v, u');

lambda = L*(-pv*u2/(4*h)+(2*ph-pv'*v/h)*v+ptau*u*u2/((-2*h)^(3/2)));

a = amax*lambda/norm(lambda);

duds=v;
dhds=2*v'*L'*a;
dvds=-u/4-(u2)/(4*h)*L'*a-(dhds/(2*h))*v;
dtauds=(mug+4*(u'*v)*dhds+u2*u'*L'*a)/((-2*h)^(3/2));

H=-1+pu'*duds+pv'*dvds+ph'*dhds+ptau'*dtauds;

dpuds=-gradient(H, u');
dpvds=-gradient(H, v');
dphds=-gradient(H, h);
dptauds=-gradient(H, tau);

y = [u', v', h, tau, pu', pv', ph, ptau];

f = [duds', dvds', dhds, dtauds, dpuds', dpvds', dphds,  dptauds]';

%J = jacobian(f, y);

symF = matlabFunction(f);

a = [0 0 0 0]';

duds=v;
dhds=2*v'*L'*a;
dvds=-u/4-(u2)/(4*h)*L'*a-(dhds/(2*h))*v;
dtauds=(mug+4*(u'*v)*dhds+u2*u'*L'*a)/((-2*h)^(3/2));

H=-1+pu'*duds+pv'*dvds+ph'*dhds+ptau'*dtauds;

dpuds=-gradient(H, u');
dpvds=-gradient(H, v');
dphds=-gradient(H, h);
dptauds=-gradient(H, tau);

y = [u', v', h, tau, pu', pv', ph, ptau];

f = [duds', dvds', dhds, dtauds, dpuds', dpvds', dphds,  dptauds]';

symF_a0 = matlabFunction(f);
%symJ = matlabFunction(J)