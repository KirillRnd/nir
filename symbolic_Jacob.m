mug = 132712.43994*(10^6)*(10^(3*3));

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

a=L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v));

duds=v;
dhds=2*v'*L'*a;
dvds=-u/4-(u2)/(4*h)*L'*a-(dhds/h)*v;
H=-a'*a/2+pu'*duds+pv'*dvds+ph'*dhds;

dpuds=-gradient(H, u');
dpvds=-gradient(H, v');
dphds=-gradient(H, h);
dtauds=(mug+4*u2*dhds+u2*u'*L'*a)/((-2*h)^(3/2));
dptauds=0;

y = [u', v', h, pu', pv', ph, tau, ptau];

f = [duds', dvds', dhds, dpuds', dpvds', dphds, dtauds, dptauds]';

J = jacobian(f, y);

symF = matlabFunction(f)
%symJ = matlabFunction(J)