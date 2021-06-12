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
%���������� �����������
bil=@(x)[x(4); -x(3); x(2); -x(1)];

u2 = u'*u;
L = [[u(1) -u(2) -u(3) u(4)];
    [u(2) u(1) -u(4) -u(3)];
    [u(3) u(4) u(1) u(2)];
    [u(4) -u(3) u(2) -u(1)]]; 
dLvdu = jacobian(L*v, u');
dLudu = jacobian(L*u, u');

dtds=u2/sqrt(-2*h);
%a = sym('a', [1 4],'real')';
%������ ������ ���������� ��������� � ��������� ������
Lambda=L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*((u2)*u+8*(u'*v)*v)/((-2*h)^(3/2)));
m=Lambda(4)/u2;
a=(Lambda-m*L*bil(u))/dtds;
a=simplify(a);
g=bil(u)'*L'*a;
%a(4)=0;
duds=v;
dhds=2*v'*L'*a;
dvds=-u/4-(u2)/(4*h)*L'*a-(dhds/(2*h))*v;
dtauds=(mug+4*(u'*v)*dhds+u2*u'*L'*a)/((-2*h)^(3/2));

H=-dtds*(a'*a)/2+pu'*duds+pv'*dvds+ph'*dhds+ptau'*dtauds-m*g;

%a_solved = solve(gradient(H,a)==0,a);
%a_S=simplify([a_solved.a1;a_solved.a2;a_solved.a3;a_solved.a4;]);
%a_S(4)=0;
matlabFunction(a,'File','a_reactive','Optimize',true, 'Vars', {u,v,h,pu,pv,ph,ptau});

%duds=v;
%dhds=2*v'*L'*a_S;
%dvds=-u/4-(u2)/(4*h)*L'*a_S-(dhds/(2*h))*v;
%dtauds=(mug+4*(u'*v)*dhds+u2*u'*L'*a_S)/((-2*h)^(3/2));

%H=-dtds*(a_S'*a_S)/2+pu'*duds+pv'*dvds+ph'*dhds+ptau'*dtauds;


dpuds=-simplify(gradient(H, u'));
dpvds=-simplify(gradient(H, v'));
dphds=-simplify(gradient(H, h));
dptauds=-simplify(gradient(H, tau));

y = [u', v', h, tau, pu', pv', ph, ptau];

f = [duds', dvds', dhds, dtauds, dpuds', dpvds', dphds,  dptauds]';
%symF = matlabFunction(f,'File','symF','Optimize',true, 'Vars', {u,v,h,pu,pv,ph,ptau});

J = jacobian(f, y);
symJ = matlabFunction(J,'File','symJ','Optimize',false, 'Vars', {u,v,h,pu,pv,ph,ptau});

y = 2*(u(1)*u(2)-u(3)*u(4));
z = 2*(u(1)*u(3)+u(2)*u(4));
g = [u(1)^2+u(4)^2; u(2)^2+u(3)^2; atan2(u(4),u(1))+atan2(u(3),u(2))];
dgdu=jacobian(g,u);
symFdgdu = matlabFunction(dgdu,'File','get_dgdu','Optimize',true, 'Vars', {u});

V = sym('V', [1 4],'real')';
V(4)=0;
F = [[V(1) V(2) V(3) V(4)];
    [V(2) -V(1) -V(4) V(3)];
    [V(3) V(4) -V(1) -V(2)];
    [-V(4) V(3) -V(2) V(1)]];
Fv = F^-1*v;
g=[L*u;2*sqrt(-2*h)*L*v/(u'*u);h];
g=[g(1:3);g(5:7);g(9)];

dgduvh=jacobian(g,[u;v;h]);
ortdgduvh=null(dgduvh);


matlabFunction(g,'File','get_target_gh','Optimize', true, 'Vars', {u,v,h});
matlabFunction(dgduvh,'File','get_dgduvh','Optimize', true, 'Vars', {u,v,h});
matlabFunction(ortdgduvh,'File','get_ortdgduvh','Optimize', true, 'Vars', {u,v,h});

g=g(1:3);

dgdu=jacobian(g,u);
ortdgdu=null(dgdu);


matlabFunction(g,'File','get_target_g_u','Optimize', true, 'Vars', {u});
matlabFunction(dgdu,'File','get_dgdu','Optimize', true, 'Vars', {u});
matlabFunction(ortdgdu,'File','get_ortdgdu','Optimize', true, 'Vars', {u});
