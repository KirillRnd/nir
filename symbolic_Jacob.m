%mug = 132712.43994*(10^6)*(10^(3*3));
mug=1;
u = sym('u', [1 4],'real')';
w = sym('w', [1 4],'real')';
h = sym('h','real');
tau = sym('tau','real');
amax = sym('amax','real')';
pu = sym('pu', [1 4],'real')';
pw = sym('pw', [1 4],'real')';
ph = sym('ph','real');
ptau = sym('ptau','real');
%билинейное соотношение
bil=@(x)[x(4); -x(3); x(2); -x(1)];

u2 = u'*u;
L = [[u(1) -u(2) -u(3) u(4)];
    [u(2) u(1) -u(4) -u(3)];
    [u(3) u(4) u(1) u(2)];
    [u(4) -u(3) u(2) -u(1)]]; 
dLvdu = jacobian(L*w, u');
dLudu = jacobian(L*u, u');

dtds=u2/sqrt(-2*h);
%a = sym('a', [1 4],'real')';
%Данная запись управления совпадает с находимой дальше
Lambda=L*(-(u2)*pw/(4*h) + w*(2*ph-(1/h)*pw'*w)+ptau*((u2)*u+8*(u'*w)*w)/((-2*h)^(3/2)));
m=Lambda(4)/u2;
LambdaTilde=Lambda-m*L*bil(u);
LambdaTilde=simplify(LambdaTilde);
a=amax*LambdaTilde/norm(LambdaTilde);

g=bil(u)'*L'*a;
%a(4)=0;
duds=w;
dhds=2*w'*L'*a;
dvds=-u/4-(u2)/(4*h)*L'*a-(dhds/(2*h))*w;
dtauds=(mug+4*(u'*w)*dhds+u2*u'*L'*a)/((-2*h)^(3/2));

%H=-dtds*(a'*a)/2+pu'*duds+pw'*dvds+ph'*dhds+ptau'*dtauds-m*g;

H_opt=-dtds + LambdaTilde'*a+pu'*w-pw'*u/4+ptau*mug/((-2*h)^(3/2));
%a_solved = solve(gradient(H,a)==0,a);
%a_S=simplify([a_solved.a1;a_solved.a2;a_solved.a3;a_solved.a4;]);
%a_S(4)=0;
matlabFunction(a,'File','a_reactive','Optimize',true, 'Vars', {u,w,h,pu,pw,ph,ptau,amax});

%duds=v;
%dhds=2*v'*L'*a_S;
%dvds=-u/4-(u2)/(4*h)*L'*a_S-(dhds/(2*h))*v;
%dtauds=(mug+4*(u'*v)*dhds+u2*u'*L'*a_S)/((-2*h)^(3/2));

%H=-dtds*(a_S'*a_S)/2+pu'*duds+pv'*dvds+ph'*dhds+ptau'*dtauds;


dpuds=-simplify(gradient(H_opt, u'));
dpvds=-simplify(gradient(H_opt, w'));
dphds=-simplify(gradient(H_opt, h));
dptauds=-simplify(gradient(H_opt, tau));

y = [u', w', h, tau, pu', pw', ph, ptau];

f = [duds', dvds', dhds, dtauds, dpuds', dpvds', dphds,  dptauds]';
symF = matlabFunction(f,'File','symF','Optimize',true, 'Vars', {u,w,h,pu,pw,ph,ptau,amax});

J = jacobian(f, y);
symJ = matlabFunction(J,'File','symJ','Optimize',true, 'Vars', {u,w,h,pu,pw,ph,ptau,amax});

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
Fv = F^-1*w;
g=[L*u;2*sqrt(-2*h)*L*w/(u'*u);h];
g=[g(1:3);g(5:7);g(9)];

dgduvh=jacobian(g,[u;w;h]);
ortdgduvh=null(dgduvh);


matlabFunction(g,'File','get_target_gh','Optimize', true, 'Vars', {u,w,h});
matlabFunction(dgduvh,'File','get_dgduvh','Optimize', true, 'Vars', {u,w,h});
matlabFunction(ortdgduvh,'File','get_ortdgduvh','Optimize', true, 'Vars', {u,w,h});

g=g(1:3);

dgdu=jacobian(g,u);
ortdgdu=null(dgdu);


matlabFunction(g,'File','get_target_g_u','Optimize', true, 'Vars', {u});
matlabFunction(dgdu,'File','get_dgdu','Optimize', true, 'Vars', {u});
matlabFunction(ortdgdu,'File','get_ortdgdu','Optimize', true, 'Vars', {u});
