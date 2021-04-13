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
%билинейное соотношение
bil=@(x)[x(4); -x(3); x(2); -x(1)];

u2 = u'*u;
L = [[u(1) -u(2) -u(3) u(4)];
    [u(2) u(1) -u(4) -u(3)];
    [u(3) u(4) u(1) u(2)];
    [u(4) -u(3) u(2) -u(1)]]; 
dLvdu = jacobian(L*v, u');
dLudu = jacobian(L*u, u');

dtds=u2/sqrt(-2*h);
a = sym('a', [1 4],'real')';
%a=L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*(u2)*u/((-2*h)^(3/2)))/dtds;

duds=v;
dhds=2*v'*L'*a;
dvds=-u/4-(u2)/(4*h)*L'*a-(dhds/(2*h))*v;
dtauds=(mug+4*(u'*v)*dhds+u2*u'*L'*a)/((-2*h)^(3/2));

H=-dtds*(a'*a)/2+pu'*duds+pv'*dvds+ph'*dhds+ptau'*dtauds;

a_solved = solve(gradient(H,a)==0,a);
a_S=simplify([a_solved.a1;a_solved.a2;a_solved.a3;a_solved.a4;]);
matlabFunction(a_S,'File','a_reactive','Optimize',true, 'Vars', {u,v,h,pu,pv,ph,ptau});

duds=v;
dhds=2*v'*L'*a_S;
dvds=-u/4-(u2)/(4*h)*L'*a_S-(dhds/(2*h))*v;
dtauds=(mug+4*(u'*v)*dhds+u2*u'*L'*a_S)/((-2*h)^(3/2));

H=-dtds*(a_S'*a_S)/2+pu'*duds+pv'*dvds+ph'*dhds+ptau'*dtauds;


dpuds=-simplify(gradient(H, u'));
dpvds=-simplify(gradient(H, v'));
dphds=-simplify(gradient(H, h));
dptauds=-simplify(gradient(H, tau));

y = [u', v', h, tau, pu', pv', ph, ptau];

f = [duds', dvds', dhds, dtauds, dpuds', dpvds', dphds,  dptauds]';

J = jacobian(f, y);

symF = matlabFunction(f,'File','symF','Optimize',true, 'Vars', {u,v,h,pu,pv,ph,ptau});
symJ = matlabFunction(J,'File','symJ','Optimize',true, 'Vars', {u,v,h,pu,pv,ph,ptau});

y = 2*(u(1)*u(2)-u(3)*u(4));
z = 2*(u(1)*u(3)+u(2)*u(4));
g = [u(1)^2+u(4)^2; u(2)^2+u(3)^2; atan2(z,y)];
dgdu=jacobian(g,u);
symFdgdu = matlabFunction(dgdu,'File','get_dgdu','Optimize',true, 'Vars', {u});

V = sym('V', [1 4],'real')';
V(4)=0;
F = [[V(1) V(2) V(3) V(4)];
    [V(2) -V(1) -V(4) V(3)];
    [V(3) V(4) -V(1) -V(2)];
    [-V(4) V(3) -V(2) V(1)]];
Fv = F^-1*v;
gv(1)=v'*v;
%gv(1)=simplify(Fv(1)^2+Fv(4)^2)*((V'*V)^2);
%gv(2)=simplify(Fv(2)^2+Fv(3)^2)*((V'*V)^2);
gv(2)=simplify(Fv(1)^2-Fv(2)^2-Fv(3)^2+Fv(4)^2);
gv(3)=atan2(simplify(Fv(4)*(V'*V)), simplify(Fv(1)*(V'*V)))+...
    atan2(simplify(Fv(3)*(V'*V)), simplify(Fv(2)*(V'*V)));

dgdv=jacobian(gv,v);
%не получится записать V(4)=0 
V = sym('V', [1 3],'real')';
symgv = matlabFunction(gv,'File','get_gv','Optimize', true, 'Vars', {v,V});
symFdgdv = matlabFunction(dgdv,'File','get_dgdv','Optimize', true, 'Vars', {v,V});


%Вторая формула для параметрической скорости
L123=L(1:3,:);
GR=L123*L123';
v_proj_coef=GR\L123*v;
v_proj=L123'*v_proj_coef;
v_ort=simplify(v-v_proj);
v_ort_norm=simplify(norm(v_ort)^2);
v_b=[v(4); -v(3); v(2); -v(1)];
u_v=simplify(F^-1*v);