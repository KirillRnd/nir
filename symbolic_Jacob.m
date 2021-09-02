%mug = 132712.43994*(10^6)*(10^(3*3));
mug=1;
u = sym('u', [1 4],'real')';
w = sym('w', [1 4],'real')';
%h = sym('h','real');
h=-mug/(u'*u+4*(w'*w));
tau = sym('tau','real');

pu = sym('pu', [1 4],'real')';
pw = sym('pw', [1 4],'real')';
%ph = sym('ph','real');
%билинейное соотношение
bil=@(x)[x(4); -x(3); x(2); -x(1)];

u2 = u'*u;
L = L_KS(u);
%dLvdu = jacobian(L*w, u');
%dLudu = jacobian(L*u, u');

dtds=u2/sqrt(-2*h);
%a = sym('a', [1 4],'real')';
%Данная запись управления совпадает с находимой дальше
Lambda=L*(-(u2)*pw/(4*h) - w*(1/h)*(pw'*w));
m=Lambda(4)/u2;
LambdaTilde=Lambda-m*L*bil(u);
LambdaTilde=simplify(LambdaTilde);
a=LambdaTilde/dtds;

g=bil(u)'*L'*a;
%a(4)=0;
duds=w;
dhds=2*w'*(L'*a);
dwds=-u/4-(L'*a)*(u2)/(4*h)-(dhds/(2*h))*w;
dtauds=(mug+4*(u'*w)*dhds+u2*u'*(L'*a))/((-2*h)^(3/2));

H=-dtds*(a'*a)/2+pu'*duds+pw'*dwds-m*g;
H_opt=H;
%H_opt=LambdaTilde'*LambdaTilde/(dtds*2)+pu'*w-pw'*u/4;

%a_solved = solve(gradient(H,a)==0,a);
%a_S=simplify([a_solved.a1;a_solved.a2;a_solved.a3;a_solved.a4;]);
%a_S(4)=0;
matlabFunction(a,'File','a_reactive','Optimize',true, 'Vars', {u,w,pu,pw});
matlabFunction(H_opt,'File','calculateHamiltonian','Optimize',true, 'Vars', {u,w,pu,pw});
%duds=v;
%dhds=2*v'*L'*a_S;
%dvds=-u/4-(u2)/(4*h)*L'*a_S-(dhds/(2*h))*v;
%dtauds=(mug+4*(u'*v)*dhds+u2*u'*L'*a_S)/((-2*h)^(3/2));

%H=-dtds*(a_S'*a_S)/2+pu'*duds+pv'*dvds+ph'*dhds+ptau'*dtauds;


dpuds=-gradient(H_opt, u');
dpwds=-gradient(H_opt, w');
%dphds=-simplify(gradient(H_opt, h));

y = [u', w', pu', pw', tau];

f = [duds', dwds', dpuds', dpwds',dtauds]';
symF = matlabFunction(f,'File','symF','Optimize',true, 'Vars', {u,w,pu,pw});

%J = jacobian(f, y);
%symJ = matlabFunction(J,'File','symJ','Optimize',true, 'Vars', {u,w,pu,pw});

%y = 2*(u(1)*u(2)-u(3)*u(4));
%z = 2*(u(1)*u(3)+u(2)*u(4));
%g = [u(1)^2+u(4)^2; u(2)^2+u(3)^2; atan2(u(4),u(1))+atan2(u(3),u(2))];
%dgdu=jacobian(g,u);
%symFdgdu = matlabFunction(dgdu,'File','get_dgdu','Optimize',true, 'Vars', {u});
%%
% V = sym('V', [1 4],'real')';
% V(4)=0;
% F = [[V(1) V(2) V(3) V(4)];
%     [V(2) -V(1) -V(4) V(3)];
%     [V(3) V(4) -V(1) -V(2)];
%     [-V(4) V(3) -V(2) V(1)]];
% Fv = F^-1*w;
Lt=L(1:3,:);
F=[Lt*u;2*sqrt(-2*h)*Lt*w/(u'*u)];
%g=[g(1:3);g(5:7)];

pr = sym('pr', [1 3],'real')';
pV = sym('pV', [1 3],'real')';

dgduw=jacobian(F,[u;w]);
ortdgduw=null(dgduw);
prpV2pupw=[pr;pV]'*dgduw;

matlabFunction(F,'File','get_target_g','Optimize', true, 'Vars', {u,w});
matlabFunction(dgduw,'File','get_dgduv','Optimize', true, 'Vars', {u,w});
matlabFunction(ortdgduw,'File','get_ortdgduv','Optimize', true, 'Vars', {u,w});
matlabFunction(prpV2pupw,'File','prpV2pupw','Optimize', true, 'Vars', {u,w,pr,pV});

F=F(1:3);

dgdu=jacobian(F,u);
ortdgdu=null(dgdu);

matlabFunction(F,'File','get_target_g_u','Optimize', true, 'Vars', {u});
matlabFunction(dgdu,'File','get_dgdu','Optimize', true, 'Vars', {u});
matlabFunction(ortdgdu,'File','get_ortdgdu','Optimize', true, 'Vars', {u});
%%
a_pV = simplify(subs(a, [pu;pw], prpV2pupw'))