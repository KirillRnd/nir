mug=1;
%Для Шперлинга-Бюрде
r = sym('r', [1 3],'real')';
Vs = sym('Vs', [1 3],'real')';
A = sym('A', [1 3],'real')';

%сопряжённые переменные

pX = sym('pX', [1 3],'real')';
pVX = sym('pVX', [1 3],'real')';
h = sym('h', 'real')';

%h=-mug*norm(r)/(r'*r+Vs'*Vs);

dtds=norm(r)/sqrt(-2*h);

%оптимальное управление в координатах Шперлинга-Бюрде
a = (pVX*(r'*r) + (Vs*Vs')*pVX)/dtds/(-2*h);
%уравнения движения
drds=Vs;
dhds=Vs'*a;
dVsds=-r+(-A+(r'*r)*a+(Vs'*a)*Vs)/(-2*h);
dtauds=(mug+2*(r'*Vs)/norm(r)*dhds+norm(r)*(r'*a))/(-2*h)^(3/2);
dAds=2*(Vs'*a)*r-(r'*a)*Vs-(Vs'*r)*a;

H=-dtds*(a'*a)/2+pX'*drds+pVX'*dVsds;

dpXdt=-gradient(H, r');
dpVXdt=-gradient(H, Vs');
%Правая часть для интегрирования

f = [drds;dVsds;dpXdt;dpVXdt;dtauds;dAds;dhds];
symFShperlingBurde = matlabFunction(f,'File','symFShperlingBurde','Optimize',true, 'Vars', {r,Vs,pX,pVX,A,h});

h_rv=-mug*norm(r)/(r'*r+Vs'*Vs);
dtds=norm(r)/sqrt(-2*h_rv);
F = [r;Vs/dtds];

pr = sym('pr', [1 3],'real')';
pV = sym('pV', [1 3],'real')';
prpV2pXpVX=[pr',pV']*jacobian(F,[r;Vs]);

prpV2pXpVXShperlingBurde = matlabFunction(prpV2pXpVX,'File','prpV2pXpVXShperlingBurde','Optimize',true, 'Vars', {r,Vs,pr,pV});
%Должно быть равно pV
a_pV = simplify(subs(a, [pX; pVX], prpV2pXpVX'));
a_pV = simplify(subs(a_pV, h, h_rv))