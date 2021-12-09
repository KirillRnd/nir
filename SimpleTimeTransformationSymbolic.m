mug=1;
%Для Шперлинга-Бюрде
rs = sym('rs', [1 3],'real')';
Vs = sym('Vs', [1 3],'real')';

%сопряжённые переменные

prs = sym('prs', [1 3],'real')';
pVs = sym('pVs', [1 3],'real')';

%h=-mug*norm(rs)/(rs'*rs+Vs'*Vs);
rs_norm=sqrt(rs'*rs);
dtds=rs_norm;
k=dtds;
%оптимальное управление в координатах rs vs
%a = pVX*dtds;
a = k*pVs;
%уравнения движения
drds=Vs;
dVsds=dtds^2*(a-mug*rs/rs_norm^3)+Vs*(rs'*Vs)/(rs_norm^2);

H=(k*(-a'*a)/2+(prs'*drds+pVs'*dVsds))*dtds^(-1);

dprsdt=-dtds^(1)*gradient(H, rs');
%Совпадает
%dpVsdt=-prs-(pVs'*Vs)*rs/(rs_norm^2)-(rs'*Vs)*pVs/(rs_norm^2);
dpVsdt=-dtds^(1)*gradient(H, Vs');
%Правая часть для интегрирования

f = [drds;dVsds;dprsdt;dpVsdt;dtds;];
symFSimpleTimeTransformation = matlabFunction(f,'File','symFSimpleTimeTransformation','Optimize',true, 'Vars', {rs,Vs,prs,pVs});

F = [rs;Vs/dtds];

pr = sym('pr', [1 3],'real')';
pV = sym('pV', [1 3],'real')';

r = sym('r', [1 3],'real')';
V = sym('V', [1 3],'real')';

prpV2pXpVX=[pr',pV']*jacobian(F,[rs;Vs]);

prpV2pXpVXSimple = matlabFunction(prpV2pXpVX,'File','prpV2pXpVXSimple','Optimize',true, 'Vars', {rs,Vs,pr,pV});
%Должно быть равно pV
a_pV = simplify(subs(a, [prs; pVs], prpV2pXpVX'))
