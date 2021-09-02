mug=1;
%цилиндрические координаты
rho = sym('rho','real')';
phi = sym('phi', 'real')';
z = sym('z','real')';

X = [rho; phi; z];
VX = sym('V', [1 3],'real')';

x = sym('x','real')';
y = sym('y','real')';
% обычные координаты
r = [x;y;z];



%сопряжённые переменные

pX = sym('pX', [1 3],'real')';
pVX = sym('pVX', [1 3],'real')';

y = [X;VX;pX;pVX];

%оптимальное управление в цилиндрических координатах
a = [pVX(1)*cos(phi)-pVX(2)*sin(phi)/rho; pVX(1)*sin(phi)+pVX(2)*cos(phi)/rho; pVX(3)];
%уравнения движения
dVrhodt=rho*(VX(2)^2)-mug*rho/((rho^2+z^2)^(3/2))+a(1)*cos(phi)+a(2)*sin(phi);
dVphidt=(-2*VX(1)*VX(2)-a(1)*sin(phi)+a(2)*cos(phi))/rho;
dVzdt=-mug*z/((rho^2+z^2)^(3/2))+a(3);

dVXdt=[dVrhodt;dVphidt;dVzdt];


H=-(a'*a)/2+pX'*VX+pVX'*dVXdt;

dpXdt=-simplify(gradient(H, X'));
dpVXdt=-simplify(gradient(H, VX'));
%Правая часть для интегрирования

f = [VX;dVXdt;dpXdt;dpVXdt];
symFCylindric = matlabFunction(f,'File','symFCylindric','Optimize',true, 'Vars', {X,VX,pX,pVX});

F = [rho*cos(phi); rho*sin(phi); z;...
    VX(1)*cos(phi)-rho*sin(phi)*VX(2); VX(1)*sin(phi)+rho*cos(phi)*VX(2); VX(3)];

pr = sym('pr', [1 3],'real')';
pV = sym('pV', [1 3],'real')';
prpV2pXpVX=[pr',pV']*jacobian(F,[X;VX]);

prpV2pXpVXCylindric = matlabFunction(prpV2pXpVX,'File','prpV2pXpVXCylindric','Optimize',true, 'Vars', {X,VX,pr,pV});
%Должно быть равно pV
simplify(subs(a, pVX, prpV2pXpVX(4:6)'))
