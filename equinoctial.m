%модифицированные равноденственные элементы, символьный скрипт
% Иванов Д.С., Трофимов С.П., Широбоков М.Г. Численное моделирование орбитального и
% углового движения космических аппаратов. / под общ. ред. М.Ю. Овчинникова. М.: ИПМ
% им.М.В.Келдыша, 2016. 118 с.
% doi:10.20948/mono-2016-trofimov
% URL: http://keldysh.ru/e-biblio/trofimov
%гравитационная постоянная
mug = sym('mug','real')';
%Определяем переменные
p = sym('p','real')';
ex = sym('ex','real')';
ey = sym('ey','real')';
ix = sym('ix','real')';
iy = sym('iy','real')';
L = sym('L','real')';

X = [p; ex; ey; ix;iy;L];
%ускорение
S = sym('S','real')'; %радиальный
T = sym('T','real')'; %трансверсальный
W = sym('W','real')'; %бинормальный

a = [S;T;W];

%вспомогательные выражения

phi = 1 + ix^2 + iy^2;
eta = ix*sin(L) - iy*cos(L);
sigma = 1 + ex*cos(L) + ey*sin(L);

%производные по времени
dpdt  = sqrt(p/mug)*(2*p/sigma)*T;
dexdt = sqrt(p/mug)*( sin(L)*S + ((1+1/sigma)*cos(L)+ex/sigma)*T + (-ey*eta/sigma)*W);
deydt = sqrt(p/mug)*(-cos(L)*S + ((1+1/sigma)*sin(L)+ey/sigma)*T +  (ex*eta/sigma)*W);
dixdt = sqrt(p/mug)*(phi*cos(L)/(2*sigma))*W;
diydt = sqrt(p/mug)*(phi*sin(L)/(2*sigma))*W;
dLdt  = sqrt(mug/p)*(sigma^2/p)+sqrt(p/mug)*(eta/sigma)*W;

dXdt = [dpdt; dexdt; deydt; dixdt;diydt;dLdt];
%Определяем сопряжённые переменные
Pp = sym('Pp','real')';
Pex = sym('Pex','real')';
Pey = sym('Pey','real')';
Pix = sym('Pix','real')';
Piy = sym('Piy','real')';
PL = sym('PL','real')';

P = [Pp; Pex; Pey; Pix;Piy;PL];
%Функция Гамильтона-Понтрягина
H = -a'*a/2 + P'*dXdt;

%оптимальное управление
Sopt = sqrt(p/mug)*(Pex*sin(L)-Pey*cos(L));
Topt_p  = Pp*2*p/sigma;
Topt_ex = Pex*((1+1/sigma)*cos(L)+ex/sigma);
Topt_ey = Pey*((1+1/sigma)*sin(L)+ey/sigma);
Topt = sqrt(p/mug)*(Topt_p+Topt_ex+Topt_ey);
Wopt = sqrt(p/mug)*((PL-Pex*ey+Pey*ex)*eta/sigma + (Pix*cos(L)+Piy*sin(L))*phi/(2*sigma));

aopt = [Sopt;Topt;Wopt];
Hopt = subs(H, a, aopt);
Hopt = simplify(Hopt);
dXdt_opt = subs(dXdt, a, aopt);
%производные сопряжённых переменных
dPdt=-gradient(Hopt, X');

X_extended = [X; P];
dXdt_extended = [dXdt_opt; dPdt];

dXdt_extended_mat = matlabFunction(dXdt_extended,'File','equitoctial_adjoint_for_integration','Optimize',true, 'Vars', {X_extended,mug});
%% с функционалом
a_norm = simplify(norm(a));
J_integrand = 0.5*(a_norm^2);
J_integrand_opt = subs(J_integrand, a, aopt);
dXdt_extended_J = [dXdt_opt; dPdt;J_integrand_opt];
dXdt_extended__J_mat = matlabFunction(dXdt_extended_J,'File','equitoctial_adjoint_for_integration_withJ','Optimize',true, 'Vars', {X_extended,mug});
%% переход из модифицированныех равноденственных в декартовые
% https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf

f = ex;
g = ey;
h = ix;
k = iy;
s2 = 1+h^2+k^2;
a2 = h^2-k^2;
w = 1+f*cos(L)+g*sin(L);
r = p/w;

rx = (r/s2)*(cos(L)+a2*cos(L)+2*h*k*sin(L));
ry = (r/s2)*(sin(L)-a2*sin(L)+2*h*k*cos(L));
rz = (2*r/s2)*(h*sin(L) - k*cos(L));

vx = -(sqrt(mug/p)/s2)*( sin(L)+a2*sin(L)-2*h*k*cos(L)+g-2*f*h*k+a2*g);
vy = -(sqrt(mug/p)/s2)*(-cos(L)+a2*cos(L)+2*h*k*sin(L)-f+2*g*h*k+a2*f);
vz = (2*sqrt(mug/p)/s2)*(h*cos(L)+k*sin(L)+f*h+g*k);

R = [rx; ry; rz];
V = [vx; vy; vz];

x = [R;V];
x_mat = matlabFunction(x,'File','equitoctial2decart','Optimize',true, 'Vars', {X,mug});
%% матрица частных производных
dxdX = jacobian(x, X');
dxdX_mat = matlabFunction(dxdX,'File','equitoctial2decart_jacobian','Optimize',true, 'Vars', {X,mug});