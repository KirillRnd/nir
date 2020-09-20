function [Jt, dis, s, y] = trajectorySearch(n, angle, rad, case_traj, symF, eta)
%trajectorySearch Summary of this function goes here
%Обёртка над fmincon, также вычисляем расстояние в виде суммы квпадратов

%Начальные условия
x0=[0 0 0 0 0 0 0 0 0 0 0];
A = [];
b = [];
Aeq = [];
beq = [];
ae = 149597870700;
mug = 132712.43994*(10^6)*(10^(3*3));
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;

modifier=1e-8;
modifier_p=1e-15; 

koef = 2;

s_a = (n + angle-rad)*pi*koef;
s_b = (n + angle+rad)*pi*koef;

x0(11)=(n + angle)*pi*koef;

t_f = T_earth*(n + angle);

n_M = floor(t_f/T_mars);
angle_M = t_f/T_mars-n_M;
t_Mars_0 = (angle-angle_M)*T_mars;
lb = -[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]*1e+05;
ub = -lb;

lb(11) = s_a;
ub(11) = s_b;
%домножаем на коэффициент 1е-12, чтобы fmincon работал с более крупными
%величинами и не выдавал лишних ворнингов

fun=@(x)fun2min([x(1:10)*modifier_p x(11)], case_traj, symF, t_Mars_0);
x = fmincon(fun, x0, A, b, Aeq, beq, lb, ub);
px = x(1:10)*modifier_p;
s_f = x(11);
%задаем начальные условия
r0 = [1*ae 0 0 0]';
V0 = [0 (mug/(1*ae))^(1/2) 0 0]';

t0=0;
u0 = [0 0 0 0]';
h0=(norm(V0)^2)/2-mug/norm(r0);

u0(4) = 0;
u0(1) = sqrt((norm(r0)+r0(1))/2);
u0(2) = r0(2)/(2*u0(1));
u0(3) = r0(3)/(2*u0(1));

L = L_KS(u0); 
v0 = L'*V0/(2*sqrt(-2*h0));
tau0=0;
y0 = cat(1, u0, v0, h0, tau0,  px')';

int_s0sf = linspace(0, s_f, (n+1)*1e+4);
%options = odeset('Events', @(s, y) eventIntegrationTraj(s, y, tf));
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);
%Интегрируем, используя сопряженные переменные из fmincon

[s,y] = ode113(@(s,y) integrateTraectory(s,y,symF),int_s0sf, y0, options);
Jt = integrateFunctional(s, y, symF, eta);
u=y(end, 1:4)';
v=y(end, 5:8)';
h=y(end, 9)';
tau=y(end, 10)';
t_end = tau-2*(u'*v)/(-2*h);
r_end=KS(u)';

n_M = floor((t_end+t_Mars_0)/T_mars);
angle_M = ((t_end+t_Mars_0)/T_mars-n_M)*2*pi;

rf = 1.52*ae*[cos(angle_M) sin(angle_M) 0 0];
Vf = ((mug/(1.52*ae))^(1/2))*[cos(angle_M+pi/2) sin(angle_M+pi/2) 0 0];

%Сумма квадратов невязок для задачи пролёта
if case_traj == 1
    pv=y(end, 15:18);
    dis = norm((rf-r_end)/norm(r0))^2 + (norm(pv)^2)*1e+20;
elseif case_traj == 2
    v = y(end, 5:8)';
    h = y(end, 9);
    Lend = L_KS(u);
    V = 2*sqrt(-2*h)*Lend*v/(norm(u)^2);
    dis = norm((rf-r_end)/norm(r0))^2 + (norm((V'-Vf)/norm(V0))^2);
end
end

