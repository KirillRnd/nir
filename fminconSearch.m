clear;
clc;
symbolic_Jacob


%условия на fmincon
%ЗАДАЧА ПРОЛЁТА case_traj=1; ЗАДАЧА сопровождения case_traj=2;
case_traj=1;
%Количество витков
%n = 4;
%angle = 6*pi/6;
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

n=3;
angle=0.5;
rad=0.01;

modifier=1e-8;
modifier_p=1e-15; 
tf_a = T_earth*(n + angle-rad)*modifier;
tf_b = T_earth*(n + angle+rad)*modifier;

x0(11)=T_earth*(n + angle)*modifier;

n_M = floor((x0(11)/modifier)/T_mars);
angle_M = (x0(11)/modifier)/T_mars-n_M;
t_Mars_0 = (angle-angle_M-0.03)*T_mars;
lb = -[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]*10000;
ub = -lb;

lb(11) = tf_a;
ub(11) = tf_b;
%домножаем на коэффициент 1е-12, чтобы fmincon работал с более крупными
%величинами и не выдавал лишних ворнингов

fun=@(x)fun2min([x(1:10)*modifier_p x(11)/modifier], case_traj, symF, t_Mars_0);
x = fmincon(fun, x0, A, b, Aeq, beq, lb, ub);
px = x(1:10)*modifier_p;
tf = x(11)/modifier;
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

L = [[u0(1) -u0(2) -u0(3) u0(4)];
    [u0(2) u0(1) -u0(4) -u0(3)];
    [u0(3) u0(4) u0(1) u0(2)];
    [u0(4) -u0(3) u0(2) -u0(1)]]; 
v0 = L'*V0/(2*sqrt(-2*h0));
tau0=0;
y0 = cat(1, u0, v0, h0, tau0,  px')';


n = floor(tf/T_earth);
angle = (tf/T_earth-n)*2*pi;

n_M = floor((tf+t_Mars_0)/T_mars);
angle_M = ((tf+t_Mars_0)/T_mars-n_M)*2*pi;

sf = (n*2*pi+angle)*1.5;
int_s0sf = linspace(0, sf, (n+1)*1e+4);
options = odeset('Events', @(s, y) eventIntegrationTraj(s, y, tf));
options = odeset(options,'AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);
%Интегрируем, используя сопряженные переменные из fmincon

[s,y] = ode113(@(s,y) integrateTraectory(s,y,symF),int_s0sf,y0, options);
functional = integrateFunctional(s, y);
functional = functional(end);

u = y(:, 1:4);
r=zeros(length(u),4);
a=zeros(length(u),4);
t=zeros(length(u),1);
for i = 1:length(u)
    rr = u(i,:)';
    L = [[rr(1) -rr(2) -rr(3) rr(4)];
    [rr(2) rr(1) -rr(4) -rr(3)];
    [rr(3) rr(4) rr(1) rr(2)];
    [rr(4) -rr(3) rr(2) -rr(1)]];
    r(i,:)=L*rr;
    u2=norm(u)^2;
    v=y(i, 5:8)';
    h=y(i, 9)';
    tau=y(i ,10)';
    pu=y(i, 11:14)';
    pv=y(i, 15:18)';
    ph=y(i, 19)';
    ptau=y(i, 20)';
    aa=L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*(rr'*rr)*rr/(-2*h)^(3/2));
    
    La = [[aa(1) -aa(2) -aa(3) aa(4)];
    [aa(2) aa(1) -aa(4) -aa(3)];
    [aa(3) aa(4) aa(1) aa(2)];
    [aa(4) -aa(3) aa(2) -aa(1)]];
    a(i, :)=La*aa;
    t(i) = tau-2*(rr'*v)/(-2*h);
end



figure(2);
plot(t, vecnorm(a, 2, 2));
title('Зависимость ускорения силы тяги от времени')
xlabel('t, время')
ylabel('a, силы тяги')

figure(3);
plot(t, s);
title('Зависимость мнимого времени от обычного')
xlabel('t, время')
ylabel('s, мнимое время')

%Проверка "на глаз"
figure(1);
plot(0, 0,'y--o')
hold on;
th = 0:pi/50:2*pi;
plot(ae*cos(th),ae*sin(th),'k');
plot(1.52*ae*cos(th),1.52*ae*sin(th),'r');
plot(r(:, 1), r(:, 2),'b')
a_scale=1e+10;

d = 24*3600;
idxes=1;
for i=1:ceil(t(end)/d)
    ix = find(t>d*i*10, 1);
    idxes=[idxes, ix];
end    
for i = idxes
    plot([r(i, 1), r(i, 1)+a_scale*a(i, 1)], [r(i, 2), r(i, 2)+a_scale*a(i, 2)],'k')
end
plot(r(end, 1), r(end, 2),'bO')
plot(1.52*ae*cos(angle_M), 1.52*ae*sin(angle_M),'rO')
axis equal

title('Траектория КА')
xlabel('x')
ylabel('y')

hold off;
