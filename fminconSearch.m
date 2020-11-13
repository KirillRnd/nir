N=1350;
m0=367;
eta=0.45;
amax=2.4e-4;

%условия на fmincon
%ЗАДАЧА ПРОЛЁТА case_traj=1; ЗАДАЧА сопровождения case_traj=2;
case_traj=2;

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

n=0;

angle=0.5;
rad=0.3;

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
lb = -[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]*1e+03;
ub = -lb;

lb(11) = s_a;
ub(11) = s_b;
%домножаем на коэффициент 1е-12, чтобы fmincon работал с более крупными
%величинами и не выдавал лишних ворнингов

fun=@(x)fun2min([x(1:10)*modifier_p x(11)], case_traj, t_Mars_0, amax);

options = optimoptions('fmincon','UseParallel', true);
options = optimoptions(options, 'Display', 'iter');
%options = optimoptions(options, 'OptimalityTolerance', 1e-07);
%options = optimoptions(options, 'MaxIterations', 75);

x = fmincon(fun, x0, A, b, Aeq, beq, lb, ub,[], options)
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

[s,y] = ode113(@(s,y) integrateTraectory(s,y, amax),int_s0sf, y0, options);
Jt = integrateFunctional(s, y);
functional = Jt(end);

uu = y(:, 1:4);
rr=zeros(length(uu),4);
a=zeros(length(uu),4);
t=zeros(length(uu),1);
VV=zeros(length(uu),4);
for i = 1:length(uu)
    u = uu(i,:)';
    r=KS(u);
    rr(i,:)=r;
    L=L_KS(u);
    u2=norm(u)^2;
    v=y(i, 5:8)';
    h=y(i, 9)';
    tau=y(i ,10)';
    pu=y(i, 11:14)';
    pv=y(i, 15:18)';
    ph=y(i, 19)';
    ptau=y(i, 20)';
    %aa=L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*(rr'*rr)*rr/(-2*h)^(3/2));
    %res=symF(amax, h,ph,ptau,pu(1),pu(2),pu(3),pu(4),pv(1),pv(2),pv(3),pv(4),u(1),u(2),u(3),u(4),v(1),v(2),v(3),v(4));
    %dvds=res(5:8);
    %dhds=res(9);
    %V = 2*sqrt(-2*h)*L*v/(u2);
    %VV(i, :)=V;
    %lambda = L*(-pv*u2/(4*h)+(2*ph-pv'*v/h)*v+ptau*u*u2/((-2*h)^(3/2)));
    %a = amax*lambda/norm(lambda);%a(i, :)=KS(aa);
    t(i) = tau-2*(u'*v)/sqrt(-2*h);
end

t_end=t(end);

n = floor(t_end/T_earth);
angle = (t_end/T_earth-n)*2*pi;

n_M = floor((t_end+t_Mars_0)/T_mars);
angle_M = ((t_end+t_Mars_0)/T_mars-n_M)*2*pi;

figure(3);
plot(t/(24*3600), s);
title('Зависимость мнимого времени от обычного')
xlabel('t, время в днях')
ylabel('s, мнимое время')
figure(4);

m=massLP(t*(amax^2), m0, N);
plot(t/(24*3600), m);
title('Зависимость массы от времени')
xlabel('t, время в днях')
ylabel('m, масса, кг')

%Проверка "на глаз"
figure(1);
plot(0, 0,'y--o')
hold on;
th = 0:pi/50:2*pi;
plot(ae*cos(th),ae*sin(th),'k');
plot(1.52*ae*cos(th),1.52*ae*sin(th),'r');
plot(rr(:, 1), rr(:, 2),'b')
%a_scale=3e+10/mean(vecnorm(a, 2, 2));
a_scale=0;
d = 24*3600;
idxes=1;
for i=1:ceil(t(end)/d)
    ix = find(t>d*i*10, 1);
    idxes=[idxes, ix];
end    
for i = idxes
    plot([rr(i, 1), rr(i, 1)+a_scale*a(i, 1)], [rr(i, 2), rr(i, 2)+a_scale*a(i, 2)],'k')
end
plot(rr(end, 1), rr(end, 2),'bO')
plot(1.52*ae*cos(angle_M), 1.52*ae*sin(angle_M),'rO')
axis equal

title('Траектория КА')
xlabel('x, м')
ylabel('y, м')

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

hold off;