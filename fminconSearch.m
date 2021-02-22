%clearvars -except symF
%clc;
%5if exist('symF','var') ~= 1
%    symbolic_Jacob
%end
%t_start = juliandate(2001,0,0);
t_start=0;

N=1350;
m0=367;
eta=0.45;
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
mug_0 = 132712.43994*(10^6)*(10^(3*3));
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;
T_mars_days = 365.256363004*1.8808476;



r_norm=ae;
V_unit=sqrt(mug_0/ae);
T_unit = T_earth/(2*pi);
planet_start = 'Earth';
planet_end = 'Mars';
%[r0, V0] = planetEphemeris(t_start,'SolarSystem',planet_start,'430');

r0 = [1, 0, 0, 0]';
V0 = [0, 1, 0, 0]';

mug=1;

n=1;
angle=0.5;
rad=1/32;

%t_f = T_earth*(n + angle);

%d_mars=-0.25;
%n_M = floor(t_f/T_mars);
%angle_M = t_f/T_mars-n_M;
%t_Mars_0 = (d_mars+angle-angle_M);

t_Mars_0=-0.5;

modifier_p=1e-04;
modifier_f=1e+04;
phi = n + angle;

s_a = phi - rad;
s_b = phi + rad;

x0(11) = phi;

lb = -[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]*1e+13;
ub = -lb;

lb(11) = s_a;
ub(11) = s_b;
%домножаем на коэффициент 1е-12, чтобы fmincon работал с более крупными
%величинами и не выдавал лишних ворнингов
tic;
fun=@(x)fun2min([x(1:10)*modifier_p x(11)], case_traj, t_start, r0, V0, planet_end, t_Mars_0*T_mars_days, modifier_f);

options = optimoptions('fmincon','UseParallel', true);
options = optimoptions(options, 'Display', 'iter');
options = optimoptions(options, 'OptimalityTolerance', 1e-10);
options = optimoptions(options, 'MaxFunctionEvaluations', 1e+10);
options = optimoptions(options, 'StepTolerance', 1e-12);
options = optimoptions(options, 'ConstraintTolerance', 1e-10);

%options = optimoptions(options, 'Algorithm', 'sqp');
%options = optimoptions(options, 'HessianApproximation','lbfgs');
%options = optimoptions(options, 'ScaleProblem', 'obj-and-constr');
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub,[], options);
toc
px = x(1:10)*modifier_p;
s_f = x(11)*2*pi;
%задаем начальные условия

t0=0;
u0 = [0 0 0 0]';
h0=(norm(V0)^2)/2-mug/norm(r0);

u0(4) = 0;
u0(1) = sqrt((norm(r0)+r0(1))/2);
u0(2) = r0(2)/(2*u0(1));
u0(3) = r0(3)/(2*u0(1));

L = L_KS(u0); 
v0 = L'*V0/(2*sqrt(-2*h0));
%tau0= getEccentricAnomaly(r0(1:3),V0(1:3),mug);
tau0=0;
y0 = cat(1, u0, v0, 0, tau0,  px')';

int_s0sf = linspace(0, s_f, (n+1)*1e+4);
%options = odeset('Events', @(s, y) eventIntegrationTraj(s, y, tf));
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);
%Интегрируем, используя сопряженные переменные из fmincon

[s,y] = ode113(@(s,y) integrateTraectory(s,y,h0),int_s0sf, y0, options);
Jt = integrateFunctional(s, y, eta, h0);
functional = Jt(end);

uu = y(:, 1:4);
vv = y(:, 5:8);
rr=zeros(length(uu),4);
a=zeros(length(uu),4);
a_ks=zeros(length(uu),4);
t=zeros(length(uu),1);
VV=zeros(length(uu),4);
for i = 1:length(uu)
    u = uu(i,:)';
    r=KS(u);
    rr(i,:)=r;
    L=L_KS(u);
    u2=norm(u)^2;
    v=y(i, 5:8)';
    h=y(i, 9)'+h0;
    tau=y(i ,10)';
    pu=y(i, 11:14)';
    pv=y(i, 15:18)';
    ph=y(i, 19)';
    ptau=y(i, 20)';
    dtds=u2/sqrt(-2*h);
    aa_ks=L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*(u2)*u/((-2*h)^(3/2)))/dtds;
    a_ks(i, :)=aa_ks/(ae/sqrt(mug_0)).^2;
    res=symF(u,v,h,pu,pv,ph,ptau);
    dvds=res(5:8);
    dhds=res(9);
    V = 2*sqrt(-2*h)*L*v/(u2);
    VV(i, :)=V;
    %a(i, :)=((-2*h/(norm(r)^2))*(2*(L_KS(v)*v+L_KS(u)*dvds)-(2*u'*v/(sqrt(-2*h)) + norm(r)*dhds/((-2*h)^(3/2)))*V)+mug*r/(norm(r)^3))/(ae/sqrt(mug_0)).^2;
    
    %a(i, :)=KS(aa);
    t(i) = T_unit*(tau-2*(u'*v)/sqrt(-2*h));
end
t = t - t(1);
t_end=t(end);

figure(2);
plot(t/(24*3600), vecnorm(a_ks, 2, 2)*1e+03, 'LineWidth', 3);
%title('Зависимость ускорения силы тяги от времени')
xlabel('t, время, дни','FontSize',14)
ylabel('Реактивное ускорение, мм/c^2','FontSize',14)
box off;
set(gca,'FontSize',14)

figure(3);
plot(t/(24*3600), s, 'LineWidth', 3);
title('Зависимость мнимого времени от физического')
xlabel('t, время, дни')
ylabel('s, мнимое время')
box off;
set(gca,'FontSize',14)

figure(4);
m=massLP(Jt, m0, N);
plot(t/(24*3600), m, 'LineWidth', 3);
title('Зависимость массы от времени')
xlabel('t, время, дни')
ylabel('m, масса, кг')
box off;
set(gca,'FontSize',14)

%Проверка "на глаз"
figure(1);
plot3(0, 0, 0, 'y--o')
set(gca,'FontSize',14)
hold on;

th = linspace(0 ,2*pi,100)';
mars_traj = 1.52*[cos(th), sin(th), zeros(100,1)];
earth_traj  = [cos(th), sin(th), zeros(100,1)];
%plot(cos(th),sin(th),'k');
%plot(1.52*cos(th),1.52*sin(th),'r');
plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k')
plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'r')

n_M = floor((t_end+t_Mars_0*T_mars)/T_mars);
angle_M = ((t_end+t_Mars_0*T_mars)/T_mars-n_M)*2*pi;

mars_r_f = 1.52*[cos(angle_M) sin(angle_M) 0 0]';
mars_v_f = ((mug/(1.52))^(1/2))*[cos(angle_M+pi/2) sin(angle_M+pi/2) 0 0]';
%mars_r_f=planetEphemeris([t_start, t_end/(24*3600)],'SolarSystem',planet_end,'430','AU');

plot3(rr(:, 1), rr(:, 2), rr(:, 3), 'b', 'LineWidth', 2.5);
a_scale=3e-01/mean(vecnorm(a_ks, 2, 2));
%a_scale=0;
d = 24*3600;
idxes=1;
for i=1:ceil(t(end)/d)
    ix = find(t>d*i*10, 1);
    idxes=[idxes, ix];
end    
for i = idxes
    plot3([rr(i, 1), rr(i, 1)+a_scale*a_ks(i, 1)], [rr(i, 2), rr(i, 2)+a_scale*a_ks(i, 2)],[rr(i, 3), rr(i, 3)+a_scale*a_ks(i, 3)],'k')
end

plot3(rr(end, 1), rr(end, 2), rr(end, 3),'bO')
plot3(mars_r_f(1), mars_r_f(2),mars_r_f(3),'rO')
axis equal

%title('Траектория КА')
xlabel('x, a.e.')
ylabel('y, a.e.')

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
hold off;

%Проверка "на глаз"
figure(5);
plot3(0, 0, 0, 'y--o')
set(gca,'FontSize',14)
hold on;

th = linspace(0 ,4*pi,100)';
mars_traj = 1.52*[cos(th), sin(th), zeros(100,1)];
earth_traj  = [cos(th), sin(th), zeros(100,1)];

mars_traj_ks = arrayfun(@(r1, r2, r3, th) rToU([r1,r2,r3],th), mars_traj(:, 1),mars_traj(:, 2),mars_traj(:, 3), th,'UniformOutput',false);
mars_traj_ks = cell2mat(mars_traj_ks')';
earth_traj_ks = arrayfun(@(r1, r2, r3, th) rToU([r1,r2,r3],th), earth_traj(:, 1),earth_traj(:, 2),earth_traj(:, 3), th,'UniformOutput',false);
earth_traj_ks = cell2mat(earth_traj_ks')';
plot3(earth_traj_ks(:, 1), earth_traj_ks(:, 2), earth_traj_ks(:, 3), 'k')
plot3(mars_traj_ks(:, 1), mars_traj_ks(:, 2), mars_traj_ks(:, 3), 'r')

plot3(uu(:, 1), uu(:, 2), uu(:, 3), 'b', 'LineWidth', 2.5);
%a_scale=3e-01/mean(vecnorm(a_ks, 2, 2));
a_scale=0;
d = 24*3600;
idxes=1;
for i=1:ceil(t(end)/d)
    ix = find(t>d*i*10, 1);
    idxes=[idxes, ix];
end    
for i = idxes
    %plot3([uu(i, 1), uu(i, 1)+a_scale*a_ks(i, 1)], [uu(i, 2), uu(i, 2)+a_scale*a_ks(i, 2)],[uu(i, 3), uu(i, 3)+a_scale*a_ks(i, 3)],'k')
end

%plot3(uu(end, 1), uu(end, 2), uu(end, 3),'bO')
%plot3(mars_r_f(1), mars_r_f(2),mars_r_f(3),'rO')
axis equal

title('Траектория КА KS')
xlabel('x, a.e.')
ylabel('y, a.e.')

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
hold off;
%[mars_r_f, mars_v_f]=planetEphemeris([t_start, t_end/(24*3600)],'SolarSystem',planet_end,'430');

disp(['Расход массы ', num2str(m(1)-m(end)), 'кг'])
disp(['Невязка координаты ', num2str(ae*(norm(rr(end, 1:3)-mars_r_f(1:3)')),'%10.2e\n'),',м'])
disp(['Невязка скорости ', num2str(V_unit*(norm(VV(end, 1:3)-mars_v_f(1:3)')),'%10.2e\n'),',м/с'])
% относительное число обусловленности
disp(['Относительно е число обусловленности ', num2str(norm(x)*norm(grad)/fval,'%10.2e\n')])
% абсолютное число обусловленности
%disp(['Абсолютное число обусловленности ', num2str(1/norm(grad),'%10.2e\n')])