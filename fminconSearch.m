%clearvars -except symF
%clc;
%5if exist('symF','var') ~= 1
%    symbolic_Jacob
%end
t_start = juliandate(2001,12,1);
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

r_norm=ae;
V_unit=sqrt(mug_0/ae);
T_unit = T_earth/(2*pi);

[r0, V0] = planetEphemeris(t_start,'SolarSystem','Earth','430');
[rf, Vf] = planetEphemeris(t_start,'SolarSystem','Mars','430');

z0 = 

[t, y] = ode45(@(y)integrateTraectory(y), z0);

tic;
[tau, z] = ode45(@(y)integrateExternal(y, b), z0);
toc
y0 = cat(1, u0, v0, 0, tau0,  z')';

int_s0sf = linspace(0, s_f, (n+1)*1e+4);
%options = odeset('Events', @(s, y) eventIntegrationTraj(s, y, tf));
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);
%Интегрируем, используя сопряженные переменные из fmincon

[s,y] = ode113(@(s,y) integrateTraectory(s,y,h0),int_s0sf, y0, options);
Jt = integrateFunctional(s, y, eta, h0);
functional = Jt(end);

figure(2);
plot(t/(24*3600), vecnorm(a, 2, 2)*1e+03, 'LineWidth', 3);
%title('Зависимость ускорения силы тяги от времени')
xlabel('t, время, дни','FontSize',14)
ylabel('Реактивное ускорение, мм/c^2','FontSize',14)
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
th = 0:pi/50:2*pi;

t_orbit = linspace(t_start,t_start+T_earth/(24*3600), 1000);
earth_traj = planetEphemeris(t_orbit','SolarSystem','Earth','430', 'AU');

t_orbit = linspace(t_start,t_start+T_mars/(24*3600), 1000);
mars_traj = planetEphemeris(t_orbit','SolarSystem','Mars','430', 'AU');

%plot(cos(th),sin(th),'k');
%plot(1.52*cos(th),1.52*sin(th),'r');
plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k')
plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'r')

mars_r_f=planetEphemeris([t_start, t_end/(24*3600)],'SolarSystem','Mars','430','AU');

plot3(rr(:, 1), rr(:, 2), rr(:, 3), 'b', 'LineWidth', 2.5);
a_scale=3e-01/mean(vecnorm(a, 2, 2));
%a_scale=0;
d = 24*3600;
idxes=1;
for i=1:ceil(t(end)/d)
    ix = find(t>d*i*10, 1);
    idxes=[idxes, ix];
end    
for i = idxes
    plot3([rr(i, 1), rr(i, 1)+a_scale*a(i, 1)], [rr(i, 2), rr(i, 2)+a_scale*a(i, 2)],[rr(i, 3), rr(i, 3)+a_scale*a(i, 3)],'k')
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

[mars_r_f, mars_v_f]=planetEphemeris([t_start, t_end/(24*3600)],'SolarSystem','Mars','430');

disp(['Расход массы ', num2str(m(1)-m(end)), 'кг'])
disp(['Невязка координаты ', num2str(norm(rr(end, 1:3)*r_norm-mars_r_f*1e+03),'%10.2e\n'),',м'])
disp(['Невязка скорости ', num2str(norm(VV(end, 1:3)*V_unit-mars_v_f*1e+03),'%10.2e\n'),',м/с'])
% относительное число обусловленности
disp(['Относительное число обусловленности ', num2str(norm(x)*norm(grad)/fval,'%10.2e\n')])
% абсолютное число обусловленности
%disp(['Абсолютное число обусловленности ', num2str(1/norm(grad),'%10.2e\n')])