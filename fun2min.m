function dis = fun2min(x, case_traj, t_Mars_0, amax)
%UNTITLED Summary of this function goes here
% Функция расстояния до Марса, в квадратах координаты-скорости.
% Зависит от сопряжённых переменных в начальный момент времени

mug = 132712.43994*(10^6)*(10^(3*3));
ae = 149597870700;
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;

pu0=x(1:4)';
pv0=x(5:8)';
ph0=x(9);
pt0=x(10);
s_f=x(11);

r0 = [1*ae 0 0 0]';
V0 = [0 (mug/(1*ae))^(1/2) 0 0]';

u0 = [0 0 0 0]';
h0 = (norm(V0)^2)/2-mug/norm(r0);

u0(4) = 0;
u0(1) = sqrt((norm(r0)+r0(1))/2);
u0(2) = r0(2)/(2*u0(1));
u0(3) = r0(3)/(2*u0(1));

L = L_KS(u0); 

v0 = L'*V0/(2*sqrt(-2*h0));
t0 = 0;

y0 = cat(1, u0, v0, h0, t0, pu0, pv0, ph0, pt0)';

%
%StartTime=clock;
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);
%options = odeset(options, 'Events', @(s, y) eventIntegrationTraj(s, y, StartTime));

%options = odeset(options,'MinStep',1e-10); 
[s,y] = ode113(@(s,y) integrateTraectory(s,y, amax),[0 s_f],y0,options);

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

%ЗАДАЧА ПРОЛЁТА
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

