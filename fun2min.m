function dis = fun2min(x, case_traj, t_start, r0, V0)
%UNTITLED Summary of this function goes here
% Функция расстояния до Марса, в квадратах координаты-скорости.
% Зависит от сопряжённых переменных в начальный момент времени

mug_0 = 132712.43994*(10^6)*(10^(3*3));
ae = 149597870700;
V_unit=sqrt(mug_0/ae);
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;

T_unit = T_earth/(2*pi);

mug=1;

pu0=x(1:4)';
pv0=x(5:8)';
ph0=x(9);
pt0=x(10);
s_f=x(11);



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

%options = odeset('Events', @(s, y) eventIntegrationTraj(s, y,  t_f));
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);
options = odeset(options,'NonNegative', 10);

[s,y] = ode113(@(s,y) integrateTraectory(s, y),[0 s_f],y0, options);

u=y(end, 1:4)';
v=y(end, 5:8)';
h_end=y(end, 9)';
ph=y(end, 19)';
tau=y(end, 10)';
ptau=y(end, 20)';
t_end = T_unit*(tau-2*(u'*v)/sqrt(-2*h_end))/(24*60*60);
r_end=KS(u);
L_end = L_KS(u);
V_end = 2*sqrt(-2*h_end)*L_end*v/(norm(u)^2);

[rf, Vf] = planetEphemeris(t_end+t_start,'SolarSystem','Mars','430');

rf = [rf/ae, 0]'*1e+03;
Vf = [Vf/V_unit, 0]'*1e+03;
hf=(norm(Vf)^2)/2-mug/norm(rf);

uf(4) = 0;
uf(1) = sqrt((norm(r0)+r0(1))/2);
uf(2) = r0(2)/(2*uf(1));
uf(3) = r0(3)/(2*uf(1));

L = L_KS(uf); 
vf = L'*Vf/(2*sqrt(-2*hf));

%ЗАДАЧА ПРОЛЁТА
if case_traj == 1
    pv=y(end, 15:18);
    dis = norm((rf-r_end))^2 + (norm(pv)^2)*1e+5;
elseif case_traj == 2
    dis_p = [rf-r_end; Vf-V_end; ptau;];
    dis = norm(dis_p)^2;
end
end

