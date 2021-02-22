function dis = fun2min(x, case_traj, t_start, r0, V0, planet_end,t_Mars_0, modifier_f)
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
s_f=x(11)*2*pi;



u0 = rToU(r0,0);
h0 = (norm(V0)^2)/2-mug/norm(r0);

L = L_KS(u0); 

v0 = L'*V0/(2*sqrt(-2*h0));
%t0 = getEccentricAnomaly(r0(1:3),V0(1:3),mug);
t0=0;
y0 = cat(1, u0, v0, 0, t0, pu0, pv0, ph0, pt0)';

time0 = tic;
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);
options = odeset(options,'NonNegative', 10);
options = odeset(options, 'Events',@(s, y) eventIntegrationTraj(s, y, time0));

[s,y] = ode113(@(s,y) integrateTraectory(s, y, h0), [0 s_f], y0, options);
u=y(end, 1:4)';
v=y(end, 5:8)';
h_end=y(end, 9)'+h0;
ph=y(end, 19)';
tau=y(end, 10)';
ptau=y(end, 20)';
pv=y(end, 15:18);
t_start_fix=T_unit*(y(1, 10)-2*(y(1, 1:4)*y(1, 5:8)')/sqrt(-2*(y(1, 9)'+h0)))/(24*60*60);
t_end = T_unit*(tau-2*(u'*v)/sqrt(-2*h_end))/(24*60*60)-t_start_fix;
r_end=KS(u);
L_end = L_KS(u);
V_end = 2*sqrt(-2*h_end)*L_end*v/(norm(u)^2);

T_mars_days = 365.256363004*1.8808476;
n_M = floor((t_end+t_Mars_0)/T_mars_days);
angle_M = ((t_end+t_Mars_0)/T_mars_days-n_M)*2*pi;

rf = 1.52*[cos(angle_M) sin(angle_M) 0 0]';
Vf = ((mug/(1.52))^(1/2))*[cos(angle_M+pi/2) sin(angle_M+pi/2) 0 0]';
th = angle_M+n_M*2*pi;
uf=rToU(rf, th);
vf=vFromV(Vf,rf,mug,th);
%ЗАДАЧА ПРОЛЁТА
if case_traj == 1
    dis_p = [uf+u; pv';];
elseif case_traj == 2
    dis_p = [uf+u; vf+v;];
end
% if case_traj == 1
%     dis_p = [rf-r_end; pv';];
% elseif case_traj == 2
%     dis_p = [rf-r_end; Vf-V_end;];
% end
dis = modifier_f*norm(dis_p)^2;
end

