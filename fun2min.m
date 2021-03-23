function dis = fun2min(x, case_traj, t_start, r0, V0, planet_end, modifier_f, UorR,direction)
%UNTITLED Summary of this function goes here
% Функция расстояния до Марса, в квадратах координаты-скорости.
% Зависит от сопряжённых переменных в начальный момент времени

%Определяем безразмерные переменные
mug_0 = 132712.43994*(10^6)*(10^(3*3));
ae = 149597870700;
V_unit=sqrt(mug_0/ae);
T_earth = 365.256363004*3600*24;
T_unit = T_earth/(2*pi);

%Переходим к едининой гравитационной постоянной
mug=1;

%Задаём начальные условия на левом конце 
pu0=x(1:4)';
pv0=x(5:8)';
ph0=x(9);
pt0=x(10);
s_f=x(11)*2*pi;
phi=x(12)*2*pi;
u0 = rToU(r0, 0);
h0 = (norm(V0)^2)/2-mug/norm(r0);
v0 = vFromV(V0,r0,mug,0);
t0 = getEccentricAnomaly(r0(1:3),V0(1:3),mug);
y0 = cat(1, u0, v0, 0, t0, pu0, pv0, ph0, pt0)';

%Определяем параметры для оптимизатора
time0 = tic;
acc=1e-10;
options = odeset('AbsTol',acc);
options = odeset(options,'RelTol',acc);
options = odeset(options,'NonNegative', 10);
options = odeset(options, 'Events',@(s, y) eventIntegrationTraj(s, y, time0));
warning('off','all');
[s,y] = ode113(@(s,y) integrateTraectory(s, y, h0), [0 s_f], y0, options);

%Разбираем результат в конечный момент на переменные
u=y(end, 1:4)';
u2=u'*u;
v=y(end, 5:8)';
h_end=y(end, 9)'+h0;
tau=y(end, 10)';
pv=y(end, 15:18)';
ph=y(end, 19)';
ptau=y(end, 20)';
t_start_fix=T_unit*(y(1, 10)-2*(y(1, 1:4)*y(1, 5:8)')/sqrt(-2*(y(1, 9)'+h0)))/(24*60*60);
t_end = T_unit*(tau-2*(u'*v)/sqrt(-2*h_end))/(24*60*60)-t_start_fix;
r_end=KS(u);
L_end = L_KS(u);
V_end = 2*sqrt(-2*h_end)*L_end*v/(norm(u)^2);
aa_ks_end=L_end*(-(u2)*pv/(4*h_end) + v*(2*ph-(1/h_end)*pv'*v)+ptau*(u2)*u/((-2*h_end)^(3/2)));

%Получаем координату и скорость планеты в эфемеридах и поворачиваем систеу
%координат
[rf, Vf] = planetEphemeris(t_end+t_start,'SolarSystem',planet_end,'430');

eul = [0 pi/4 0];
rotmZYX = eul2rotm(eul);
rf = [rotmZYX*rf'; 0]/ae*1e+03;
Vf = [rotmZYX*Vf'; 0]/V_unit*1e+03;
%Получаем параметрические координату и скорость планеты
uf=rToU(rf, phi);
vf=vFromV(Vf,rf,mug,phi);

%Оптимизриуем по параметрическим координатам или по физическим
if UorR == 'u'
    %ЗАДАЧА ПРОЛЁТА или ЗАДАЧА СОПРОВОЖДЕНИЯ
    %direction - выбор положительного или отрицательного семейства
    if case_traj == 1
        dis_p = [uf+direction*u; aa_ks_end;];
    elseif case_traj == 2
        dis_p = [uf+direction*u; vf+direction*v;];
    end
elseif  UorR == 'r'
    %ЗАДАЧА ПРОЛЁТА или ЗАДАЧА СОПРОВОЖДЕНИЯ
    if case_traj == 1
        dis_p = [rf-r_end; pv';];
    elseif case_traj == 2
        dis_p = [rf-r_end; Vf-V_end;];
    end
end
%Сумма квадратов невязок, modifier_f влияет на сходимость
dis = modifier_f*norm(dis_p)^2;
end

