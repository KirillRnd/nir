function [c, ceq] = fun2min(x, case_traj, t_start, r0, V0, planet_end, modifier_f, UorR,direction,terminal_state)
%UNTITLED Summary of this function goes here
% Функция расстояния до Марса, в квадратах координаты-скорости.
% Зависит от сопряжённых переменных в начальный момент времени
c =[];
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
if terminal_state == 's'
    s_f=x(11)*2*pi;
elseif terminal_state == 't'
    t_end=x(11)*365.256363004;
    s_f=1.5*x(11)*2*pi;
end
phi=x(12)*2*pi;
u0 = rToU(r0, 0);
h0 = (norm(V0)^2)/2-mug/norm(r0);
v0 = vFromV(V0,r0,mug,0);
t0 = getEccentricAnomaly(r0(1:3),V0(1:3),mug);
y0 = cat(1, u0, v0, 0, t0, pu0, pv0, ph0, pt0)';
t_start_fix=T_unit*(y0(10)-2*(y0(1:4)*y0(5:8)')/sqrt(-2*(y0(9)'+h0)))/(24*60*60);
%Определяем параметры для оптимизатора
time0 = tic;
acc=1e-10;
options = odeset('AbsTol',acc);
options = odeset(options,'RelTol',acc);
options = odeset(options,'NonNegative', 10);
if terminal_state == 's'
    options = odeset(options, 'Events',@(s, y) eventIntegrationTraj(s, y, time0));
elseif terminal_state == 't'
    options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopTime(s, y, time0, t_end, h0, t_start_fix));
end

warning('off','all');
int_s0sf = linspace(0, s_f, 10);
[s,y] = ode113(@(s,y) integrateTraectory(s, y, h0), int_s0sf, y0, options);

ub_matrix=[y(:, 4), -y(:, 3), y(:, 2), -y(:, 1)];
v_matrix=[y(:, 5), y(:, 6), y(:, 7), y(:, 8)];
pv_matrix=[y(:, 15), y(:, 16), y(:, 17), y(:, 18)];
ubTv=diag(ub_matrix*v_matrix');
ubTpv=diag(ub_matrix*pv_matrix');
%Разбираем результат в конечный момент на переменные
u_end=y(end, 1:4)';
u2=u_end'*u_end;
v_end=y(end, 5:8)';
h_end=y(end, 9)'+h0;
tau=y(end, 10)';
pu_end=y(end, 11:14)';
pv_end=y(end, 15:18)';
ph=y(end, 19)';
ptau=y(end, 20)';

t_end = T_unit*(tau-2*(u_end'*v_end)/sqrt(-2*h_end))/(24*60*60)-t_start_fix;
r_end=KS(u_end);
L_end = L_KS(u_end);
V_end = 2*sqrt(-2*h_end)*L_end*v_end/(norm(u_end)^2);
a_ks_end=L_end*(-(u2)*pv_end/(4*h_end) + v_end*(2*ph-(1/h_end)*pv_end'*v_end)+ptau*(u2)*u_end/((-2*h_end)^(3/2)));

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
hf=norm(Vf)^2/2-mug/norm(rf);
x=rf(1);
y=rf(2);
z=rf(3);

R1=sqrt(0.5*(norm(rf)+x));
R2=sqrt(0.5*(norm(rf)-x));

phi_end=atan2(u_end(4),u_end(1));
gamma_end=atan2(u_end(3),u_end(2));
theta_end=phi_end+gamma_end;
if theta_end<0
    theta_end=theta_end+2*pi;
end
thetaf = atan2(z,y);
if thetaf<0
    thetaf=thetaf+2*pi;
end
bil=@(x)[x(4); -x(3); x(2); -x(1)];
C1=@(x)R1*x(1)+R2*(x(2)*cos(thetaf)+x(3)*sin(thetaf));
C2=@(x)R1*x(4)+R2*(x(2)*sin(thetaf)-x(3)*cos(thetaf));
%вычисляем невязку до u_hat
gu_left=[u_end(1)^2+u_end(4)^2;
    u_end(2)^2+u_end(3)^2;
    theta_end];
gu_right=[R1^2;R2^2;thetaf];

dgdu = get_dgdu(uf);
ort_u=null(dgdu);

pu_proj_coef=dgdu*dgdu'\dgdu*pu_end;
pu_proj=dgdu'*pu_proj_coef;
%вычисляем проекцию на линейное подпространство 
pu_ort=pu_end-pu_proj;
%вычисляем невязку до v_hat
gv_left=get_gv(v_end,V_end)';

theta_end_v=gv_left(3);
if theta_end_v<0
    theta_end_v=theta_end_v+2*pi;
end
gv_left(3)=theta_end_v;
gv_right=[Vf'*Vf*norm(rf)/(-8*hf); x/(-8*hf); thetaf];
dgdv=get_dgdv(vf,Vf);
ort_v=null(dgdv);
pv_proj_coef=dgdv*dgdv'\dgdv*pv_end;
pv_proj=dgdv'*pv_proj_coef;
%вычисляем проекцию на линейное подпространство 
pv_ort=pv_end-pv_proj;

F = [[Vf(1) Vf(2) Vf(3) Vf(4)];
    [Vf(2) -Vf(1) -Vf(4) Vf(3)];
    [Vf(3) Vf(4) -Vf(1) -Vf(2)];
    [-Vf(4) Vf(3) -Vf(2) Vf(1)]];

pu_ort_eq=C1(bil(pu_end))^2+C2(bil(pu_end))^2;
pv_ort_eq=C1(bil(F'*pv_end))^2+C2(bil(F'*pv_end))^2;
%Оптимизриуем по параметрическим координатам или по физическим
modifier_f_2=0;
if strcmp(UorR,'u_hat')
    %ЗАДАЧА ПРОЛЁТА или ЗАДАЧА СОПРОВОЖДЕНИЯ
    %direction - выбор положительного или отрицательного семейства
    if case_traj == 1
        dis_p = [gu_left-gu_right; a_ks_end];
    elseif case_traj == 2
        %dis_p = [gu_left-gu_right; gv_left-gv_right;pu_ort_eq;pv_ort_eqt];
        dis_p_eqs = [gu_left-gu_right; gv_left-gv_right];
        %dis_p_tr = [C1(bil(pu_end));C2(bil(pu_end));C1(bil(F'*pv_end));C2(bil(F'*pv_end))];
        dis_p_tr=[pu_end'*ort_u;pv_end'*ort_v];
        dis_p=[modifier_f*dis_p_eqs;modifier_f_2*dis_p_tr];
    end
elseif  strcmp(UorR,'u')
        %ЗАДАЧА ПРОЛЁТА или ЗАДАЧА СОПРОВОЖДЕНИЯ
    %direction - выбор положительного или отрицательного семейства
    if case_traj == 1
        dis_p = [uf-u_end; a_ks_end];
    elseif case_traj == 2
        dis_p = modifier_f*[uf-u_end; vf-v_end];
    end
elseif  strcmp(UorR,'r')
    %ЗАДАЧА ПРОЛЁТА или ЗАДАЧА СОПРОВОЖДЕНИЯ
    if case_traj == 1
        dis_p = [rf-r_end; pv_end;];
    elseif case_traj == 2
        dis_p = [rf-r_end; Vf-V_end;];
    end
end
%Сумма квадратов невязок, modifier_f влияет на сходимость
%dis = modifier_f*norm(dis_p)^2;

ceq = dis_p;
end

