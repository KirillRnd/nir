function [c, ceq] = fun2min(x, case_traj, t_start, r0, V0, planet_end, modifier_f, UorR,decreaseNonPsysical,terminal_state, integration_acc,amax)
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
pw0=x(5:8)';
ph0=x(9);
pt0=x(10);
if terminal_state == 's'
    s_f=x(11)*2*pi;
elseif terminal_state == 't'
    s_f=15*x(11)*2*pi;
end
t_end_0=x(11)*365.256363004;
phi=x(12)*2*pi;
phi0=0;
u0 = rToU(r0, phi0);
u_b0=[u0(4); -u0(3);u0(2);-u0(1)];
h0 = (norm(V0)^2)/2-mug/norm(r0);
w0 = vFromV(V0,r0,mug, phi0);
%h0=-mug/(u0'*u0+4*v0'*v0)
%t0 = getEccentricAnomaly(r0(1:3),V0(1:3),mug);
%tau0=0;
tau0=2*u0'*w0/sqrt(-2*h0);
y0 = cat(1, u0, w0, h0, tau0, pu0, pw0, ph0, pt0)';
%t_start_fix=T_unit*(y0(10)-2*(y0(1:4)*y0(5:8)')/sqrt(-2*(y0(9)')))/(24*60*60);
%Определяем параметры для оптимизатора
time0 = tic;
%acc=1e-14;
options = odeset('AbsTol',integration_acc);
options = odeset(options,'RelTol',integration_acc);
options = odeset(options,'NonNegative', 10);
%максимальное время интегрирования
maxtime=10;
if terminal_state == 's'
    options = odeset(options, 'Events',@(s, y) eventIntegrationTraj(s, y, time0,maxtime));
elseif terminal_state == 't'
    options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopTime(s, y, time0,maxtime, t_end_0, h0, t_start_fix));
end

warning('off','all');

[s,y] = ode113(@(s,y) integrateTraectory(s, y, amax), [0 s_f], y0, options);
int_s0sf = linspace(0, s(end), 100);
% time0 = tic;
% %максимальное время интегрирования
% maxtime=10;
% if terminal_state == 's'
%     options = odeset(options, 'Events',@(s, y) eventIntegrationTraj(s, y, time0,maxtime));
% elseif terminal_state == 't'
%     options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopTime(s, y, time0,maxtime, t_end_0, h0, t_start_fix));
% end
% 
% [s,y] = ode113(@(s,y) integrateTraectory(s, y, h0), int_s0sf, y0, options);

% ub_matrix=[y(:, 4), -y(:, 3), y(:, 2), -y(:, 1)];
% v_matrix=[y(:, 5), y(:, 6), y(:, 7), y(:, 8)];
% pv_matrix=[y(:, 15), y(:, 16), y(:, 17), y(:, 18)];
% ubTv=diag(ub_matrix*v_matrix');
% ubTpv=diag(ub_matrix*pv_matrix');
%Разбираем результат в конечный момент на переменные
u_end=y(end, 1:4)';
u_b_end=[u_end(4); -u_end(3);u_end(2);-u_end(1)];
u2=u_end'*u_end;
v_end=y(end, 5:8)';
h_end=y(end, 9)';
tau_end=y(end, 10)';
pu_end=y(end, 11:14)';
pv_end=y(end, 15:18)';
ph_end=y(end, 19)';
ptau_end=y(end, 20)';

t_end = T_unit*(tau_end-2*(u_end'*v_end)/sqrt(-2*h_end))/(24*60*60);
r_end=KS(u_end);
L_end = L_KS(u_end);
V_end = 2*sqrt(-2*h_end)*L_end*v_end/(norm(u_end)^2);
%a_ks_end=L_end*(-(u2)*pv_end/(4*h_end) + v_end*(2*ph_end-(1/h_end)*pv_end'*v_end)+ptau_end*(u2)*u_end/((-2*h_end)^(3/2)));

%Получаем координату и скорость планеты в эфемеридах и поворачиваем систеу
%координат
if terminal_state == 's'
    [rf, Vf] = planetEphemeris(t_start+t_end,'SolarSystem',planet_end,'430');
elseif terminal_state == 't'
    [rf, Vf] = planetEphemeris(t_start+t_end_0,'SolarSystem',planet_end,'430');
end


eul = [0 pi/4 0];
rotmZYX = eul2rotm(eul);
rf = [rotmZYX*rf'; 0]/ae*1e+03;
Vf = [rotmZYX*Vf'; 0]/V_unit*1e+03;
%Получаем параметрические координату и скорость планеты
uf=rToU(rf, phi);
vf=vFromV(Vf,rf,mug,phi);
hf=norm(Vf)^2/2-mug/norm(rf);
if case_traj == 1
    g_left=get_target_g_u(u_end);
    g_right=rf(1:3);
    ortdgdu=get_ortdgdu(u_end);
elseif case_traj == 2
    g_left=get_target_gh(u_end,v_end,h_end);
    g_right=[rf(1:3);Vf(1:3);hf];
    ortdgduvh=get_ortdgduvh(u_end,v_end,h_end);
end

%Оптимизриуем по параметрическим координатам или по физическим

if strcmp(UorR,'u_hat')
    %ЗАДАЧА ПРОЛЁТА или ЗАДАЧА СОПРОВОЖДЕНИЯ
    %direction - выбор положительного или отрицательного семейства
    if case_traj == 1
        dis_p_eqs = g_left-g_right;
        dis_p_tr = [pu_end'*ortdgdu]';
        dis_p = [dis_p_eqs(1:3); dis_p_tr; pv_end];
    elseif case_traj == 2
        dis_p_eqs = g_left-g_right;
        dis_p_tr=[[pu_end;pv_end;ph_end]'*ortdgduvh]';
        dis_p=[dis_p_eqs;dis_p_tr];        
    end
elseif  strcmp(UorR,'u')
        %ЗАДАЧА ПРОЛЁТА или ЗАДАЧА СОПРОВОЖДЕНИЯ
    %direction - выбор положительного или отрицательного семейства
    if case_traj == 1
        dis_p = [uf-u_end; pv_end];
    elseif case_traj == 2
        dis_p = [uf-u_end; vf-v_end;];
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
if decreaseNonPsysical == 1
    % time0 = tic;
    %максимальное время интегрирования
    maxtime=10;
    if terminal_state == 's'
        options = odeset(options, 'Events',@(s, y) eventIntegrationTraj(s, y, time0,maxtime));
    elseif terminal_state == 't'
        options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopTime(s, y, time0,maxtime, t_end_0, h0, t_start_fix));
    end
    
    [s,y] = ode113(@(s,y) integrateTraectory(s, y, h0), int_s0sf, y0, options);

    ub_matrix=[y(:, 4), -y(:, 3), y(:, 2), -y(:, 1)];
    v_matrix=[y(:, 5), y(:, 6), y(:, 7), y(:, 8)];
    pv_matrix=[y(:, 15), y(:, 16), y(:, 17), y(:, 18)];
    ubTv=diag(ub_matrix*v_matrix');
    ubTpv=diag(ub_matrix*pv_matrix');
    dis_ubT = zeros(100,1);
    dis_ubT(1:length(ubTv))=ubTpv;
    %dis_ubT(101:length(ubTpv)+100)=ubTpv;
    dis_p=[dis_p;dis_ubT];
end   
ceq = modifier_f*[dis_p;];
end

