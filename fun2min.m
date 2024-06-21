function [c, ceq] = fun2min(x, st_f2m)
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
%st_f2m.case_traj, st_f2m.t_start, st_f2m.r0, st_f2m.V0, st_f2m.planet_end, st_f2m.modifier_f, st_f2m.UorR,st_f2m.decreaseNonPsysical,st_f2m.terminal_state, st_f2m.integration_acc, st_f2m.orbits, st_f2m.omega
%Задаём начальные условия на левом конце 
pu_0=x(1:4)';
pw_0=x(5:8)';
%ph_0=x(9);
if st_f2m.terminal_state == 's'
    s_f=x(9)*2*pi;
elseif st_f2m.terminal_state == 't'
    s_f=15*x(9)*2*pi;
end
t_end_0=x(9)*365.256363004;
phi=x(10)*2*pi;
%phi0=0;

if strcmp(st_f2m.UorR,'u_hat_v_hyp')
     hi = st_f2m.hi + x(11);
     dV_add = st_f2m.dV_value*st_f2m.rotmZYX*[cos(hi), sin(hi), 0]';
     st_f2m.V0 = st_f2m.V0+[dV_add; 0];
end


phi0=phi;
u_0 = rToU(st_f2m.r0, phi0);
u_b0=[u_0(4); -u_0(3);u_0(2);-u_0(1)];
h_0 = (norm(st_f2m.V0)^2)/2-mug/norm(st_f2m.r0);
w_0 = vFromV(st_f2m.V0,st_f2m.r0,mug, phi0);
%Условия трансверсальности на левом конце
f_left=get_target_g(u_0,w_0);
f_right=[st_f2m.r0(1:3);st_f2m.V0(1:3)];

if strcmp(st_f2m.UorR,'u_hat_v_hyp')
     vecPV_cartesian = a_reactive(u_0,w_0,pu_0,pw_0);
     if norm(vecPV_cartesian) > 0 && norm(dV_add)>0
         u_hat_v_hyp_left = 1e-09*cross(vecPV_cartesian(1:3)/norm(vecPV_cartesian(1:3)), dV_add(1:3)/norm(dV_add(1:3)));
     else
         u_hat_v_hyp_left = [0;0;0];
     end
end



%     ortdgduv_left=get_ortdgduv_Vhyp(u_0,w_0, st_f2m.Vp);
%     ortdgduv_left(1:4,1:4) = ortdgduv_left(1:4,1:4)*1e-13;
% else
%     ortdgduv_left=get_ortdgduv(u_0,w_0);
% end
ortdgduv_left=get_ortdgduv(u_0,w_0);
dis_p_tr_left=[[pu_0;pw_0]'*ortdgduv_left]';



%h0=-mug/(u0'*u0+4*st_f2m.V0'*st_f2m.V0)
%t0 = getEccentricAnomaly(st_f2m.r0(1:3),st_f2m.V0(1:3),mug);
%tau0=0;
tau0=1+2*u_0'*w_0/sqrt(-2*h_0);
y0 = cat(1, u_0, w_0, pu_0, pw_0, tau0)';
%st_f2m.t_start_fix=T_unit*(y0(10)-2*(y0(1:4)*y0(5:8)')/sqrt(-2*(y0(9)')))/(24*60*60);
%Определяем параметры для оптимизатора
time0 = tic;
%acc=1e-14;
options = odeset('AbsTol',st_f2m.integration_acc);
options = odeset(options,'RelTol',st_f2m.integration_acc);
options = odeset(options,'NonNegative', 17);
%максимальное время интегрирования
maxtime=10;
if st_f2m.terminal_state == 's'
    options = odeset(options, 'Events',@(s, y) eventIntegrationTraj(s, y, time0,maxtime));
elseif st_f2m.terminal_state == 't'
    options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopTime(s, y, time0,maxtime, t_end_0));
end

warning('off','all');

[s,y] = ode113(@(s,y) integrateTraectory(s, y), [0 s_f], y0, options);
int_s0sf = linspace(0, s(end), 100);
% time0 = tic;
% %максимальное время интегрирования
% maxtime=10;
% if st_f2m.terminal_state == 's'
%     options = odeset(options, 'Events',@(s, y) eventIntegrationTraj(s, y, time0,maxtime));
% elseif st_f2m.terminal_state == 't'
%     options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopTime(s, y, time0,maxtime, t_end_0, h0, st_f2m.t_start_fix));
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
w_end=y(end, 5:8)';
%h_end=y(end, 9)';
h_end=-mug/(u_end'*u_end+4*w_end'*w_end);
tau_end=y(end, 17)'-1;
pu_end=y(end, 9:12)';
pw_end=y(end, 13:16)';
%ph_end=y(end, 19)';
%ptau_end=y(end, 20)';

t_end = T_unit*(tau_end-2*(u_end'*w_end)/sqrt(-2*h_end))/(24*60*60);

if t_end > 10000
    t_end=10000;
end

L_end = L_KS(u_end);
r_end = L_end*u_end;
V_end = 2*sqrt(-2*h_end)*L_end*w_end/(norm(u_end)^2);
%a_ks_end=L_end*(-(u2)*pv_end/(4*h_end) + v_end*(2*ph_end-(1/h_end)*pv_end'*v_end)+ptau_end*(u2)*u_end/((-2*h_end)^(3/2)));

%Получаем координату и скорость планеты в эфемеридах и поворачиваем систеу
%координат
st.planet = st_f2m.planet_end;
st.mode = st_f2m.orbits;
st.delta_omega = st_f2m.omega;
st.a_rel = st_f2m.a_rel;
if st_f2m.terminal_state == 's'
    st.t = st_f2m.t_start+t_end;
    [rf, Vf] = planetModel(st);
elseif st_f2m.terminal_state == 't'
    st.t = st_f2m.t_start+t_end_0;
    [rf, Vf] = planetModel(st);
end



eul = [pi/12 pi/4 pi/12];
rotmZYX = eul2rotm(eul);
rf = [rotmZYX*rf'; 0]/ae*1e+03;
Vf = [rotmZYX*Vf'; 0]/V_unit*1e+03;

%Положение и скорость Земли для отладки
% [rf_e, Vf_e] = planetEphemeris(st_f2m.t_start+t_end,'SolarSystem','Earth','430');
% rf_e = [rotmZYX*rf_e'; 0]/ae*1e+03;
% Vf_e = [rotmZYX*Vf_e'; 0]/V_unit*1e+03;

%Получаем параметрические координату и скорость планеты
uf=rToU(rf, phi);
wf=vFromV(Vf,rf,mug,phi);
hf=norm(Vf)^2/2-mug/norm(rf);
if st_f2m.case_traj == 1
    g_left=get_target_g_u(u_end);
    g_right=rf(1:3);
    ortdgdu=get_ortdgdu(u_end);
elseif st_f2m.case_traj == 2
    g_left=get_target_g(u_end,w_end);
    g_right=[rf(1:3);Vf(1:3)];
    ortdgduv_right=get_ortdgduv(u_end,w_end);
end


%Оптимизриуем по параметрическим координатам или по физическим

if strcmp(st_f2m.UorR,'u_hat')
    %ЗАДАЧА ПРОЛЁТА или ЗАДАЧА СОПРОВОЖДЕНИЯ
    %direction - выбор положительного или отрицательного семейства
    if st_f2m.case_traj == 1
        dis_p_eqs = g_left-g_right;
        dis_p_tr = [pu_end'*ortdgdu; pw_end; ph_end];
        dis_p = [dis_p_eqs; dis_p_tr;];
    elseif st_f2m.case_traj == 2
        dis_p_eqs_right = g_left-g_right;
        %dis_p_eqs_left = f_left-f_right;
        dis_p_tr_right=[[pu_end;pw_end]'*ortdgduv_right]';
        dis_p=[dis_p_eqs_right;dis_p_tr_left];  
        %dis_p=dis_p_eqs_right;
    end
elseif  strcmp(st_f2m.UorR,'u')
        %ЗАДАЧА ПРОЛЁТА или ЗАДАЧА СОПРОВОЖДЕНИЯ
    %direction - выбор положительного или отрицательного семейства
    if st_f2m.case_traj == 1
        dis_p = [uf-u_end; pw_end];
    elseif st_f2m.case_traj == 2   
        dis_p = [uf-u_end; wf-w_end;dis_p_tr_left];
        %dis_p = [uf-u_end; wf-w_end];
    end
elseif  strcmp(st_f2m.UorR,'r')
    %ЗАДАЧА ПРОЛЁТА или ЗАДАЧА СОПРОВОЖДЕНИЯ
    if st_f2m.case_traj == 1
        dis_p = [rf-r_end; pw_end;];
    elseif st_f2m.case_traj == 2
        dis_p = [rf-r_end; Vf-V_end;];
    end
elseif  strcmp(st_f2m.UorR,'u_hat_v_hyp')
    if st_f2m.case_traj == 1
        dis_p_eqs = g_left-g_right;
        dis_p_tr = [pu_end'*ortdgdu; pw_end; ph_end];
        dis_p = [dis_p_eqs; dis_p_tr;];
    elseif st_f2m.case_traj == 2
        dis_p_eqs_right = g_left-g_right;
        %dis_p_eqs_left = f_left-f_right;
        dis_p=[dis_p_eqs_right;dis_p_tr_left;u_hat_v_hyp_left];  
        %dis_p=dis_p_eqs_right;
    end
end
%Сумма квадратов невязок, st_f2m.modifier_f влияет на сходимость
%dis = st_f2m.modifier_f*norm(dis_p)^2;
   
ceq = st_f2m.modifier_f*[dis_p;];
end

