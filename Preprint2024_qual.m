% -*- coding: utf-8 -*-
% в этом скрипте изучаем пространственные перелёты с помощью афинных
% преобразований
clear;
load('mat-files/Preprint2024_1_Mars.mat')
load('mat-files/Preprint2024_1_Venus.mat', 't_end_Venus')
%%
k = 1.5;
%dvx_coef_range=-1e-5:1e-6:1e-5;
dvx_coef_range = [logspace(-12, -0, 100)];
%dvy_coef_range=-0.10:0.01:0.10;
N_count_1 = length(dvx_coef_range);
r_err_massive = zeros(N_count_1,1);
v_err_massive = zeros(N_count_1,1);

r_err_massive_rot = zeros(N_count_1,1);
v_err_massive_rot = zeros(N_count_1,1);

rf_target = [k 0 0];
vf_target = [0 1/sqrt(k) 0];

[pr_0,pv_0] = get_adjoint(k);
t_end = get_t_end(k);
t0 = 0;
start_pos = [1, 0, 0];
start_vel = [0, 1, 0];

for i = 1:N_count_1
    start_pos_2 = start_pos;
    start_vel_2 = start_vel;
    start_vel_2(1) = dvx_coef_range(i);
    %R-преобразование
    R = calculateRotMatrix(start_pos_2,start_vel_2);
    
    pr = R^(-1)*pr_0;
    pv = R^(-1)*pv_0;
    y0 = cat(2,start_pos_2,start_vel_2,pr',pv')';
    tspan = linspace(t0,t0+t_end, 10);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    
    rf = y(end,1:3);
    vf = y(end,4:6);
    
    res_r = norm(rf_target-rf);
    res_v = norm(vf_target-vf);
    
    %обычный поворот
    y0 = cat(2,start_pos_2,start_vel_2,pr_0',pv_0')';
    tspan = linspace(t0,t0+t_end, 10);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    
    rf = y(end,1:3);
    vf = y(end,4:6);
    
    res_r_rot = norm(rf_target-rf);
    res_v_rot = norm(vf_target-vf);


    r_err_massive(i) = res_r;
    v_err_massive(i) = res_v;

    r_err_massive_rot(i) = res_r_rot;
    v_err_massive_rot(i) = res_v_rot;
end
%%
figure(1)
plot(dvx_coef_range, r_err_massive,'r--', 'LineWidth', 1.0,'DisplayName', 'невязка положения, R-преобразование')
hold on;
plot(dvx_coef_range, r_err_massive_rot,'b--', 'LineWidth', 1.0,'DisplayName', 'невязка положения, матрица поворота')
plot(dvx_coef_range, v_err_massive,'r', 'LineWidth', 1.0,'DisplayName', 'невязка скорости, R-преобразование')
plot(dvx_coef_range, v_err_massive_rot,'b', 'LineWidth', 1.0,'DisplayName', 'невязка скорости, матрица поворота')
hold off;
legend('Location','best');
grid;
xlabel('\delta v_x, безразм.')
ylabel('Невязка, безразм.')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca, 'FontSize', 11)
%%

%dvx_coef_range=-1e-5:1e-6:1e-5;
dvy_coef_range = [logspace(-12, -0, 100)];
%dvy_coef_range=-0.10:0.01:0.10;
N_count_1 = length(dvy_coef_range);
r_err_massive = zeros(N_count_1,1);
v_err_massive = zeros(N_count_1,1);

r_err_massive_rot = zeros(N_count_1,1);
v_err_massive_rot = zeros(N_count_1,1);

[pr_0,pv_0] = get_adjoint(k);
t_end = get_t_end(k);
t0 = 0;
start_pos = [1, 0, 0];
start_vel = [0, 1, 0];

for i = 1:N_count_1
    start_pos_2 = start_pos;
    start_vel_2 = start_vel;
    start_vel_2(2) = 1+dvy_coef_range(i);
    %R-преобразование
    R = calculateRotMatrix(start_pos_2,start_vel_2);
    
    pr = R^(-1)*pr_0;
    pv = R^(-1)*pv_0;
    y0 = cat(2,start_pos_2,start_vel_2,pr',pv')';
    tspan = linspace(t0,t0+t_end, 10);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    
    rf = y(end,1:3);
    vf = y(end,4:6);
    
    res_r = norm(rf_target-rf);
    res_v = norm(vf_target-vf);
    
    %обычный поворот
    y0 = cat(2,start_pos_2,start_vel_2,pr_0',pv_0')';
    tspan = linspace(t0,t0+t_end, 10);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    
    rf = y(end,1:3);
    vf = y(end,4:6);
    
    res_r_rot = norm(rf_target-rf);
    res_v_rot = norm(vf_target-vf);


    r_err_massive(i) = res_r;
    v_err_massive(i) = res_v;

    r_err_massive_rot(i) = res_r_rot;
    v_err_massive_rot(i) = res_v_rot;
end
%%
figure(2)
plot(dvy_coef_range, r_err_massive,'r--', 'LineWidth', 1.0,'DisplayName', 'невязка положения, R-преобразование')
hold on;
plot(dvy_coef_range, r_err_massive_rot,'b--', 'LineWidth', 1.0,'DisplayName', 'невязка положения, матрица поворота')
plot(dvy_coef_range, v_err_massive,'r', 'LineWidth', 1.0,'DisplayName', 'невязка скорости, R-преобразование')
plot(dvy_coef_range, v_err_massive_rot,'b', 'LineWidth', 1.0,'DisplayName', 'невязка скорости, матрица поворота')
hold off;
legend('Location','best');
grid;
xlabel('\delta v_y, безразм.')
ylabel('Невязка, безразм.')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca, 'FontSize', 11)
%%
%dvx_coef_range=-1e-5:1e-6:1e-5;
drx_coef_range = [logspace(-12, -0, 100)];
%dvy_coef_range=-0.10:0.01:0.10;
N_count_1 = length(drx_coef_range);
r_err_massive = zeros(N_count_1,1);
v_err_massive = zeros(N_count_1,1);

r_err_massive_rot = zeros(N_count_1,1);
v_err_massive_rot = zeros(N_count_1,1);

[pr_0,pv_0] = get_adjoint(k);
t_end = get_t_end(k);
t0 = 0;
start_pos = [1, 0, 0];
start_vel = [0, 1, 0];

for i = 1:N_count_1
    start_pos_2 = start_pos;
    start_pos(1) = 1 + drx_coef_range(i);
    start_vel_2 = start_vel;
    %start_vel_2(1) = dvx_coef_range(i);
    %R-преобразование
    R = calculateRotMatrix(start_pos_2,start_vel_2);
    
    pr = R^(-1)*pr_0;
    pv = R^(-1)*pv_0;
    y0 = cat(2,start_pos_2,start_vel_2,pr',pv')';
    tspan = linspace(t0,t0+t_end, 10);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    
    rf = y(end,1:3);
    vf = y(end,4:6);
    
    res_r = norm(rf_target-rf);
    res_v = norm(vf_target-vf);
    
    %обычный поворот
    y0 = cat(2,start_pos_2,start_vel_2,pr_0',pv_0')';
    tspan = linspace(t0,t0+t_end, 10);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    
    rf = y(end,1:3);
    vf = y(end,4:6);
    
    res_r_rot = norm(rf_target-rf);
    res_v_rot = norm(vf_target-vf);


    r_err_massive(i) = res_r;
    v_err_massive(i) = res_v;

    r_err_massive_rot(i) = res_r_rot;
    v_err_massive_rot(i) = res_v_rot;
end
%%
figure(3)
plot(drx_coef_range, r_err_massive,'r--', 'LineWidth', 1.0,'DisplayName', 'невязка положения, R-преобразование')
hold on;
plot(drx_coef_range, r_err_massive_rot,'b--', 'LineWidth', 1.0,'DisplayName', 'невязка положения, матрица поворота')
plot(drx_coef_range, v_err_massive,'r', 'LineWidth', 1.0,'DisplayName', 'невязка скорости, R-преобразование')
plot(drx_coef_range, v_err_massive_rot,'b', 'LineWidth', 1.0,'DisplayName', 'невязка скорости, матрица поворота')
hold off;
legend('Location','best');
grid;
xlabel('\delta r_x, безразм.')
ylabel('Невязка, безразм.')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca, 'FontSize', 11)
%%
%dvx_coef_range=-1e-5:1e-6:1e-5;
dry_coef_range = [logspace(-12, -1, 100)];
%dvy_coef_range=-0.10:0.01:0.10;
N_count_1 = length(dry_coef_range);
r_err_massive = zeros(N_count_1,1);
v_err_massive = zeros(N_count_1,1);

r_err_massive_no = zeros(N_count_1,1);
v_err_massive_no = zeros(N_count_1,1);

r_err_massive_rot = zeros(N_count_1,1);
v_err_massive_rot = zeros(N_count_1,1);

[pr_0,pv_0] = get_adjoint(k);
t_end = get_t_end(k);
t0 = 0;
start_pos = [1, 0, 0];
start_vel = [0, 1, 0];

for i = 1:N_count_1
    start_pos_2 = start_pos;
    start_pos(2) = dry_coef_range(i);
    start_vel_2 = start_vel;
    %start_vel_2(1) = dvx_coef_range(i);
    %R-преобразование
    R = calculateRotMatrix(start_pos_2,start_vel_2);
    
    pr = R^(-1)*pr_0;
    pv = R^(-1)*pv_0;
    y0 = cat(2,start_pos_2,start_vel_2,pr',pv')';
    tspan = linspace(t0,t0+t_end, 10);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    
    rf = y(end,1:3);
    vf = y(end,4:6);
    
    res_r = norm(rf_target-rf);
    res_v = norm(vf_target-vf);
    
    %обычный поворот


    r_norm = start_pos_2/norm(start_pos_2);
    v_proj = start_vel_2-(start_vel_2*r_norm')*r_norm;
    v_norm = v_proj/norm(v_proj);
    R_simple = calculateRotMatrix(r_norm,v_norm);

    pr = R_simple^(-1)*pr_0;
    pv = R_simple^(-1)*pv_0;

    y0 = cat(2,start_pos_2,start_vel_2,pr',pv')';
    tspan = linspace(t0,t0+t_end, 10);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    
    rf = y(end,1:3);
    vf = y(end,4:6);
    
    res_r_rot = norm(rf_target-rf);
    res_v_rot = norm(vf_target-vf);

    %нет коррекции
    y0 = cat(2,start_pos_2,start_vel_2,pr_0',pv_0')';
    tspan = linspace(t0,t0+t_end, 10);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    
    rf = y(end,1:3);
    vf = y(end,4:6);
    
    res_r_no = norm(rf_target-rf);
    res_v_no = norm(vf_target-vf);

    %сохраняем результат
    r_err_massive(i) = res_r;
    v_err_massive(i) = res_v;

    r_err_massive_no(i) = res_r_no;
    v_err_massive_no(i) = res_v_no;

    r_err_massive_rot(i) = res_r_rot;
    v_err_massive_rot(i) = res_v_rot;
end
%%
figure(4)
plot(dry_coef_range, r_err_massive,'r--', 'LineWidth', 1.0,'DisplayName', 'невязка положения, R-преобразование')
hold on;
plot(dry_coef_range, r_err_massive_rot,'b--', 'LineWidth', 1.0,'DisplayName', 'невязка положения, матрица поворота')
plot(dry_coef_range, r_err_massive_no,'k--', 'LineWidth', 1.0,'DisplayName', 'невязка положения, без коррекции')
plot(dry_coef_range, v_err_massive,'r', 'LineWidth', 1.0,'DisplayName', 'невязка скорости, R-преобразование')
plot(dry_coef_range, v_err_massive_rot,'b', 'LineWidth', 1.0,'DisplayName', 'невязка скорости, матрица поворота')
plot(dry_coef_range, v_err_massive_no,'k', 'LineWidth', 1.0,'DisplayName', 'невязка скорости, без коррекции')
hold off;
legend('Location','best');
grid;
xlabel('\delta r_y, безразм.')
ylabel('Невязка, безразм.')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca, 'FontSize', 11)
%%
%dvx_coef_range=-1e-5:1e-6:1e-5;
drz_coef_range = [logspace(-12, -0, 100)];
%dvy_coef_range=-0.10:0.01:0.10;
N_count_1 = length(drz_coef_range);
r_err_massive = zeros(N_count_1,1);
v_err_massive = zeros(N_count_1,1);

r_err_massive_no = zeros(N_count_1,1);
v_err_massive_no = zeros(N_count_1,1);

r_err_massive_rot = zeros(N_count_1,1);
v_err_massive_rot = zeros(N_count_1,1);

[pr_0,pv_0] = get_adjoint(k);
t_end = get_t_end(k);
t0 = 0;
start_pos = [1, 0, 0];
start_vel = [0, 1, 0];

for i = 1:N_count_1
    start_pos_2 = start_pos;
    start_pos(3) = drz_coef_range(i);
    start_vel_2 = start_vel;
    %start_vel_2(1) = dvx_coef_range(i);
    %R-преобразование
    R = calculateRotMatrix(start_pos_2,start_vel_2);
    
    pr = R^(-1)*pr_0;
    pv = R^(-1)*pv_0;
    y0 = cat(2,start_pos_2,start_vel_2,pr',pv')';
    tspan = linspace(t0,t0+t_end, 10);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    
    rf = y(end,1:3);
    vf = y(end,4:6);
    
    res_r = norm(rf_target-rf);
    res_v = norm(vf_target-vf);
    
    %обычный поворот


    r_norm = start_pos_2/norm(start_pos_2);
    v_proj = start_vel_2-(start_vel_2*r_norm')*r_norm;
    v_norm = v_proj/norm(v_proj);
    R_simple = calculateRotMatrix(r_norm,v_norm);

    pr = R_simple^(-1)*pr_0;
    pv = R_simple^(-1)*pv_0;

    y0 = cat(2,start_pos_2,start_vel_2,pr',pv')';
    tspan = linspace(t0,t0+t_end, 10);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    
    rf = y(end,1:3);
    vf = y(end,4:6);
    
    res_r_rot = norm(rf_target-rf);
    res_v_rot = norm(vf_target-vf);

    %нет коррекции
    y0 = cat(2,start_pos_2,start_vel_2,pr_0',pv_0')';
    tspan = linspace(t0,t0+t_end, 10);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    
    rf = y(end,1:3);
    vf = y(end,4:6);
    
    res_r_no = norm(rf_target-rf);
    res_v_no = norm(vf_target-vf);

    %сохраняем результат
    r_err_massive(i) = res_r;
    v_err_massive(i) = res_v;

    r_err_massive_no(i) = res_r_no;
    v_err_massive_no(i) = res_v_no;

    r_err_massive_rot(i) = res_r_rot;
    v_err_massive_rot(i) = res_v_rot;
end
%%
figure(5)
plot(drz_coef_range, r_err_massive,'r--', 'LineWidth', 1.0,'DisplayName', 'невязка положения, R-преобразование')
hold on;
plot(drz_coef_range, r_err_massive_rot,'b--', 'LineWidth', 1.0,'DisplayName', 'невязка положения, матрица поворота')
plot(drz_coef_range, r_err_massive_no,'g--', 'LineWidth', 1.0,'DisplayName', 'невязка положения, без коррекции')
plot(drz_coef_range, v_err_massive,'r', 'LineWidth', 1.0,'DisplayName', 'невязка скорости, R-преобразование')
plot(drz_coef_range, v_err_massive_rot,'b', 'LineWidth', 1.0,'DisplayName', 'невязка скорости, матрица поворота')
plot(drz_coef_range, v_err_massive_no,'g', 'LineWidth', 1.0,'DisplayName', 'невязка скорости, без коррекции')
hold off;
legend('Location','best');
grid;
xlabel('\delta r_z, безразм.')
ylabel('Невязка, безразм.')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca, 'FontSize', 11)
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off'); 
%%
%dvx_coef_range=-1e-5:1e-6:1e-5;
dvz_coef_range = [logspace(-12, -0, 100)];
%dvy_coef_range=-0.10:0.01:0.10;
N_count_1 = length(dvz_coef_range);
r_err_massive = zeros(N_count_1,1);
v_err_massive = zeros(N_count_1,1);

r_err_massive_no = zeros(N_count_1,1);
v_err_massive_no = zeros(N_count_1,1);

r_err_massive_rot = zeros(N_count_1,1);
v_err_massive_rot = zeros(N_count_1,1);

[pr_0,pv_0] = get_adjoint(k);
t_end = get_t_end(k);
t0 = 0;
start_pos = [1, 0, 0];
start_vel = [0, 1, 0];

for i = 1:N_count_1
    start_pos_2 = start_pos;
    %start_pos(3) = drz_coef_range(i);
    start_vel_2 = start_vel;
    start_vel_2(3) = dvz_coef_range(i);
    %R-преобразование
    R = calculateRotMatrix(start_pos_2,start_vel_2);
    
    pr = R^(-1)*pr_0;
    pv = R^(-1)*pv_0;
    y0 = cat(2,start_pos_2,start_vel_2,pr',pv')';
    tspan = linspace(t0,t0+t_end, 10);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    
    rf = y(end,1:3);
    vf = y(end,4:6);
    
    res_r = norm(rf_target-rf);
    res_v = norm(vf_target-vf);
    
    %обычный поворот


    r_norm = start_pos_2/norm(start_pos_2);
    v_proj = start_vel_2-(start_vel_2*r_norm')*r_norm;
    v_norm = v_proj/norm(v_proj);
    R_simple = calculateRotMatrix(r_norm,v_norm);

    pr = R_simple^(-1)*pr_0;
    pv = R_simple^(-1)*pv_0;

    y0 = cat(2,start_pos_2,start_vel_2,pr',pv')';
    tspan = linspace(t0,t0+t_end, 10);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    
    rf = y(end,1:3);
    vf = y(end,4:6);
    
    res_r_rot = norm(rf_target-rf);
    res_v_rot = norm(vf_target-vf);

    %нет коррекции
    y0 = cat(2,start_pos_2,start_vel_2,pr_0',pv_0')';
    tspan = linspace(t0,t0+t_end, 10);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    
    rf = y(end,1:3);
    vf = y(end,4:6);
    
    res_r_no = norm(rf_target-rf);
    res_v_no = norm(vf_target-vf);

    %сохраняем результат
    r_err_massive(i) = res_r;
    v_err_massive(i) = res_v;

    r_err_massive_no(i) = res_r_no;
    v_err_massive_no(i) = res_v_no;

    r_err_massive_rot(i) = res_r_rot;
    v_err_massive_rot(i) = res_v_rot;
end
%%
figure(6)
plot(dvz_coef_range, r_err_massive,'r--', 'LineWidth', 1.0,'DisplayName', 'невязка положения, R-преобразование')
hold on;
plot(dvz_coef_range, r_err_massive_rot,'b--', 'LineWidth', 1.0,'DisplayName', 'невязка положения, матрица поворота')
plot(dvz_coef_range, r_err_massive_no,'g--', 'LineWidth', 1.0,'DisplayName', 'невязка положения, без коррекции')
plot(dvz_coef_range, v_err_massive,'r', 'LineWidth', 1.0,'DisplayName', 'невязка скорости, R-преобразование')
plot(dvz_coef_range, v_err_massive_rot,'b', 'LineWidth', 1.0,'DisplayName', 'невязка скорости, матрица поворота')
plot(dvz_coef_range, v_err_massive_no,'g', 'LineWidth', 1.0,'DisplayName', 'невязка скорости, без коррекции')
hold off;
legend('Location','best');
grid;
xlabel('\delta v_z, безразм.')
ylabel('Невязка, безразм.')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca, 'FontSize', 11)
function res = fsolve_find_i(i_matrix,i_target, start_pos,start_vel,zf,t0,t_end,O_target)
%изменение наклонения
delta_i = i_matrix;
n_rot_i = [1; 0; 0];
q_rot_i = [cos(delta_i/2); n_rot_i*sin(delta_i/2)];
R_i0 = quat2rotm(q_rot_i');
R_i0(1:3,3)=R_i0(1:3,3)*R_i0(2,2);
R_i0(1:3,2)=R_i0(1:3,2)*R_i0(2,2);
%восходящий узел
%O_target = 0*pi/6;
delta_O = O_target-pi;
n_rot = [0; 0; 1];
q_rot = [cos(delta_O/2); n_rot*sin(delta_O/2)];
R_rot = quat2rotm(q_rot');
R_iO=R_rot*R_i0*R_rot';
[i,O,J,AN] = integrate_iO(R_iO, start_pos,start_vel,zf,t0,t_end);
res = i-i_target;
end
function [i,O,p,a,J,AN] = integrate_iO(R_iO, start_pos,start_vel,zf,t0,t_end)
pr = R_iO^(-1)*zf(1:3);
pv = R_iO^(-1)*zf(4:6);
y0 = cat(2,start_pos,start_vel,pr',pv',0,0)';
tspan = linspace(t0,t0+t_end, 1000);
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
[t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);

%J = y(end,13);       %функционал
%AN = y(end,14);%угловая дальность
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(y(end, 1:3)',y(end, 4:6)',1);
[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
end
function [ex,ey,lonPer,p,a] = integrate_iO_equinoctial(R_iO, start_pos,start_vel,zf,t0,t_end)
pr = R_iO^(-1)*zf(1:3);
pv = R_iO^(-1)*zf(4:6);

load('mat-files/Preprint2024_1_Mars.mat','dxdX');
P = dxdX'*[pr;pv];
y0 = [1;0;0;0;0;0;P];
tspan = linspace(t0,t0+t_end, 1000);
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
mug = 1;
[t,y] = ode113(@(t,y) equitoctial_adjoint_for_integration(y,mug), tspan,y0,options);

traj = arrayfun(@(p,ex,ey,ix,iy,L) equitoctial2decart([p;ex;ey;ix;iy;L], mug),...
    y(:, 1),y(:, 2),y(:, 3),y(:, 4),y(:, 5),y(:, 6),'UniformOutput',false);
traj = cell2mat(traj')';

[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(traj(end, 1:3)',traj(end, 4:6)',1);
[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
end
function t_end = get_t_end(k)
load('mat-files/Preprint2024_1_universal_2.mat','k_coef_for_spline','T_40_array')
t_end = spline(k_coef_for_spline,T_40_array',k);
end
function [pr_0,pv_0] = get_adjoint(k)
load('mat-files/Preprint2024_1_universal_2.mat','k_coef_for_spline','P_40_array')
z = spline(k_coef_for_spline,P_40_array',k);
pr_0 = z(1:3);
pv_0 = z(4:6);
end
function out1 = Ca0fromIK(i, k)
  out1 = k + i.*(i.*k.*(i.*(0.75708324 - 0.22989681*k).*exp(k) + 1.637357) + 0.031784486);
end
function out1 = Ca1fromIK(i, k)
  out1 = i.^2.*(k + 2.518653).*(0.050697707*k.*(-6.1506586*i + k) - 0.8663172);
end
function out1 = Ci0fromK(x0)

  out1 = 0.415541771668952*x0.*exp(-1.0*x0).*log(x0).^2 + 0.25060964*log(x0) + 0.0796329*exp(-5.9961753*x0).*log(x0).^2;

end
function out1 = Ci1fromK(x0)

  out1 = 0.16930881*x0.*(x0 - 1.1203467).*(x0 - 0.8851997).*exp(x0).*log(x0)./(x0.^2.*(x0 - 0.8851997).*exp(x0) + 1.3398796*exp(x0) + log(2*x0));

end
