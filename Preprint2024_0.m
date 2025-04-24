%% В этом скрипте сравниваем задачу встречи с задачей перелёта
% делал для препринта весной 2024
%внешние константы
clc;
clear;
ae = 149597870700;
mug_0 = 132712.43994*(10^6)*(10^(3*3));
mug=1;
T_earth_days = 365.256363004;
T_earth = T_earth_days*3600*24;
T_mars_days = T_earth_days*1.8808476;
T_mars=T_earth*1.8808476;

r_unit=ae;
V_unit=sqrt(mug_0/ae);
T_unit = T_earth/(2*pi);
%гиперпараметры
t_start=0;
orbits='Flat';
N=1350;
m0=367;
eta=0.45;
planet_start = 'Earth';
planet_end = 'Mars';
omega = -pi;

%определяем стартовое положение
st.t = t_start;
st.planet = planet_start;
st.mode = orbits;
st.delta_omega = omega;

[start_pos, start_vel] = planetModel(st);
start_pos=start_pos*1e+03/ae;
start_vel=start_vel*1e+03/V_unit;


%%
T_range = 1:0.1:11;
N_count = length(T_range) ;
PRevery_transfer=zeros(N_count, 3); %базис-вектор
PVevery_transfer=zeros(N_count, 3); %базис-вектор
Jevery_transfer=zeros(N_count, 1);     %функционал
Tevery_transfer = zeros(N_count, 1);   %время перелёта
ANevery_transfer = zeros(N_count, 1);  %угловая дальность
OMevery_transfer = zeros(N_count, 1);  %разность фаз

%%
savefilename = 'mat-files/Preprint2024_0.mat';

load(savefilename)
for i=1:N_count
    if Jevery_transfer(i)>0
        continue
    end
    disp(i)
    %интегрируем
    t0 = t_start;
    T_i = 2*pi*T_range(i);
    d_coef = 0;
    a_rel = 1.52;
    AN_i = ANfromT(T_i/(2*pi), a_rel);
    %T_i = 2*pi*TfromAN(AN_i, a_rel);
    OM_i = 2*pi*OMfromAN(AN_i, a_rel);
    start_vel_2 = start_vel;
    start_pos_2 = start_pos;
    
    [pr_0,pv_0] = get_initial_adjoint(AN_i,a_rel,d_coef);
    pr_0 = (pr_0)';
    pv_0 = (pv_0)';
    t_end = T_i;
    y0 = [start_pos,start_vel,pr_0,pv_0];
    z0 = [pr_0,pv_0];
    
    minimize_delta_omega = @(delta_omega) find_close_solution(delta_omega, OM_i, y0, z0, t0,t_end);
    angle_range = pi/8;
    [delta_omega, Jt_end_min] = fminbnd(minimize_delta_omega, -angle_range, angle_range);
    omega_min = delta_omega+OM_i;
    st.t = [t0, T_earth_days*t_end/(2*pi)];
    st.planet = 'Mars';
    st.mode = 'Flat';
    st.delta_omega = omega_min;
    
    [mars_r_f, mars_v_f]=planetModel(st);
    mars_r_f=mars_r_f'*1e+03/r_unit;
    mars_v_f=mars_v_f'*1e+03/V_unit;
    
    yf = [mars_r_f;mars_v_f]';
    
    options_fsolve = optimoptions('fsolve','Display','off');
    fsolve_traj_fun=@(z)fsolve_traj(z,y0,yf,t0,t_end);
    zf = fsolve(fsolve_traj_fun, z0, options_fsolve);
    
    y0 = cat(2,start_pos_2,start_vel_2,pr_0,pv_0,0,0)';
    y0(7:9)=zf(1:3);
    y0(10:12)=zf(4:6);
    y0(13:14)=0;
    tspan = linspace(t0,t0+T_i, AN_i*400);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);
    
    J = y(end,13);       %функционал
    AN_final = y(end,14);%угловая дальность
    res_r_f = norm(y(end, 1:3)-mars_r_f');
    res_v_f = norm(y(end, 4:6)-mars_v_f');

    PRevery_transfer(i,:)=zf(1:3); %базис-вектор
    PVevery_transfer(i,:)=zf(4:6); %базис-вектор
    Jevery_transfer(i)=J;     %функционал
    Tevery_transfer(i) = t_end;   %время перелёта
    ANevery_transfer(i) = AN_final;  %угловая дальность
    OMevery_transfer(i) = omega_min;  %разность фаз
    save(savefilename);
end
%% второй проход
savefilename = 'mat-files/Preprint2024_1.mat';

load(savefilename)
%%
for i=1:50
    disp(i)
    %интегрируем
    t0 = t_start;
    T_i = 2*pi*T_range(i);
    d_coef = 0;
    a_rel = 1.52;
    AN_i = ANfromT(T_i/(2*pi), a_rel);
    %T_i = 2*pi*TfromAN(AN_i, a_rel);
    OM_i = OMevery_transfer(i);
    start_vel_2 = start_vel;
    start_pos_2 = start_pos;
    
    [pr_0,pv_0] = get_initial_adjoint(AN_i,a_rel,d_coef);
    pr_0 = (pr_0)';
    pv_0 = (pv_0)';
    pr_0 = PRevery_transfer(i,1:3);
    pv_0 = PVevery_transfer(i,1:3);
    t_end = T_i;
    y0 = [start_pos,start_vel,pr_0,pv_0];
    z0 = [pr_0,pv_0];
    
    minimize_delta_omega = @(delta_omega) find_close_solution(delta_omega, OM_i, y0, z0, t0,t_end);
    angle_range = pi/8;
    [delta_omega, Jt_end_min] = fminbnd(minimize_delta_omega, -angle_range, angle_range);
    omega_min = delta_omega+OM_i;
    st.t = [t0, T_earth_days*t_end/(2*pi)];
    st.planet = 'Mars';
    st.mode = 'Flat';
    st.delta_omega = omega_min;
    
    [mars_r_f, mars_v_f]=planetModel(st);
    mars_r_f=mars_r_f'*1e+03/r_unit;
    mars_v_f=mars_v_f'*1e+03/V_unit;
    
    yf = [mars_r_f;mars_v_f]';
    
    options_fsolve = optimoptions('fsolve','Display','off');
    fsolve_traj_fun=@(z)fsolve_traj(z,y0,yf,t0,t_end);
    zf = fsolve(fsolve_traj_fun, z0, options_fsolve);
    
    y0 = cat(2,start_pos_2,start_vel_2,pr_0,pv_0,0,0)';
    y0(7:9)=zf(1:3);
    y0(10:12)=zf(4:6);
    y0(13:14)=0;
    tspan = linspace(t0,t0+T_i, AN_i*400);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);
    
    J = y(end,13);       %функционал
    AN_final = y(end,14);%угловая дальность
    res_r_f = norm(y(end, 1:3)-mars_r_f');
    res_v_f = norm(y(end, 4:6)-mars_v_f');

    PRevery_transfer(i,:)=zf(1:3); %базис-вектор
    PVevery_transfer(i,:)=zf(4:6); %базис-вектор
    Jevery_transfer(i)=J;     %функционал
    Tevery_transfer(i) = t_end;   %время перелёта
    ANevery_transfer(i) = AN_final;  %угловая дальность
    OMevery_transfer(i) = omega_min;  %разность фаз
    save(savefilename);
end
%%

load('D:\MATLAB\ipm\mat-files\16-Jun-2024.mat',"Tevery_a","Jevery_a","ANevery_a")
Tevery_a_MARS_true = Tevery_a(:,113);
Jevery_a_MARS_true = Jevery_a(:,113);
ANevery_a_MARS_true = ANevery_a(:,113);
J_unit = 176.62545106129272;
figure(2)
plot(Tevery_transfer/(2*pi),Jevery_transfer)
hold on;
plot(Tevery_a_MARS_true/T_earth_days,Jevery_a_MARS_true/J_unit)
hold off;
grid on;
figure(4)
plot(ANevery_transfer/(2*pi),Jevery_transfer)
hold on;
plot(ANevery_a_MARS_true/(2*pi),Jevery_a_MARS_true/J_unit)
hold off;
grid on;
%%
%рисуем
t0 = t_start;
t_orbit = linspace(t0,t0+T_earth/(24*3600), 1000);

st.t = t_orbit';
st.planet = 'Earth';
st.mode = orbits;
st.delta_omega = omega;

earth_traj = planetModel(st);
earth_traj=1e+3*earth_traj/ae;

t_orbit = linspace(t0,t0+T_mars/(24*3600), 1000);

st.t = t_orbit';
st.planet = planet_end;
st.mode = orbits;
st.delta_omega = 0;

mars_traj = planetModel(st);
mars_traj=1e+3*mars_traj/ae;
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
y0 = cat(2,y(1, 1:3),y(1, 4:6),[0,0,0],[0,0,0])';
[t,orbit_initial] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
y0 = cat(2,y(end, 1:3),y(end, 4:6),[0,0,0],[0,0,0])';
[t,orbit_final] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);


%выводим график
figure(5);
plot3(0, 0, 0, 'k--o');
set(gca,'FontSize',14);
hold on;

plot3(y(1, 1), y(1, 2), y(1, 3), 'O', 'LineWidth', 1);
plot3(y(:, 1), y(:, 2), y(:, 3), 'cyan', 'LineWidth', 1);
plot3(orbit_initial(:, 1), orbit_initial(:, 2), orbit_initial(:, 3), 'k');
%plot3(orbit_final(:, 1), orbit_final(:, 2), orbit_final(:, 3), 'r');
%plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k')
plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'k')
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(y(end, 1:3)',y(end, 4:6)',1);
disp(['e=', num2str(eMag)])
disp(['a=', num2str(a)])
disp(['i=', num2str(i)])
plot3(y(end, 1), y(end, 2), y(end, 3),'bO')

plot3(mars_r_f(1), mars_r_f(2), mars_r_f(3),'rO')

hold off;
axis equal

%title('Траектория КА')
xlabel('x, AU')
ylabel('y, AU')
zlabel('z, AU')
view(0,90)
%view(90,0)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
grid on;
box off;
%%
function Jt_end = find_close_solution(delta_omega, omega0, y0, z0, t0,t_end)
    mug_0 = 132712.43994*(10^6)*(10^(3*3));
    ae = 149597870700;
    r_unit=ae;
    V_unit=sqrt(mug_0/ae);
    T_earth_days = 365.256363004;
    omega_min = delta_omega+omega0;
    st.t = [t0, T_earth_days*t_end/(2*pi)];
    st.planet = 'Mars';
    st.mode = 'Flat';
    st.delta_omega = omega_min;
    
    [mars_r_f, mars_v_f]=planetModel(st);
    mars_r_f=mars_r_f'*1e+03/r_unit;
    mars_v_f=mars_v_f'*1e+03/V_unit;
    
    yf = [mars_r_f;mars_v_f]';
    options_fsolve = optimoptions('fsolve','Display','off');
    fsolve_traj_fun=@(z)fsolve_traj(z,y0,yf,t0,t_end);
    zf = fsolve(fsolve_traj_fun, z0, options_fsolve);
    
    pr0=zf(1:3);
    pv0=zf(4:6);
    y0(7:9)=pr0;
    y0(10:12)=pv0;
    y0(13:14)=0;
    tspan = linspace(t0,t0+t_end, round(t_end*100));
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);
    
    Jt_end = y(end, 13);
end
function res = fsolve_traj(z,y0,yf,t0,t_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tspan = linspace(t0,t0+t_end, 10);
pr0=z(1:3);
pv0=z(4:6);
y0(7:9)=pr0;
y0(10:12)=pv0;
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
[~,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
res = y(end,1:6)-yf(1:6);
end