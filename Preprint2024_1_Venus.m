%% В этом скрипте получаем ршение задачи встречи для 40 витков
% делал для препринта весной 2024
%внешние константы
%clc;
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
planet_end = 'Venus';
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

%интегрируем
t0 = t_start;
d_coef = 0;
a_rel = 0.72;
AN_i = 40;
T_i = 2*pi*TfromAN(AN_i, a_rel);
OM_i = 2*pi*OMfromAN(AN_i, a_rel);
start_vel_2 = start_vel;
start_pos_2 = start_pos;

[pr_0,pv_0] = get_initial_adjoint(AN_i,a_rel,d_coef);
pr_0 = (pr_0)';
pv_0 = (pv_0)';
t_end_0 = T_i-1.15;
y0 = [start_pos,start_vel,pr_0,pv_0];
z0 = [pr_0,pv_0];

[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(start_pos_2',start_vel_2',mug);

[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
X0 = [p;ex;ey;ix;iy;L];
dxdX = equitoctial2decart_jacobian(X0,mug);
Z0 = dxdX'*[pr_0';pv_0'];
Y0 = [X0;Z0];

mars_r_f = 0.72*[cos(AN_i*2*pi);sin(AN_i*2*pi);0];
mars_v_f = (1/sqrt(0.72))*[cos(AN_i*2*pi+pi/2);sin(AN_i*2*pi+pi/2);0];

yf = [mars_r_f;mars_v_f]';

[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(mars_r_f,mars_v_f,mug);
[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
YF = [p,ex,ey,ix,iy,L];

minimize_delta_t_end = @(delta_t_end) find_close_solution(delta_t_end, t_end_0, y0, Z0, t0,AN_i);
angle_range = pi/2;
[delta_t_end, Jt_end_min] = fminbnd(minimize_delta_t_end, -angle_range, angle_range);
t_end = delta_t_end+t_end_0;
%t_end = t_end_0;

options_fsolve = optimoptions('fsolve','Display','off');
fsolve_traj_fun=@(z)fsolve_traj(z,Y0,yf,t0,t_end);
ZF = fsolve(fsolve_traj_fun, Z0, options_fsolve);
zf = inv(dxdX)'*ZF;
%zf = z0';
y0 = cat(2,start_pos_2,start_vel_2,zf(1:3)',zf(4:6)',0,0)';
tspan = linspace(t0,t0+t_end, AN_i*400);
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
[t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);

J = y(end,13);       %функционал
AN_final = y(end,14);%угловая дальность
res_r_f = norm(y(end, 1:3)-mars_r_f');
res_v_f = norm(y(end, 4:6)-mars_v_f');
Z0_Venus = ZF;
t_end_Venus = t_end;
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
disp(['J=', num2str(J)])
%disp(['i=', num2str(i)])
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
function Jt_end = find_close_solution(delta_t_end, t_end_0, y0, Z0, t0,AN_i)
    
    mars_r_f = 0.72*[cos(AN_i*2*pi);sin(AN_i*2*pi);0];
    mars_v_f = (1/sqrt(0.72))*[cos(AN_i*2*pi+pi/2);sin(AN_i*2*pi+pi/2);0];
    t_end = t_end_0+delta_t_end;
    start_pos = y0(1:3);
    start_vel = y0(4:6);
    mug = 1;
    [~,eMag,i,O,~,nu,~,~,lonPer,p] = rv2orb(start_pos',start_vel',mug);

    [p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
    X0 = [p;ex;ey;ix;iy;L];
    dxdX = equitoctial2decart_jacobian(X0,mug);
    Y0 = [X0;Z0];

    yf = [mars_r_f;mars_v_f]';
    options_fsolve = optimoptions('fsolve','Display','off');
    fsolve_traj_fun=@(z)fsolve_traj(z,Y0,yf,t0,t_end);
    ZF = fsolve(fsolve_traj_fun, Z0, options_fsolve);
    Y0(7:12)=ZF;
    Y0(13)=0;
    % zf = inv(dxdX)'*ZF;
    % 
    % pr0=zf(1:3);
    % pv0=zf(4:6);
    % y0(7:9)=pr0;
    % y0(10:12)=pv0;
    % y0(13:14)=0;
    tspan = linspace(t0,t0+t_end, round(t_end*10));
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    % [t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);
    % 
    % Jt_end = y(end, 13);
    [t,y] = ode113(@(t,y) equitoctial_adjoint_for_integration_withJ(y,mug), tspan,Y0,options);

    traj = arrayfun(@(p,ex,ey,ix,iy,L) equitoctial2decart([p;ex;ey;ix;iy;L], mug),...
    y(:, 1),y(:, 2),y(:, 3),y(:, 4),y(:, 5),y(:, 6),'UniformOutput',false);
    traj = cell2mat(traj')';
    res = traj(end,1:6)-yf(1:6);
    Jt_end = y(end, 13);
    
end
function res = fsolve_traj(z,y0,yf,t0,t_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
y0(7:12)=z;
tspan = linspace(t0,t0+t_end, 10);
mug = 1;
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   

[t,y] = ode113(@(t,y) equitoctial_adjoint_for_integration(y,mug), tspan,y0,options);

traj = arrayfun(@(p,ex,ey,ix,iy,L) equitoctial2decart([p;ex;ey;ix;iy;L], mug),...
    y(:, 1),y(:, 2),y(:, 3),y(:, 4),y(:, 5),y(:, 6),'UniformOutput',false);
traj = cell2mat(traj')';
res = (traj(end,1:6)-yf(1:6));
end
function [c, ceq] = fun2min_border(z,param)
y0 = param.y0;
yf = param.yf;
t0 = param.t0;
t_end_0 = param.t_end_0;
t_end = t_end_0+z(7);
y0(7:12)=param.p_scaler*z(1:6);
tspan = linspace(t0,t0+t_end, 10);
mug = 1;
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   

[~,y] = ode113(@(t,y) equitoctial_adjoint_for_integration(y,mug), tspan,y0,options);

traj = arrayfun(@(p,ex,ey,ix,iy,L) equitoctial2decart([p;ex;ey;ix;iy;L], mug),...
    y(:, 1),y(:, 2),y(:, 3),y(:, 4),y(:, 5),y(:, 6),'UniformOutput',false);
traj = cell2mat(traj')';
c =[];
ceq = (traj(end,1:6)-yf(1:6));
end
function J  = fun2min_J(z,param)
y0 = param.y0;
t0 = param.t0;
t_end_0 = param.t_end_0;
t_end = t_end_0+z(7);
y0(7:12)=param.p_scaler*z(1:6);
y0(13)=0;
tspan = linspace(t0,t0+t_end, 10);
mug = 1;
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   

[~,y] = ode113(@(t,y) equitoctial_adjoint_for_integration_withJ(y,mug), tspan,y0,options);

J = y(end,13);
end