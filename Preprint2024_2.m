% в этом скрипте изучаем пространственные перелёты с помощью афинных
% преобразований
Preprint2024_1
%% перебираем разные наклонения
zf = inv(dxdX)'*Z0_2;
i_target_range=0:pi/32:pi;
i_best=zeros(length(i_target_range),1,1);
R_i_matrices=zeros(length(i_target_range),1,3,3);
%%
for j = 7:length(i_target_range)
    disp(j)
    i_target = i_target_range(j);
    O_target = i_target_range(j);
    i_0 = i_best(j-1,1,1);

    options_fsolve = optimoptions('fsolve','Display','off');
    fsolve_fun=@(i_matrix)fsolve_find_i(i_matrix,i_target, start_pos,start_vel,zf,t0,t_end,O_target);
    i_f = fsolve(fsolve_fun, i_0, options_fsolve);

    %изменение наклонения
    delta_i = i_f;
    n_rot_i = [1; 0; 0];
    q_rot_i = [cos(delta_i/2); n_rot_i*sin(delta_i/2)];
    R_i0 = quat2rotm(q_rot_i');
    R_i0(1:3,3)=R_i0(1:3,3)*R_i0(2,2);
    R_i0(1:3,2)=R_i0(1:3,2)*R_i0(2,2);
    %восходящий узел
    O_target = 0*pi/6;
    delta_O = O_target-pi;
    n_rot = [0; 0; 1];
    q_rot = [cos(delta_O/2); n_rot*sin(delta_O/2)];
    R_rot = quat2rotm(q_rot');
    R_iO=R_rot*R_i0*R_rot';
    [i,O,J,AN] = integrate_iO(R_iO, start_pos,start_vel,zf,t0,t_end);
    i_best(j,1,1) = i_f;
    R_i_matrices(j,1,:,:) = R_i0;
end
%%
load('mat-files/Preprint2024_2.mat')
i_coef_range=0:pi/180:pi/3;
O_coef_range=0:pi/30:2*pi;
for j = 1:length(i_coef_range)
    disp(j)
    for k = 1:length(O_coef_range)
    i_coef = i_coef_range(j);
    O_coef = O_coef_range(k);

    %изменение наклонения
    delta_i = atan(i_coef);
    n_rot_i = [1; 0; 0];
    q_rot_i = [cos(delta_i/2); n_rot_i*sin(delta_i/2)];
    R_i0 = quat2rotm(q_rot_i');
    R_i0(1:3,3)=R_i0(1:3,3)*R_i0(2,2);
    R_i0(1:3,2)=R_i0(1:3,2)*R_i0(2,2);
    %восходящий узел
    delta_O = O_coef-pi;
    n_rot = [0; 0; 1];
    q_rot = [cos(delta_O/2); n_rot*sin(delta_O/2)];
    R_rot = quat2rotm(q_rot');
    R_iO=R_rot*R_i0*R_rot';
    [i,O,p,a,J,AN] = integrate_iO(R_iO, start_pos,start_vel,zf,t0,t_end);
    end
end
%%
plot(i_target_range,i_best)
%%
zf = inv(dxdX)'*Z0_2;

%пространственный поворот
n_rot_i = [1; 0; 0];
%изменение наклонения
i_target = 20*pi/180;
delta_i = atan(7.06192346456835*i_target);
%delta_i = i_f;

q_rot_i = [cos(delta_i/2); n_rot_i*sin(delta_i/2)];
R_i = quat2rotm(q_rot_i');
R_i(1:3,3)=R_i(1:3,3)*R_i(2,2);
R_i(1:3,2)=R_i(1:3,2)*R_i(2,2);
n_rot = [0; 0; 1];
%восходящий узел
O_target = 2*pi/6;
delta_O = O_target-pi;
q_rot = [cos(delta_O/2); n_rot*sin(delta_O/2)];
R_rot = quat2rotm(q_rot');
R_iO=R_rot*R_i*R_rot';

pr = R_iO^(-1)*zf(1:3);
pv = R_iO^(-1)*zf(4:6);
y0 = cat(2,start_pos_2,start_vel_2,pr',pv',0,0)';
tspan = linspace(t0,t0+t_end, AN_i*400);
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
[t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);

J = y(end,13);       %функционал
AN_final = y(end,14);%угловая дальность
res_r_f = norm(y(end, 1:3)-mars_r_f');
res_v_f = norm(y(end, 4:6)-mars_v_f');
%% матрица частных производных
ddeltady0=zeros([6,6]);
step_h=sqrt(eps);
for i=7:12
    y0 = cat(2,start_pos_2,start_vel_2,pr',pv')';
    y0_delta=zeros([12,1]);
    y0_delta(i)=step_h;
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0+y0_delta,options);

    p_plus=y(end,1:6);

    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0-y0_delta,options);
    p_minus=y(end,1:6);
    partial=(p_plus-p_minus)/(2*step_h);
    ddeltady0(:,i-6)=partial;
end
%%
zf = inv(dxdX)'*Z0_2;

%эксцентриситет
ix = 0.0;
iy = 0;
L  = 0;
ey = 0;
ex = 0.01;
p = 1-ex^2;
%p = 1;
X_2 = equitoctial2decart([p;ex;ey;ix;iy;L], mug);
start_pos_3 = X_2(1:3)';
start_vel_3 = X_2(4:6)';
R_e = calculateRotMatrix(start_pos_3,start_vel_3);

%аргумент перицентра
o_target = 2*pi/6;
delta_o = o_target/2-pi/2;
q_rot = [cos(delta_o/2); n_rot*sin(delta_o/2)];
R_rot = quat2rotm(q_rot');
R_eo=R_rot*R_e*R_rot';

pr = R_eo^(-1)*zf(1:3);
pv = R_eo^(-1)*zf(4:6);
y0 = cat(2,start_pos_2,start_vel_2,pr',pv',0,0)';
tspan = linspace(t0,t0+t_end, AN_i*400);
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
[t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);

J = y(end,13);       %функционал
AN_final = y(end,14);%угловая дальность
res_r_f = norm(y(end, 1:3)-mars_r_f');
res_v_f = norm(y(end, 4:6)-mars_v_f');
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


%% выводим график
figure(5);
plot3(0, 0, 0, 'k--o');
set(gca,'FontSize',14);
hold on;

plot3(y(1, 1), y(1, 2), y(1, 3), 'O', 'LineWidth', 1);
plot3(y(:, 1), y(:, 2), y(:, 3), 'cyan', 'LineWidth', 1);
plot3(orbit_initial(:, 1), orbit_initial(:, 2), orbit_initial(:, 3), 'k');
plot3(1.52*[cos(delta_O+pi),0*cos(delta_O)],1.52*[sin(delta_O+pi),0*sin(delta_O)],[0,0])
%plot3(orbit_final(:, 1), orbit_final(:, 2), orbit_final(:, 3), 'r');
%plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k')
plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'k')
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(y(end, 1:3)',y(end, 4:6)',1);
[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
disp(['p=', num2str(p)])
disp(['O=', num2str(O)])
disp(['o=', num2str(lonPer)])
disp(['ex=', num2str(ex)])
disp(['ey=', num2str(ey)])
disp('-------------')
%disp(['i=', num2str(i)])
%disp(['ix=', num2str(ix)])
%disp(['iy=', num2str(iy)])
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

J = y(end,13);       %функционал
AN = y(end,14);%угловая дальность
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(y(end, 1:3)',y(end, 4:6)',1);
[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
end