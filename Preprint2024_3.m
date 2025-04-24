% в этом скрипте изучаем пространственные перелёты с помощью афинных
% преобразований
clear;
load('mat-files/Preprint2024_1_Mars.mat')
load('mat-files/Preprint2024_1_Venus.mat', 't_end_Venus')
%%
savefilename = 'mat-files/Preprint2024_3.mat';
k_coef_range=0.4:0.1:2.5;
i_coef_range=0:0.1:3;
O_coef_range=0:pi/15:2*pi;
N_count_1 = length(k_coef_range);
N_count_2 = length(i_coef_range);
N_count_3 = length(O_coef_range);

i_massive = zeros(N_count_1,N_count_2,N_count_3);
O_massive = zeros(N_count_1,N_count_2,N_count_3);
p_massive = zeros(N_count_1,N_count_2,N_count_3);
a_massive = zeros(N_count_1,N_count_2,N_count_3);
load(savefilename)
for j1 = 1:N_count_1
    for j2 = 1:N_count_2
        for j3 = 1:N_count_3
        k_coef = k_coef_range(j1);
        i_coef = i_coef_range(j2);
        O_coef = O_coef_range(j3);
        
        if p_massive(j1,j2,j3) > 0
            continue
        end
        disp([num2str(j1),'   ',num2str(j2),'   ',num2str(j3)]);
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
        t_end = get_t_end(k_coef);
        [pr_0,pv_0] = get_adjoint(k_coef);
        z = [pr_0;pv_0];
        [i,O,p,a] = integrate_iO_equinoctial(R_iO, start_pos,start_vel,z,t0,t_end);
        i_massive(j1,j2,j3) = i;
        O_massive(j1,j2,j3) = O;
        p_massive(j1,j2,j3) = p;
        a_massive(j1,j2,j3) = a;
        save(savefilename)
        end
    end
end
%%
k = 0.72;

t_end = get_t_end(k);
[pr_0,pv_0] = get_adjoint(k);
t0 = 0;
start_pos = [1 0 0];
start_vel = [0 1 0];
y0 = cat(2,start_pos,start_vel,pr_0',pv_0',0,0)';
tspan = linspace(t0,t0+t_end, round(t_end*20));
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
[t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);

J = y(end,13);       %функционал
AN_final = y(end,14);%угловая дальность
%рисуем
t_orbit = linspace(t0,t0+T_earth/(24*3600), 1000);

st.t = t_orbit';
st.planet = 'Earth';
st.mode = orbits;
st.delta_omega = 0;

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
%plot3(1.52*[cos(delta_O+pi),0*cos(delta_O)],1.52*[sin(delta_O+pi),0*sin(delta_O)],[0,0])
%plot3(orbit_final(:, 1), orbit_final(:, 2), orbit_final(:, 3), 'r');
%plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k')
plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'k')
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(y(end, 1:3)',y(end, 4:6)',1);
[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
disp(['p=', num2str(p)])
%disp(['O=', num2str(O)])
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

%J = y(end,13);       %функционал
%AN = y(end,14);%угловая дальность
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(y(end, 1:3)',y(end, 4:6)',1);
[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
end
function [i,O,p,a] = integrate_iO_equinoctial(R_iO, start_pos,start_vel,zf,t0,t_end)
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
load('mat-files/Preprint2024_1_Mars.mat','t_end_Mars')
load('mat-files/Preprint2024_1_Venus.mat','t_end_Venus')
k1=0.72;
k2=1.52;
CT1 = ((t_end_Venus-40*2*pi)/log(k1)-(t_end_Mars-40*2*pi)/log(k2))/(k1-k2);
CT0 = (t_end_Venus-40*2*pi)/log(k1)-CT1*k1;
t_end = 40*2*pi+(CT1*k+CT0)*log(k);
end
function [pr_0,pv_0] = get_adjoint(k)
%Adjoint variables approximation
% a - angular distance, revolutions
% k - radii ratio
% d - multiplier of D coefficients, 0 or 1
%MARS
load('mat-files/Preprint2024_1_Mars.mat','Z0_Mars','dxdX')
load('mat-files/Preprint2024_1_Venus.mat','Z0_Venus')
z_Mars = inv(dxdX)'*Z0_Mars;
z_Venus = inv(dxdX)'*Z0_Venus;
pr_x_M = z_Mars(1);
pr_x_V = z_Venus(1);
pr_y_M = z_Mars(2);
pr_y_V = z_Venus(2);

pv_x_M = z_Mars(4);
pv_x_V = z_Venus(4);
pv_y_M = z_Mars(5);
pv_y_V = z_Venus(5);

k1=0.72;
k2=1.52;

%prx
p1 = pr_x_V;
p2 = pr_x_M;
g = (p2*power(k2,1/4)/log(k2) - p1*k1*power(k2,-3/4)/log(k1))/(1-power(k1/k2,3/4));
f = p1*k1/log(k1)-g*power(k1,3/4);
pr_x = @(k)(log(k)*(f/k+g*power(k,-1/4)));
%pry
p1 = pr_y_V;
p2 = pr_y_M;
g = (p2*k2^(3/2)/log(k2) - p1*k1*sqrt(k2)/log(k1))/(1-sqrt(k2/k1));
f = p1*k1/log(k1)-g/sqrt(k1);
pr_y = @(k)(log(k)*(f/k+g/(k^(3/2))));
%pvx
p1 = pv_x_V;
p2 = pv_x_M;
g = (p2*k2^(3/2)/log(k2) - p1*k1*sqrt(k2)/log(k1))/(1-sqrt(k2/k1));
f = p1*k1/log(k1)-g/sqrt(k1);
pv_x = @(k)(log(k)*(f/k+g/(k^(3/2))));
%pvy
p1 = pv_y_V;
p2 = pv_y_M;
g = (p2*power(k2,1/4)/log(k2) - p1*k1*power(k2,-3/4)/log(k1))/(1-power(k1/k2,3/4));
f = p1*k1/log(k1)-g*power(k1,3/4);
pv_y = @(k)(log(k)*(f/k+g*power(k,-1/4)));

pr_0 = [pr_x(k);pr_y(k);0.];
pv_0 = [pv_x(k);pv_y(k);0.]; 
end