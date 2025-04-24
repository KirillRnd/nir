%внешние константы
clc;
clear;
ae = 149597870700;
mug_0 = 132712.43994*(10^6)*(10^(3*3));
mug=1;
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;
T_mars_days = 365.256363004*1.8808476;

r_unit=ae;
V_unit=sqrt(mug_0/ae);
T_unit = T_earth/(2*pi);
%гиперпараметры
t_start=0;
orbits='Simple';
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

B=0.2721831;
A=0.5433279;
T_a  = @(a, a_rel)(a_rel-a_rel*B*log(a_rel))*a+(a_rel*B^2*log(a_rel));
a_OM = @(OM, a_rel,k)(OM+k+A*exp(-a_rel)*log(a_rel)*log(a_rel))/(A*(1+exp(-a_rel))*log(a_rel));
OM_a = @(a, a_rel)A*(a*(exp(a_rel)+1)-log(a_rel))*exp(-a_rel)*log(a_rel);

%%
%интегрируем
a_max=0.04;
tau=1;

t0 = t_start;
AN_i = 1;
d_coef = 1;
a_rel = 1.52;
T_i = 2*pi*T_a(AN_i, a_rel);
OM_i = 2*pi*OM_a(AN_i, a_rel);
start_vel_2 = start_vel;
start_pos_2 = start_pos;
%start_vel_2(1) = 0.4;
%start_pos_2(1) = 1.2;
%start_pos_2(2) = 0.2;
%start_vel_2(2) = 0.7;
%пространственный поворот
n_rot_i = [1; 0; 0];
delta_i = 0*pi/180;
q_rot_i = [cos(delta_i/2); n_rot_i*sin(delta_i/2)];
R_1 = quat2rotm(q_rot_i');
R_1(1:3,3)=R_1(1:3,3)*R_1(2,2);
R_1(1:3,2)=R_1(1:3,2)*R_1(2,2);
%R_1(1:3,1)=R_1(1:3,1)*R_1(1,1);
n_rot = [0; 0; 1];
delta_i = 0*pi/6;
q_rot = [cos(delta_i/2); n_rot*sin(delta_i/2)];
R_rot = quat2rotm(q_rot');
R_1=R_rot*R_1*R_rot';
 
ix = 0.0;
iy = 0;
L  = 0;
ey = 0;
ex = -0.0;
p = 1-ex^2;
X_2 = equitoctial2decart([p;ex;ey;ix;iy;L], mug);
start_pos_3 = X_2(1:3)';
start_vel_3 = X_2(4:6)';
R_2 = calculateRotMatrix(start_pos_3,start_vel_3);
delta_phi = 0*pi/24;
R_rot = [cos(delta_phi),-sin(delta_phi),0;
        sin(delta_phi),cos(delta_phi),0;
        0,0,1];
R_2=R_rot*R_2*R_rot';
R_2=R_2^(-1);

R = calculateRotMatrix(start_pos_2,start_vel_2);

[pr_0,pv_0] = get_initial_adjoint(AN_i,a_rel,d_coef);
pr_0 = (R_1^(-1)*R^(-1)*pr_0)';
pv_0 = (R_1^(-1)*R^(-1)*pv_0)';
y0 = cat(2,start_pos_2,start_vel_2,pr_0,pv_0,0,0,0)';

st.t = [t_start, T_i*T_unit/(3600*24)];
st.planet = planet_end;
st.mode = orbits;
st.delta_omega = OM_i;

[mars_r_f, mars_v_f]=planetModel(st);
mars_r_f_target=mars_r_f'*1e+03/r_unit;
mars_v_f_target=mars_v_f'*1e+03/V_unit;
yf = [mars_r_f_target;mars_v_f_target]';

fmincon_traj_fun=@(z)fmincon_traj(z,y0,yf,t0,T_i, a_max, tau);
options = optimoptions('fmincon','Display','iter');
options = optimoptions(options,'UseParallel', true);
options = optimoptions(options, 'OptimalityTolerance', 1e-10);
options = optimoptions(options, 'MaxFunctionEvaluations', 1e+10);
options = optimoptions(options, 'StepTolerance', 1e-10);
options = optimoptions(options, 'ConstraintTolerance', 1e-10);
options = optimoptions(options, 'MaxIterations', 400);
%options = optimoptions(options, 'FiniteDifferenceType', 'central');
options = optimoptions(options, 'EnableFeasibilityMode', true);

%pr_0 = [0.9862, -0.0107, 0];
%pv_0 = [-0.0107, 0.9850, 0];
pr_0 = pr_0/norm(pr_0);
pv_0 = pv_0/norm(pv_0);
p0 = [pr_0,pv_0]';
p_limit = 10;
A = [];
b = [];
Aeq = [];
beq = [];
lb = p0-[1; 1; 0; 1; 1; 0]*p_limit;
ub = p0+[1; 1; 0; 1; 1; 0]*p_limit;
pf = fmincon(@(x)1, p0, A, b, Aeq, beq, lb, ub, fmincon_traj_fun, options);
y0(7:12)=pf;
%%
tspan = linspace(t0,t0+T_i, AN_i*1000);
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
[t,y] = ode15s(@(t,y) HomotopyIntegration3D_withJ_TH(t,y,a_max, tau), tspan,y0,options);

J = y(end,13);       %функционал квадратичный
AN_final = y(end,14);%угловая дальность
J_linear = y(:,15);  %функционал линейный

st.t = [t_start, T_i*T_unit/(3600*24)];
%st.t = [t_start, t(end)*T_unit/(3600*24)];
st.planet = planet_end;
st.mode = orbits;
st.delta_omega = OM_i;

[mars_r_f, mars_v_f]=planetModel(st);
mars_r_f=mars_r_f'*1e+03/r_unit;
mars_v_f=mars_v_f'*1e+03/V_unit;
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

%% break points
diffs = abs(diff(J_linear))/max(diff(J_linear));
% Порог для определения точки разрыва
threshold = 1e-5;
threshold2 = 1e-4;
% Индексы, где разность превышает порог
nonzero_derivative = diffs > threshold;
nonzero_derivative(end+1) = nonzero_derivative(end);

diffs2 = abs(diff(nonzero_derivative));

% Индексы, где разность превышает порог
break_points = find(diffs2 > threshold2);
%выводим график
figure(5);
plot3(0, 0, 0, 'k--o', 'HandleVisibility','off');
set(gca,'FontSize',14);
hold on;

plot3(y(1, 1), y(1, 2), y(1, 3), 'O', 'LineWidth', 1, 'HandleVisibility','off');
plot3(y(:, 1), y(:, 2), y(:, 3), 'cyan', 'LineWidth', 1, 'HandleVisibility','off');
plot3(orbit_initial(:, 1), orbit_initial(:, 2), orbit_initial(:, 3), 'k--', 'DisplayName','Орбита Земли');
%plot3(orbit_final(:, 1), orbit_final(:, 2), orbit_final(:, 3), 'r');
%plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k')
plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'k', 'DisplayName','Орбита Марса')
%plot3(1.52*[cos(delta_i+pi),cos(delta_i)],1.52*[sin(delta_i+pi),sin(delta_i)],[0,0])
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(y(end, 1:3)',y(end, 4:6)',1);
% disp(['e=', num2str(eMag)])
% disp(['a=', num2str(a)])
% disp(['i=', num2str(i)])
plot3(y(end, 1), y(end, 2), y(end, 3),'bO', 'HandleVisibility','off')
plot3(y(1,1),y(1,2),y(1,3), 'magenta', 'LineWidth',2, 'DisplayName','Активный участок')
plot3(y(1,1),y(1,2),y(1,3), 'green', 'LineWidth',2, 'DisplayName','Пассивный участок')

plot3(mars_r_f_target(1), mars_r_f_target(2), mars_r_f_target(3),'rO', 'HandleVisibility','off')


color_count = nonzero_derivative(1);
for i = 1:length(break_points) + 1
    if color_count == 0
        color = 'green';
        disp_name = 'Пассивный участок';
    else
        color = 'magenta';
        disp_name = 'Активный участок';
    end
    color_count = ~color_count;
    if length(break_points) == 0
        plot3(y(1:end,1), y(1:end,2), y(1:end,3), color, 'LineWidth', 2);
    else
        if i == 1
            % Первая часть до первого разрыва
            plot3(y(1:break_points(i),1), y(1:break_points(i),2), y(1:break_points(i),3), color, 'LineWidth', 2, 'HandleVisibility','off');
        elseif i == length(break_points) + 1
            % Последняя часть после последнего разрыва
            plot3(y(break_points(i-1)+1:end, 1), y(break_points(i-1)+1:end, 2), y(break_points(i-1)+1:end, 3), color, 'LineWidth', 2, 'HandleVisibility','off');
        else
            % Промежуточные части
            plot3(y(break_points(i-1)+1:break_points(i), 1), y(break_points(i-1)+1:break_points(i), 2),y(break_points(i-1)+1:break_points(i), 3), color, 'LineWidth', 2, 'HandleVisibility','off');
        end
    end
    
end
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
legend;
%%
function [c, ceq] = fmincon_traj(z,y0,yf,t0,t_end, a_max, tau)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tspan = linspace(t0,t0+t_end, 10);
y0(7:12)=z;
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);   
[~,y] = ode15s(@(t,y) HomotopyIntegration3D_withJ_TH(t,y, a_max, tau), tspan,y0,options);
c =[];
ceq = (y(end,1:6)-yf(1:6));
end
function dydt = HomotopyIntegration3D_withJ_TH(~,y, a_max, tau)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
mug=1;
%первый набор переменных
r=y(1:3);
v=y(4:6);
pr=y(7:9);
pv=y(10:12);
nr=norm(r);

npv = norm(pv);
if npv < 1e-12
    e_star = [0; 0; 0];
else
    e_star = pv / npv;
end

S = dot(pv, e_star)-1;
epsTol = 1e-8;
if S > epsTol
    a_star = a_max;
elseif S < -epsTol
    a_star = 0;
else
    % Сингулярный или "пограничный" случай S=0
    a_star = 0;
end

dydt = zeros(size(y));
dydt(1:3)=v;
dydt(4:6)=-mug*r/nr^3+(1-tau)*pv+tau*a_star*e_star;
dydt(7:9)=mug*pv/nr^3 - mug*(3*(r*r')/nr^(5))*pv;
dydt(10:12)=-pr;

dydt(13)=(norm(pv)^2)/2; %dJdt
dydt(14)=(r(1)*v(2) - r(2) * v(1)) / (norm(r)^2); %dTHdt
dydt(15)=a_star; %dJdt
end