%внешние константы
ae = 149597870700;
mug_0 = 132712.43994*(10^6)*(10^(3*3));
T_earth_days = 365.256363004;
T_earth = T_earth_days*3600*24;
T_mars=T_earth*1.8808476;
T_mars_days = T_earth_days*1.8808476;

%поворот кватернионом в плоскость эклиптики
n_eps = [1; 0; 0];
eps_0 = 2*pi*(23+26/60+21/3600)/360;
q_eps = [cos(eps_0/2); n_eps*sin(eps_0/2)];
%коэффициенты обезразмеривания
r_unit=ae;
V_unit=sqrt(mug_0/ae);
T_unit = T_earth/(2*pi);
%гиперпараметры
t_start = juliandate(2026,1,1);
orbits='Ephemeris';
N=1350;
m0=367;
eta=0.45;
planet_start = 'Earth';
planet_end = 'Mars';

%определяем стартовое положение
st.t = t_start;
st.planet = planet_start;
st.mode = orbits;

[start_pos, start_vel] = planetModel(st);
start_pos=start_pos*1e+03/ae;
start_vel=start_vel*1e+03/V_unit;
R = calculateRotMatrix(start_pos,start_vel);

%формулы аппроксимации
pr_x = @(a)(0.20074211987250948)/(7.780825627926492*a-0.34159993346898454+sin(2*pi*a+0.1231174475945867+0.5853586584957071/(a^2-1.5690074915306857*a+1.9819223815154108)) + (-41.26470185582318+86.35986725318767*a-51.08371071757919*a^2)*exp(-3.424188046763629*a) ); %OK
pr_y = @(a)(0.0015463527047873246 - 0.0031678016585249486 * cos(2*pi * a+(-2.7056930690628103*a^2+5.176300671261158*a-2.336378597152512)*exp(3.2003639491958467*a-3.2205106566099353*a^2)))/(a^2-0.18328217249773662*a+0.11514749689845788);%OK
pv_x = @(a)(0.0018661818762522732 - 0.0031537940448213170 * cos(2*pi * a+(-7.277086750486235*a^2+13.74021831427378*a-6.262552219285274)*exp(1.381195073394106*a-2.278956114228774*a^2)))/(a^2-0.19152227381595716*a+0.0993237695309777);
pv_y = @(a)(0.10456877009148134)/(4.050922375438426*a-0.15902042687671913+sin(2*pi*a+0.11429217276196703+0.012218372217935064/(a^2-2.8347986183365705*a+2.06664561825753)-0.0019877228343783706*a^2) +(-540.5981059117084*a^2+958.0032106019531*a-415.15730979365355)*exp(-5.877538308252977*a));

J_a  = @(a)(0.002117*exp(a))/(a*exp(a)+0.17856*(sin(7.44443*a)-a));

%% получаем решение для графика и строим график
t_end = 500/T_earth_days; %дни переводим в годы
a_rel=1.52; %соотношение полуосей Земли и Марса
B=0.2721831;
AN_mars = (t_end-a_rel*B^2*log(a_rel))/(a_rel-a_rel*B*log(a_rel)); %получаем оптимальную угловую дальность

%Thetha 0 
st.mode = orbits;
st.t = t_start;
st.planet = planet_end;

[mars_r_0, mars_v_0]=planetModel(st);
st.planet = planet_start;
[earth_r_0, earth_v_0]=planetModel(st);
mars_r_0 = quatrotate(q_eps',mars_r_0);
earth_r_0 = quatrotate(q_eps',earth_r_0);
%в проекции на эклиптику
mars_an_0 = atan2(mars_r_0(2), mars_r_0(1));
earth_an_0 = atan2(earth_r_0(2), earth_r_0(1));
delta_theta_0 = (earth_an_0-mars_an_0)/(2*pi);
%AN_mars = AN_mars-delta_theta_0*(a_rel-a_rel*B*log(a_rel));
AN_mars = AN_mars;


pr_0 = R^(-1)*[pr_x(AN_mars);pr_y(AN_mars);0];
pv_0 = R^(-1)*[pv_x(AN_mars);pv_y(AN_mars);0]; 
tspan = linspace(0,t_end*2*pi, 1000);

y0 = cat(2,start_pos,start_vel,pr_0',pv_0')';

options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10); 
options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 0.1, 10));
[t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);

st.mode = orbits;
st.t = t_start+t_end*T_earth_days;
st.planet = planet_end;

[mars_r_f, mars_v_f]=planetModel(st);
mars_r_f=mars_r_f'*1e+03/r_unit;
mars_v_f=mars_v_f'*1e+03/V_unit;

traj = y(:,1:3);
traj_v = y(:,4:6);



res_r_f = norm(y(end, 1:3)-mars_r_f');
res_v_f = norm(y(end, 4:6)-mars_v_f');
yf = [mars_r_f;mars_v_f]';

%ищем delta_theta_f
mars_r_f_new = quatrotate(q_eps',mars_r_f');
mars_an_f = atan2(mars_r_f_new(2), mars_r_f_new(1));
traj_new = quatrotate(q_eps',traj);
traj_an_f = atan2(traj_new(end,2), traj_new(end,1));
delta_theta_f = (traj_an_f-mars_an_f)/(2*pi);
if delta_theta_f > 0.5
    delta_theta_f=delta_theta_f-1;
elseif delta_theta_f < -0.5
    delta_theta_f=delta_theta_f+1;
end
AN_mars = AN_mars-delta_theta_f;

pr_0 = R^(-1)*[pr_x(AN_mars);pr_y(AN_mars);0];
pv_0 = R^(-1)*[pv_x(AN_mars);pv_y(AN_mars);0]; 
tspan = linspace(0,t_end*2*pi, 1000);
y0 = cat(2,start_pos,start_vel,pr_0',pv_0')';

opt_type = 'fsolve';
p0 = [pr_0;pv_0]';
if strcmp(opt_type,'fsolve')

    fsolve_traj_fun=@(z)fsolve_traj(z,y0,yf,0,t_end*2*pi);
    options = optimoptions('fsolve','Display','off');
    %options = optimoptions(options, 'Algorithm', 'levenberg-marquardt');

    p0 = fsolve(fsolve_traj_fun, p0, options);
end
y0 = cat(2,start_pos,start_vel,p0(1:3),p0(4:6))';

options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);   
options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 0.1, 10));
[t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
res_r_f_opt = norm(y(end, 1:3)-mars_r_f');
res_v_f_opt = norm(y(end, 4:6)-mars_v_f');
traj_opt = y(:,1:3);
disp(['------------------------------------'])
disp(['Невязка координаты аппроксимация ', num2str(res_r_f*r_unit/1000,'%10.2e\n'),',км'])
disp(['Невязка скорости аппроксимация ', num2str(res_v_f*V_unit/1000,'%10.2e\n'),',км/с'])
disp(['Невязка координаты оптимизация ', num2str(res_r_f_opt*r_unit/1000,'%10.2e\n'),',км'])
disp(['Невязка скорости оптимизация ', num2str(res_v_f_opt*V_unit/1000,'%10.2e\n'),',км/с'])

traj_New = quatrotate(q_eps',traj);
traj_opt_New = quatrotate(q_eps',traj_opt);
mars_r_f_New = quatrotate(q_eps',mars_r_f');
mars_v_f_New = quatrotate(q_eps',mars_v_f');

%рисуем
t0 = t_start;
t_orbit = linspace(t0,t0+T_earth/(24*3600), 1000);

st.t = t_orbit';
st.planet = 'Earth';
st.mode = orbits;

earth_traj = planetModel(st);
earth_traj=1e+3*earth_traj/ae;
earth_traj_New = quatrotate(q_eps',earth_traj);

t_orbit = linspace(t0,t0+T_mars/(24*3600), 1000);

st.t = t_orbit';
st.planet = planet_end;
st.mode = orbits;

mars_traj = planetModel(st);
mars_traj=1e+3*mars_traj/ae;
mars_traj_New = quatrotate(q_eps',mars_traj);

%выводим график
figure(4);
plot3(0, 0, 0, 'k--o', 'HandleVisibility','off')
set(gca,'FontSize',14)
hold on;

plot3(earth_traj_New(:, 1), earth_traj_New(:, 2), earth_traj_New(:, 3), 'k', 'HandleVisibility','off')
plot3(mars_traj_New(:, 1), mars_traj_New(:, 2), mars_traj_New(:, 3), 'k', 'HandleVisibility','off')
plot3(traj_New(:, 1), traj_New(:, 2), traj_New(:, 3), 'cyan', 'LineWidth', 1, 'DisplayName', 'Начальное приближение');
plot3(traj_opt_New(:, 1), traj_opt_New(:, 2), traj_opt_New(:, 3), 'magenta', 'LineStyle','--','LineWidth', 1, 'DisplayName', 'Оптимизированная траектория');


plot3(traj_New(end, 1), traj_New(end, 2), traj_New(end, 3),'bO', 'HandleVisibility','off')

plot3(mars_r_f_New(1), mars_r_f_New(2), mars_r_f_New(3),'rO', 'HandleVisibility','off')

hold off;
axis equal

title('Траектория КА')
xlabel('x, a.e.')
ylabel('y, a.e.')
zlabel('z, a.e.')
view(0,90)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
grid on;
box off;
legend('Location','best');
%%
function res = fsolve_traj(z,y0,yf,t0,t_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tspan = linspace(t0,t0+t_end, 10);
y0(7:12)=z;
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);   
options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 0.4, 5));
[~,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
res = y(end,1:6)-yf(1:6);
%res = [norm(y(end,1:3)-yf(1:3));norm(y(end,4:6)-yf(4:6))];
%res = [norm(y(end,1:2)-yf(1:2));y(end,3)-yf(3); norm(y(end,4:5)-yf(4:5));y(end,6)-yf(6);];
%res = sign(res).*log(sign(res).*res+1);
end

function [value, isterminal, direction] = eventIntegrationTrajStopR(~, y, r_end, r_outer)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


r=y(1:3);
value = [0,0];
value(1) = norm(r) - r_end;
if value(1) < 0
    value(1)=0;
end
value(2) = r_outer-norm(r);
if value(2) < 0
    value(2)=0;
end
isterminal = [1 1];
direction = [0 0];
end
function R = calculateRotMatrix(r0,v0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
e1 = r0/norm(r0);
e2 = v0/norm(v0);
e3 = cross(e1,e2);
R = [e1;e2;e3];
end
function res = internalIntegration3D(~,y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
mug=1;

res=zeros(12,1);
r=y(1:3);
v=y(4:6);
pr=y(7:9);
pv=y(10:12);

nr=norm(r);
res(1:3)=v;
res(4:6)=-mug*r/nr^3+pv;
res(7:9)=mug*pv/nr^3 - mug*(3*(r*r')/norm(r)^(5))*pv;
%res(7:9)=mug*ddUdrdr(r)*pv;
res(10:12)=-pr;
end