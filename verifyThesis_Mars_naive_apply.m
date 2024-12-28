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
pv_y = @(a)(0.10456877009148134)/(4.050922375438426*a-0.15902042687671913+sin(2*pi*a+0.11429217276196703+0.012218372217935064/(a^2-2.8347986183365705*a+2.06664561825753)) +(-540.5981059117084*a^2+958.0032106019531*a-415.15730979365355)*exp(-5.877538308252977*a));

J_a  = @(a)(0.002117*exp(a))/(a*exp(a)+0.17856*(sin(7.44443*a)-a));

%% перебираем разное время и дату старта
t_start_orig = juliandate(2026,1,1);
N_count = 200;
N_count_2 = 50;
r_error_approx = zeros(N_count,N_count_2);
v_error_approx = zeros(N_count,N_count_2);

r_error_opt = zeros(N_count,N_count_2);
v_error_opt = zeros(N_count,N_count_2);
J_MARS_opt = zeros(N_count,N_count_2);
time_MARS_opt = zeros(N_count,N_count_2);
delta_MARS_opt = zeros(N_count,N_count_2);

P_MARS_test = zeros(N_count,N_count_2,6);
AN_MARS_test = zeros(N_count,N_count_2);
AN_MARS_2_test = zeros(N_count,N_count_2);
J_MARS_approx= zeros(N_count,N_count_2);
T_MARS_test = linspace(400,3000, N_count);
t_start_MARS_test = linspace(t_start_orig,t_start_orig+2.135*T_earth_days, N_count_2);
savefilename = 'ThesisMarsCheck_naive_apply_6.mat';
a_rel=1.52; %соотношение полуосей Земли и Марса
B=0.2721831;
%%
load('ThesisMarsCheck_naive_apply_6.mat')
for i = 1:N_count% 
    disp([num2str(i), '/', num2str(N_count)])
    for j = 1:N_count_2 
        if AN_MARS_2_test(i,j)>0
            continue
        end
        %определяем стартовое положение
        st.t = t_start_MARS_test(j);
        st.planet = planet_start;
        st.mode = orbits;
        
        [start_pos, start_vel] = planetModel(st);
        start_pos=start_pos*1e+03/ae;
        start_vel=start_vel*1e+03/V_unit;
        R = calculateRotMatrix(start_pos,start_vel);
        t_end = T_MARS_test(i)/T_earth_days; %дни переводим в годы
        AN_mars = (t_end-a_rel*B^2*log(a_rel))/(a_rel-a_rel*B*log(a_rel)); %получаем оптимальную угловую дальность
        AN_MARS_test(i,j) = AN_mars;
        J_MARS_approx(i,j) = J_a(AN_mars);
        pr_0 = R^(-1)*[pr_x(AN_mars);pr_y(AN_mars);0];
        pv_0 = R^(-1)*[pv_x(AN_mars);pv_y(AN_mars);0]; 
        tspan = linspace(0,t_end*2*pi, round(AN_mars*1000));
        
        y0 = cat(2,start_pos,start_vel,pr_0',pv_0')';
        
        options = odeset('AbsTol',1e-10);
        options = odeset(options,'RelTol',1e-10);   
        options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 0.4,5));
        [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
        
        st.mode = orbits;
        st.t = t_start_MARS_test(j)+t_end*T_earth_days;
        %st.t = [t_start, t(end)*T_unit/(3600*24)];
        st.planet = planet_end;
        st.mode = orbits;
        
        [mars_r_f, mars_v_f]=planetModel(st);
        mars_r_f=mars_r_f'*1e+03/r_unit;
        mars_v_f=mars_v_f'*1e+03/V_unit;
        
        res_r_f = norm(y(end, 1:3)-mars_r_f');
        res_v_f = norm(y(end, 4:6)-mars_v_f');
        yf = [mars_r_f;mars_v_f]';
        r_error_approx(i,j) = res_r_f;
        v_error_approx(i,j) = res_v_f;
        
        mars_r_f_new = quatrotate(q_eps',mars_r_f');
        mars_an_f = atan2(mars_r_f_new(2), mars_r_f_new(1));
        traj_new = quatrotate(q_eps',y(:,1:3));
        traj_an_f = atan2(traj_new(end,2), traj_new(end,1));
        delta_theta_f = (traj_an_f-mars_an_f)/(2*pi);
        if delta_theta_f > 0.5
            delta_theta_f=delta_theta_f-1;
        elseif delta_theta_f < -0.5
            delta_theta_f=delta_theta_f+1;
        end
        delta_MARS_opt(i,j) = delta_theta_f;
    
        AN_mars_delta = -delta_theta_f/0.29;%более точное значение /0.29
        AN_mars_delta = clip(AN_mars_delta, -1,1);
        %AN_mars_delta = 0;
        AN_mars = AN_mars+AN_mars_delta;
        if AN_mars<=0.75
            AN_mars=0.75;
        end
        pr_0 = R^(-1)*[pr_x(AN_mars);pr_y(AN_mars);0];
        pv_0 = R^(-1)*[pv_x(AN_mars);pv_y(AN_mars);0]; 
        p0 = [pr_0;pv_0]';
        fsolve_traj_fun=@(z)fsolve_traj(z,y0,yf,0,t_end*2*pi);
        options = optimoptions('fsolve','Display','off');
        p0 = fsolve(fsolve_traj_fun, p0, options);
        y0 = cat(2,start_pos,start_vel,p0(1:3),p0(4:6))';
    
        options = odeset('AbsTol',1e-10);
        options = odeset(options,'RelTol',1e-10);   
        options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 0.4, 5));
        [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
        a_reactive = y(:,10:12);
        a_vec=vecnorm(a_reactive, 2, 2).^2;
        Jt = cumtrapz(t, a_vec)/(2);
        res_r_f_opt = norm(y(end, 1:3)-mars_r_f');
        res_v_f_opt = norm(y(end, 4:6)-mars_v_f');
        %if Jt(end) < J_MARS_opt(i)
        r_error_opt(i,j) = res_r_f_opt;
        v_error_opt(i,j) = res_v_f_opt;
        J_MARS_opt(i,j) = Jt(end);
        P_MARS_test(i,j,:) = p0;
        AN_MARS_2_test(i,j) = calculate_angular_distance(y(:,1:3))/(2*pi);
        %end
        % if delta_theta_f>=0.2 %проверяем второе семейство
        %     AN_mars = AN_mars+delta_theta_f/0.3;%более точное значение /0.28
        %     if AN_mars<=0.75
        %         AN_mars=0.75;
        %     end
        %     pr_0 = R^(-1)*[pr_x(AN_mars);pr_y(AN_mars);0];
        %     pv_0 = R^(-1)*[pv_x(AN_mars);pv_y(AN_mars);0]; 
        %     p0 = [pr_0;pv_0]';
        %     fsolve_traj_fun=@(z)fsolve_traj(z,y0,yf,0,t_end*2*pi);
        %     options = optimoptions('fsolve','Display','off');
        %     p0 = fsolve(fsolve_traj_fun, p0, options);
        %     y0 = cat(2,start_pos,start_vel,p0(1:3),p0(4:6))';
        % 
        %     options = odeset('AbsTol',1e-10);
        %     options = odeset(options,'RelTol',1e-10);   
        %     options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 0.4, 5));
        %     [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
        %     a_reactive = y(:,10:12);
        %     a_vec=vecnorm(a_reactive, 2, 2).^2;
        %     Jt = cumtrapz(t, a_vec)/(2);
        %     res_r_f_opt = norm(y(end, 1:3)-mars_r_f');
        %     res_v_f_opt = norm(y(end, 4:6)-mars_v_f');
        %     if Jt(end)<J_MARS_opt(i,j) && res_r_f_opt<1e-2
        %         %if Jt(end) < J_MARS_opt(i)
        %         r_error_opt(i,j) = res_r_f_opt;
        %         v_error_opt(i,j) = res_v_f_opt;
        %         J_MARS_opt(i,j) = Jt(end);
        %         P_MARS_test(i,j,:) = p0;
        %         AN_MARS_2_test(i,j) = calculate_angular_distance(y(:,1:3))/(2*pi);
        %     end
        % end
        save(savefilename);
    end
end
%%
figure(1)
r_error_opt_fix = r_error_opt;
r_error_opt_fix(r_error_opt_fix>1e2) = 1e2;
r_error_opt_fix(r_error_opt==0) = nan;
s = surf(T_MARS_test/T_earth_days,t_start_MARS_test-t_start_MARS_test(1), r_error_opt_fix');
s.EdgeColor = 'none';
grid;
xlabel('Время, годы.')
ylabel('Дата старта')
zlabel('Невязка на правом конце')
set(gca, 'ZScale', 'log')
title('Невязка на правом конце')
yticks([0 (t_start_MARS_test(end)-t_start_MARS_test(1))/2 t_start_MARS_test(end)-t_start_MARS_test(1)])
yticklabels({'01.01.2026','26.01.2027','20.02.2028'})
colorbar;
grid;
set(gca,'ColorScale','log')
view(0,90)
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

function [c, ceq] = fmincon_traj(z,y0,yf,t0,t_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tspan = linspace(t0,t0+t_end, 10);
y0(7:12)=z;
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);   
[~,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
c =[];
ceq = (y(end,1:6)-yf(1:6));
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
e1 = r0;
e2 = v0;
%e2_tilde = v0 - dot(v0, e1) * e1;
%e2 = e2_tilde / norm(e2_tilde);
e3 = cross(e1,e2);
%e3 = e3/norm(e3);
%e2 = cross(e3,e1);
R = [e1;e2;e3];
end

function y_g = guess(x,y0,z) % guess at solution behavior
tspan = linspace(0,x, 2);
y0(7:12)=z;
if x>0
    options = odeset('AbsTol',1e-10);
    options = odeset(options,'RelTol',1e-10);   
    options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 0.4, 5));
    [~,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    y_g = y(end,:);
else
    y_g = y0;
end
end

function res = bcfcn(ya, yb, y0, yf) % boundary conditions
res = [ya(1:6) - y0(1:6);
       yb(1:6) - yf(1:6)];
end