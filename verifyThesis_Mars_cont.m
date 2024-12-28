%внешние константы
ae = 149597870700;
mug_0 = 132712.43994*(10^6)*(10^(3*3));
T_earth_days = 365.256363004;
T_earth = T_earth_days*3600*24;
T_mars=T_earth*1.8808476;
T_mars_days = T_earth_days*1.8808476;
a_unit=(ae/sqrt(mug_0)).^2;
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
B=0.2721831;
A=0.5433279;
T_a  = @(a, a_rel)(a_rel-a_rel*B*log(a_rel))*a+(a_rel*B^2*log(a_rel));
a_OM = @(OM, a_rel,k)(OM+k+A*exp(-a_rel)*log(a_rel)*log(a_rel))/(A*(1+exp(-a_rel))*log(a_rel));
OM_a = @(a, a_rel)A*(a*(exp(a_rel)+1)-log(a_rel))*exp(-a_rel)*log(a_rel);
%% метод продолжения по грав параметру
t_start_orig = juliandate(2026,1,1);
N_count = 200;
N_count_2 = 50;
r_error_cont = zeros(N_count,N_count_2);
v_error_cont = zeros(N_count,N_count_2);

J_MARS_cont = zeros(N_count,N_count_2);
time_MARS_cont = zeros(N_count,N_count_2);
AN_MARS_2_cont = zeros(N_count,N_count_2);

T_MARS_test = linspace(400,3000, N_count);
t_start_MARS_test = linspace(t_start_orig,t_start_orig+2.135*T_earth_days, N_count_2);
savefilename = 'ThesisMarsCheck_cont.mat';
a_rel=1.52; %соотношение полуосей Земли и Марса
B=0.2721831;
case_traj = 2;
%%  метод продолжения по грав параметру
load('ThesisMarsCheck_cont.mat')
for i = 1:1% 
    disp([num2str(i), '/', num2str(N_count)])
    for j = 33:33
        % if time_MARS_cont(i,j) > 0
        %     continue
        % end
        disp(['    ', num2str(j), '/', num2str(N_count_2)])
        t_end = T_MARS_test(i)/T_earth_days; %дни переводим в годы
        AN_mars = (t_end-a_rel*B^2*log(a_rel))/(a_rel-a_rel*B*log(a_rel)); %получаем оптимальную угловую дальность
      
        n1 = floor(AN_mars);
        n2 = ceil(AN_mars);
        n3 = floor(AN_mars+1);
        n4 = ceil(AN_mars+1);
        n5 = floor(AN_mars-1);
        n6 = ceil(AN_mars-1);
        n_list = [n1,n2,n3,n4,n5,n6];
        n_list(n_list<0)=0;
        n_list=unique(n_list);

        evaluation_time_cont_best = 0;
        dr_cont_best = 0;
        dV_cont_best = 0;
        Jt_cont_best = 1e6;
        pr0_best = [0, 0, 0];
        pv0_best = [0, 0, 0];
        AN_best = 0;
        for n = n_list
            %pause(1);
            z0=zeros([1, 6]);
            tspan = linspace(0,t_end*2*pi, round(AN_mars*100));
            [rr_cont, Jt_cont, C_cont, evaluation_time_cont, dr_cont, dV_cont, pr0, pv0] =...
                checkContinuation(t_start_MARS_test(j), t_end*T_earth_days, tspan*T_unit,z0, case_traj,planet_end,1, n, orbits, 0);
            functional_cont = Jt_cont(end)/T_unit;
            if functional_cont<Jt_cont_best
                Jt_cont_best = functional_cont;
                pr0_best = pr0(1,:);
                pv0_best = pv0(1,:);
                evaluation_time_cont_best = evaluation_time_cont;
                dr_cont_best = dr_cont/ae;
                dV_cont_best = dV_cont/V_unit;
                AN_best = calculate_angular_distance(rr_cont(:,1:3))/(2*pi);
            end
        end
        if AN_mars<0.5
            [rr_cont, Jt_cont, C_cont, evaluation_time_cont, dr_cont, dV_cont, pr0, pv0] =...
                    ContinuationSimple(t_start_MARS_test(j), t_end*T_earth_days, tspan*T_unit,z0, case_traj,planet_end,1, n, orbits, 0);
            functional_cont = Jt_cont(end)/T_unit;
            if functional_cont<Jt_cont_best
                Jt_cont_best = functional_cont;
                pr0_best = pr0(1,:);
                pv0_best = pv0(1,:);
                evaluation_time_cont_best = evaluation_time_cont;
                dr_cont_best = dr_cont/ae;
                dV_cont_best = dV_cont/V_unit;
                AN_best = calculate_angular_distance(rr_cont(:,1:3))/(2*pi);
            end
        end
        r_error_cont(i,j) = dr_cont_best;
        v_error_cont(i,j) = dV_cont_best;
        J_MARS_cont(i,j) = Jt_cont_best;
        time_MARS_cont(i,j) = evaluation_time_cont_best;
        AN_MARS_2_cont(i,j) = AN_best;
        %save(savefilename);
    end
end
%% сравниваем дельты
figure(1)
%r_error_opt_fix = a_unit^2*J_MARS_cont;
%r_error_opt_fix = a_unit^2*J_MARS_cont-J_MARS_opt;
AN_diff_1 = AN_MARS_2_cont-AN_MARS_test;
AN_diff_2 = -AN_diff_1-delta_MARS_opt;
criteria_2 = abs(AN_diff_2)<0.5;
AN_diff_2(~criteria_2)=nan;
s = surf(T_MARS_test/T_earth_days,t_start_MARS_test-t_start_MARS_test(1), AN_diff_2()','HandleVisibility','off');
s.EdgeColor = 'none';
hold on;
hold off;
grid;
legend('Location','best');
xlabel('Время, годы.')
ylabel('Дата старта')
zlabel('Невязка на правом конце')
%set(gca, 'ZScale', 'log')
%title('Невязка на правом конце')
yticks([0 (t_start_MARS_test(end)-t_start_MARS_test(1))/2 t_start_MARS_test(end)-t_start_MARS_test(1)])
yticklabels({'01.01.2026','26.01.2027','20.02.2028'})
colorbar;
grid;
%set(gca,'ColorScale','log')
view(0,90)
%%
a_unit=(ae/sqrt(mug_0)).^2;
J_GOC_ephemeris = zeros(N_count,1);
t_start_GOC_ephemeris = zeros(N_count,1);
J_GOC_ephemeris_pred = zeros(N_count,1);
time_GOC_ephemeris = zeros(N_count,1);

for i = 1:N_count
    J_local = J_MARS_cont(i,:);
    J_min = min(J_local,[],'all');
    j_min = find(J_local==J_min);

    J_GOC_ephemeris(i) = a_unit^2*J_min;
    AN_i = AN_MARS_test(i,1);
    J_GOC_ephemeris_pred(i) = J_a(AN_i);
    time_GOC_ephemeris(i) = time_MARS_cont(i,j_min);
    t_start_GOC_ephemeris(i) = t_start_MARS_test(j_min)-t_start_MARS_test(1);
end

figure(4)
%plot(T_MARS_test/T_earth_days, (J_GOC_ephemeris_pred-J_GOC_ephemeris)./J_GOC_ephemeris, 'DisplayName', 'Целевой функционал')
plot(T_MARS_test/T_earth_days, J_GOC_ephemeris, 'DisplayName', 'Целевой функционал')
hold on;
plot(T_MARS_test/T_earth_days, J_GOC_ephemeris_pred, 'DisplayName', 'Аппроксимация')
hold off;
grid;
xlabel('Время, годы.')
ylabel('Функционал')
title('Сравнение функционалов')
legend('Location','best');
%%
AN_optimal_pred = zeros(N_count_2,1);
a_rel = 1.52;
for j=1:N_count_2
    t_start_j = t_start_MARS_test(j);
    st.mode = orbits;
    st.t = t_start_j;
    st.planet = planet_end;
    st.mode = orbits;
    
    [mars_r_f, mars_v_f]=planetModel(st);
    mars_r_f=mars_r_f'*1e+03/r_unit;
    mars_r_f_new = quatrotate(q_eps',mars_r_f');
    mars_an_f = atan2(mars_r_f_new(2), mars_r_f_new(1));

    st.planet = planet_start;
    [earth_r_f, earth_v_f]=planetModel(st);
    earth_r_f=earth_r_f'*1e+03/r_unit;
    earth_r_f_new = quatrotate(q_eps',earth_r_f');
    earth_an_f = atan2(earth_r_f_new(2), earth_r_f_new(1));
    an_diff = -(earth_an_f-mars_an_f)/(2*pi);%приводим к виткам
    AN_optimal_pred(j) = an_diff;
end
%%
T_optimal_pred = zeros(N_count_2,20);
for j=1:N_count_2
    i=1;
    an_diff = AN_optimal_pred(j);
    for k =-10:10
        a_OM_k = a_OM(an_diff, a_rel,k);
        if a_OM_k>0.75
            T_optimal_pred(j,i) = T_a(a_OM_k, a_rel);
            i=i+1;
        end
    end
end
%% метод продолжения по грав параметру
figure(1)
%r_error_opt_fix = a_unit^2*J_MARS_cont;
%r_error_opt_fix = a_unit^2*J_MARS_cont-J_MARS_opt;
r_error_opt_fix = r_error_opt;
r_error_opt_fix(r_error_opt_fix>1) = 1;
AN_diff = AN_MARS_2_cont-AN_MARS_2_test;
criteria = abs(AN_diff)>1e-3; %1-sum(sum(criteria))/(N_count*N_count_2)
%1-sum(sum(criteria&criteria_2))/(N_count*N_count_2)
AN_diff(AN_diff>3)=3;
AN_diff(AN_diff<-3)=-3;
%r_error_opt_fix(criteria)=nan;
s = surf(T_MARS_test/T_earth_days,t_start_MARS_test-t_start_MARS_test(1), AN_diff','HandleVisibility','off');
s.EdgeColor = 'none';
hold on;
plot3(T_MARS_test(1:36)/T_earth_days,t_start_GOC_ephemeris(1:36),J_GOC_ephemeris(1:36),'r', 'DisplayName', 'Оптимум по дате старта')
plot3(T_MARS_test(37:191)/T_earth_days,t_start_GOC_ephemeris(37:191),J_GOC_ephemeris(37:191),'r','HandleVisibility','off')
plot3(T_MARS_test(192:end)/T_earth_days,t_start_GOC_ephemeris(192:end),J_GOC_ephemeris(192:end),'r','HandleVisibility','off')

plot3(T_optimal_pred(1:19,1),t_start_MARS_test(1:19)-t_start_MARS_test(1),T_optimal_pred(1:19,1)*0+3,'k', 'DisplayName', 'Предсказанный оптимум')
plot3([T_optimal_pred(1:19,2);T_optimal_pred(20:end,1)],...
    [t_start_MARS_test(1:19),t_start_MARS_test(20:end)]-t_start_MARS_test(1),...
    [T_optimal_pred(1:19,1);T_optimal_pred(20:end,1)]*0+3,'k','HandleVisibility','off')

plot3(T_optimal_pred(20:end,2),t_start_MARS_test(20:end)-t_start_MARS_test(1),T_optimal_pred(20:end,2)*0+3,'k','HandleVisibility','off')
xlim([T_MARS_test(1), T_MARS_test(end)]/T_earth_days)
%plot3(,-t_start_MARS_test(1),*0+3,'k')
hold off;
grid;
legend('Location','best');
xlabel('Время, годы.')
ylabel('Дата старта')
zlabel('Невязка на правом конце')
%set(gca, 'ZScale', 'log')
%title('Невязка на правом конце')
yticks([0 (t_start_MARS_test(end)-t_start_MARS_test(1))/2 t_start_MARS_test(end)-t_start_MARS_test(1)])
yticklabels({'01.01.2026','26.01.2027','20.02.2028'})
colorbar;
grid;
%set(gca,'ColorScale','log')
view(0,90)

%%
figure(5)
%plot(T_MARS_test/T_earth_days, (J_GOC_ephemeris_pred-J_GOC_ephemeris)./J_GOC_ephemeris, 'DisplayName', 'Целевой функционал')
plot(T_MARS_test/T_earth_days, time_GOC_ephemeris, 'DisplayName', 'Целевой функционал')
hold on;
%plot(T_MARS_test/T_earth_days, J_GOC_ephemeris_pred, 'DisplayName', 'Аппроксимация')
hold off;
grid;
xlabel('Время, годы.')
ylabel('Время оптимизации, сек.')
title('Время оптимизации')
%legend('Location','best');
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
e1 = r0/norm(r0);
e2 = v0/norm(v0);
e3 = cross(e1,e2);
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