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

%VENUS
pr_x_V = @(a,d)(-0.2489269721475946)/(7.737402704800786*a-0.32876870077839637...
    +sin(2*pi*a+0.10377485561964216+d*0.5010768929874032/(a^2-1.9514741381243845*a+2.4016982504424953))...
    + d*(-39.740029616989275+73.28781124258543*a-40.57760476908843*a^2)*exp(-3.086139681061194*a) ); %OK
pr_y_V = @(a,d)(-0.006779491667620585 + 0.00432410667098884 * cos(2*pi * a...
    +d*(-5.520942046733101*a^2+9.93390376988929*a-4.541315615665989)*exp(1.3808955605425999*a-2.0000002462936357*a^2)))/(a^2+0.08726450366262754*a-0.20030536675890054);%OK
pv_x_V = @(a,d)(-0.006319471111610397 + 0.004278503661785389 * cos(2*pi * a...
    +d*(-1.786880982902473*a^2+3.1517137043764167*a-1.3891212954291232)*exp(2.9860241498160076*a-2.699241528126778*a^2)))/(a^2+0.13118547856852233*a-0.24254361523409607);
pv_y_V = @(a,d)(-0.12867786064634879)/(4.0002764433454105*a-0.15629991248717665...
    +sin(2*pi*a+0.08483203399561368+d*0.008896394730198013/(a^2-2.9131929806240318*a+2.1648040548553382)) +...
    d*(-1354.631280683438*a^2+2264.6420824634483*a-951.3449842396692)*exp(-6.814266658923464*a));

J_a  = @(a)(0.002117*exp(a))/(a*exp(a)+0.17856*(sin(7.44443*a)-a));
%% перебираем разное время
N_count = 300;
r_error_approx = zeros(N_count,1);
v_error_approx = zeros(N_count,1);

r_error_opt = zeros(N_count,1);
v_error_opt = zeros(N_count,1);
J_MARS_opt = zeros(N_count,1);
delta_MARS_opt = zeros(N_count,1);
delta_MARS_opt = zeros(N_count,1);
AN_MARS_test = zeros(N_count,1);
J_MARS_approx= zeros(N_count,1);
T_MARS_test = linspace(400,3000, N_count);

a_rel=1.52; %соотношение полуосей Земли и Марса
B=0.2721831;
for i = 1:N_count 
    disp([num2str(i), '/', num2str(N_count)])
    t_end = T_MARS_test(i)/T_earth_days; %дни переводим в годы
    AN_mars = (t_end-a_rel*B^2*log(a_rel))/(a_rel-a_rel*B*log(a_rel)); %получаем оптимальную угловую дальность
    AN_MARS_test(i) = AN_mars;
    J_MARS_approx(i) = J_a(AN_mars);
    pr_0 = R^(-1)*[pr_x(AN_mars);pr_y(AN_mars);0];
    pv_0 = R^(-1)*[pv_x(AN_mars);pv_y(AN_mars);0]; 
    tspan = linspace(0,t_end*2*pi, round(AN_mars*100));
    
    y0 = cat(2,start_pos,start_vel,pr_0',pv_0')';
    
    options = odeset('AbsTol',1e-10);
    options = odeset(options,'RelTol',1e-10);   
    options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 0.1));
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    st.mode = orbits;
    st.t = t_start+t_end*T_earth_days;
    %st.t = [t_start, t(end)*T_unit/(3600*24)];
    st.planet = planet_end;
    st.mode = orbits;
    
    [mars_r_f, mars_v_f]=planetModel(st);
    mars_r_f=mars_r_f'*1e+03/r_unit;
    mars_v_f=mars_v_f'*1e+03/V_unit;
    
    res_r_f = norm(y(end, 1:3)-mars_r_f');
    res_v_f = norm(y(end, 4:6)-mars_v_f');
    yf = [mars_r_f;mars_v_f]';
    r_error_approx(i) = res_r_f;
    v_error_approx(i) = res_v_f;
    
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
    delta_MARS_opt(i) = delta_theta_f;

    p0 = [pr_0;pv_0]';
    fsolve_traj_fun=@(z)fsolve_traj(z,y0,yf,0,t_end*2*pi);
    options = optimoptions('fsolve','Display','off');
    %options = optimoptions(options, 'Algorithm', 'levenberg-marquardt');

    p0 = fsolve(fsolve_traj_fun, p0, options);
    y0 = cat(2,start_pos,start_vel,p0(1:3),p0(4:6))';

    options = odeset('AbsTol',1e-10);
    options = odeset(options,'RelTol',1e-10);   
    options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 0.1));
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    a_reactive = y(:,10:12);
    a_vec=vecnorm(a_reactive, 2, 2).^2;
    Jt = cumtrapz(t, a_vec)/(2);
    res_r_f_opt = norm(y(end, 1:3)-mars_r_f');
    res_v_f_opt = norm(y(end, 4:6)-mars_v_f');
    r_error_opt(i) = res_r_f_opt;
    v_error_opt(i) = res_v_f_opt;
    J_MARS_opt(i) = Jt(end);
    if res_r_f_opt > 1e-16
        AN_mars = AN_mars-delta_theta_f;
        pr_0 = R^(-1)*[pr_x(AN_mars);pr_y(AN_mars);0];
        pv_0 = R^(-1)*[pv_x(AN_mars);pv_y(AN_mars);0]; 
        p0 = [pr_0;pv_0]';
        fsolve_traj_fun=@(z)fsolve_traj(z,y0,yf,0,t_end*2*pi);
        options = optimoptions('fsolve','Display','off');
        p0 = fsolve(fsolve_traj_fun, p0, options);
        y0 = cat(2,start_pos,start_vel,p0(1:3),p0(4:6))';
    
        options = odeset('AbsTol',1e-10);
        options = odeset(options,'RelTol',1e-10);   
        options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 0.1));
        [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
        a_reactive = y(:,10:12);
        a_vec=vecnorm(a_reactive, 2, 2).^2;
        Jt = cumtrapz(t, a_vec)/(2);
        res_r_f_opt = norm(y(end, 1:3)-mars_r_f');
        res_v_f_opt = norm(y(end, 4:6)-mars_v_f');
        if Jt(end) < J_MARS_opt(i)
            r_error_opt(i) = res_r_f_opt;
            v_error_opt(i) = res_v_f_opt;
            J_MARS_opt(i) = Jt(end);
        end
    end


    r_error_approx(i) = res_r_f;
    v_error_approx(i) = res_v_f;
end
figure(1)
plot(T_MARS_test/T_earth_days, r_error_approx,'--', 'DisplayName', 'Невязка положения аппроксимация')
hold on;
plot(T_MARS_test/T_earth_days, v_error_approx,'--', 'DisplayName', 'Невязка скорости аппроксимация')
plot(T_MARS_test/T_earth_days, r_error_opt, 'DisplayName', 'Невязка положения оптимизация')
plot(T_MARS_test/T_earth_days, v_error_opt, 'DisplayName', 'Невязка скорости оптимизация')
hold off
legend('Location','best');
grid;
xlabel('Время, годы.')
ylabel('Невязка на правом конце')
set(gca, 'YScale', 'log')
figure(2)
J_MARS_opt_fix = J_MARS_opt;
J_MARS_opt_fix(r_error_opt>1e-5) = nan;
J_MARS_opt_fix(v_error_opt>1e-5) = nan;
plot(T_MARS_test/T_earth_days, J_MARS_opt_fix, 'DisplayName', 'Целевой функционал')
hold on;
plot(T_MARS_test/T_earth_days, J_MARS_approx, 'DisplayName', 'Аппроксимация')
hold off;
grid;
xlabel('Время, годы.')
ylabel('Функционал')
figure(3)
plot(T_MARS_test/T_earth_days, delta_MARS_opt, 'DisplayName', 'Оригинал')
hold on;
%plot(T_MARS_test/T_earth_days, delta_MARS_opt_2, 'DisplayName', 'Сдвиг')
hold off;
grid;
legend;
xlabel('Время, годы.')
ylabel('Угол, витки')
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
savefilename = 'ThesisMarsCheck_naive_apply.mat';
a_rel=1.52; %соотношение полуосей Земли и Марса
B=0.2721831;
%%
for i = 1:N_count% 
    disp([num2str(i), '/', num2str(N_count)])
    for j = 1:N_count_2 
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
        tspan = linspace(0,t_end*2*pi, round(AN_mars*100));
        
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
    
        
        
        if abs(delta_theta_f) >= 0
            AN_mars = AN_mars-0.3*delta_theta_f;
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
            %options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 0.4, 5));
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
        else
            p0 = [pr_0;pv_0]';
            fsolve_traj_fun=@(z)fsolve_traj(z,y0,yf,0,t_end*2*pi);
            options = optimoptions('fsolve','Display','off');
            % options = optimoptions(options, 'Algorithm', 'levenberg-marquardt');
            % options = bvpset('AbsTol',1e-10);
            % options = bvpset(options,'RelTol',1e-10);  
            % options = bvpset(options,'Nmax',round(AN_mars*300));  
            % id = 'MATLAB:bvp5c:RelTolNotMet';
            % warning('off',id)
            % solinit = bvpinit(tspan,@(x)guess(x,y0,p0));
            tic 
            % sol  = bvp5c(@internalIntegration3D,@(ya,yb)bcfcn(ya,yb, y0,yf'),solinit,options);
            % p0 = sol.y(7:12,1)';
            p0 = fsolve(fsolve_traj_fun, p0, options);
            time_MARS_opt(i,j) = toc;
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
            r_error_opt(i,j) = res_r_f_opt;
            v_error_opt(i,j) = res_v_f_opt;
            J_MARS_opt(i,j) = Jt(end);
            P_MARS_test(i,j,:) = p0;
            AN_MARS_2_test(i,j) = calculate_angular_distance(y(:,1:3))/(2*pi);
        end
        save(savefilename);
    end
end
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
for i = 199:-1:100% 
    disp([num2str(i), '/', num2str(N_count)])
    for j = 1:N_count_2
        if time_MARS_cont(i,j) > 0
            continue
        end
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
            pause(1);
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
        
        r_error_cont(i,j) = dr_cont_best;
        v_error_cont(i,j) = dV_cont_best;
        J_MARS_cont(i,j) = Jt_cont_best;
        time_MARS_cont(i,j) = evaluation_time_cont_best;
        AN_MARS_2_cont(i,j) = AN_best;
        save(savefilename);
    end
end
%% метод продолжения по грав параметру
figure(1)
r_error_opt_fix = AN_MARS_2_cont;
s = surf(T_MARS_test/T_earth_days,t_start_MARS_test-t_start_MARS_test(1), r_error_opt_fix');
s.EdgeColor = 'none';
grid;
xlabel('Время, годы.')
ylabel('Дата старта')
zlabel('Невязка на правом конце')
%set(gca, 'ZScale', 'log')
title('Невязка на правом конце')
yticks([0 (t_start_MARS_test(end)-t_start_MARS_test(1))/2 t_start_MARS_test(end)-t_start_MARS_test(1)])
yticklabels({'01.01.2026','26.01.2027','20.02.2028'})
colorbar;
grid;
%set(gca,'ColorScale','log')
view(0,90)

%%
figure(2)
AN_delta = AN_MARS_2_test-AN_MARS_test;
criteria = abs(AN_delta)>1 | r_error_opt>0.1;
%criteria = abs(delta_MARS_opt)>0.1  | abs(AN_delta)>0.9;
AN_delta(criteria) = nan;
s = surf(T_MARS_test/T_earth_days,t_start_MARS_test-t_start_MARS_test(1), AN_delta');
s.EdgeColor = 'none';
grid;
xlabel('Время, годы.')
ylabel('Дата старта')
zlabel('Витки')
title('Разница угловой дальности')
hold on;
%s_contour = contour3(T_MARS_test/T_earth_days,t_start_MARS_test-t_start_MARS_test(1), AN_delta',[-0.1,0,0.1], 'ShowText','on', 'HandleVisibility','off', 'Color', 'k');

hold off;

yticks([0 (t_start_MARS_test(end)-t_start_MARS_test(1))/2 t_start_MARS_test(end)-t_start_MARS_test(1)])
yticklabels({'01.01.2026','26.01.2027','20.02.2028'})
colorbar;
%set(gca,'ColorScale','log')
%set(gca, 'ZScale', 'log')
view(0,90)
%%

figure(1)
r_error_opt_fix = r_error_opt;
%r_error_opt_fix(r_error_opt_fix>1e2) = 1e2;
r_error_opt_fix(criteria) = nan;
s = surf(T_MARS_test/T_earth_days,t_start_MARS_test-t_start_MARS_test(1), r_error_opt_fix');
s.EdgeColor = 'none';
grid;
xlabel('Время, годы.')
ylabel('Дата старта')
zlabel('Невязка на правом конце')
%set(gca, 'ZScale', 'log')
title('Невязка на правом конце')
yticks([0 (t_start_MARS_test(end)-t_start_MARS_test(1))/2 t_start_MARS_test(end)-t_start_MARS_test(1)])
yticklabels({'01.01.2026','26.01.2027','20.02.2028'})
colorbar;
grid;
%set(gca,'ColorScale','log')
view(0,90)

%%
figure(3)

delta_MARS_opt_contour = [-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4];
%delta_MARS_opt_contour = [-0.1, 0, 0.1];
J_MARS_opt_fix = AN_diff;
J_MARS_opt_fix(J_MARS_opt_fix>1e2) = 1e2;
%J_MARS_opt_fix(criteria) = nan;
s_contour = contour3(T_MARS_test/T_earth_days,t_start_MARS_test-t_start_MARS_test(1), delta_MARS_opt',delta_MARS_opt_contour, 'ShowText','on', 'HandleVisibility','off', 'Color', 'k');
s = surf(T_MARS_test/T_earth_days,t_start_MARS_test-t_start_MARS_test(1), J_MARS_opt_fix');
s.EdgeColor = 'none';
grid;
xlabel('Время, годы.')
ylabel('Дата старта')
zlabel('Функционал')
title('Функционал')
hold on;
K_contours = 1;
C_contours = 1;
C_contours_skip = [0, 0];
while K_contours<length(s_contour)
    H_contours = s_contour(1,K_contours);
    N_contours = s_contour(2,K_contours);
    X_contours = s_contour(1,K_contours+1:K_contours+N_contours);
    Y_contours = s_contour(2,K_contours+1:K_contours+N_contours);
    Z_contours = interp2(T_MARS_test/T_earth_days,t_start_MARS_test-t_start_MARS_test(1),J_MARS_opt_fix', X_contours, Y_contours);
    K_contours = K_contours+N_contours+1;
    C_contours = C_contours+1;
    %plot3(X_contours, Y_contours,Z_contours,'black', 'HandleVisibility','off')
    %if any(C_contours_skip(:) == C_contours)
     %   continue
    %end
    pos_i = 1;
    pos_shift = 0;
    plot3(X_contours(1:end), Y_contours(1:end),Z_contours(1:end),'black','HandleVisibility','off')
    %plot3(X_contours(1:pos_i), Y_contours(1:pos_i),Z_contours(1:pos_i),'black','HandleVisibility','off')
    %plot3(X_contours(pos_i+pos_shift:end), Y_contours(pos_i+pos_shift:end),Z_contours(pos_i+pos_shift:end),'black','HandleVisibility','off')
    text_label = num2str(H_contours);
    h = text((X_contours(pos_i)+X_contours(pos_i+pos_shift))/2,...
        -6+(Y_contours(pos_i)+Y_contours(pos_i+pos_shift))/2,...
        10+(Z_contours(pos_i)+Z_contours(pos_i+pos_shift))/2,...
        text_label,'HorizontalAlignment','center',...
     'VerticalAlignment','Bottom','FontName','consolas','FontSize',11); 
end
hold off;
yticks([0 (t_start_MARS_test(end)-t_start_MARS_test(1))/2 t_start_MARS_test(end)-t_start_MARS_test(1)])
yticklabels({'01.01.2026','26.01.2027','20.02.2028'})
colorbar;
%set(gca,'ColorScale','log')
%set(gca, 'ZScale', 'log')
view(0,90)
%%

J_GOC_ephemeris = zeros(N_count,1);
J_GOC_ephemeris_pred = zeros(N_count,1);
time_GOC_ephemeris = zeros(N_count,1);

for i = 1:N_count
    J_local = J_MARS_opt_fix(i,:);
    J_min = min(J_local,[],'all');
    j_min = find(J_local==J_min);

    J_GOC_ephemeris(i) = J_min;
    AN_i = AN_MARS_test(i,j_min);
    J_GOC_ephemeris_pred(i) = J_a(AN_i);
    time_GOC_ephemeris(i) = time_MARS_opt(i,j_min);
    
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
%% получаем решение для графика и строим график
t_end = T_MARS_test(10)/T_earth_days; %дни переводим в годы
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
%options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 0.1));
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
%options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 0.1));
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
% [a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(mars_r_f_New',mars_v_f_New',1);
% disp(['e=', num2str(eMag)])
% disp(['a=', num2str(a)])
% disp(['i=', num2str(i)])
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
