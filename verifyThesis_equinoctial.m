%внешние константы
ae = 149597870700;
mug_0 = 132712.43994*(10^6)*(10^(3*3));
mug=1;
T_earth = 365.256363004*3600*24;
T_earth_days = 365.256363004;
T_mars=T_earth*1.8808476;
T_mars_days = 365.256363004*1.8808476;


r_unit=ae;
V_unit=sqrt(mug_0/ae);
T_unit = T_earth/(2*pi);
%поворот кватернионом в плоскость эклиптики
n_eps = [1; 0; 0];
eps_0 = 2*pi*(23+26/60+21/3600)/360;
q_eps = [cos(eps_0/2); n_eps*sin(eps_0/2)];
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

B=0.2721831;
A=0.5433279;
T_a  = @(a, a_rel)(a_rel-a_rel*B*log(a_rel))*a+(a_rel*B^2*log(a_rel));
a_OM = @(OM, a_rel,k)(OM+k+A*exp(-a_rel)*log(a_rel)*log(a_rel))/(A*(1+exp(-a_rel))*log(a_rel));
OM_a = @(a, a_rel)A*(a*(exp(a_rel)+1)-log(a_rel))*exp(-a_rel)*log(a_rel);
%% перебираем разные дальности
N_count = 991;
r_error_approx = zeros(N_count,1);
v_error_approx = zeros(N_count,1);

r_error_approx_d = zeros(N_count,1);
v_error_approx_d = zeros(N_count,1);

r_opt_error = zeros(N_count,1);
v_opt_error = zeros(N_count,1);
AN_MARS_test = linspace(10,30, N_count);
d_coef = 1;
for i = 1:N_count 
    disp([num2str(i), '/', num2str(N_count)])
    %интегрируем
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   

    t0 = t_start;
    AN_i = AN_MARS_test(i);

    B=0.2721831;
    a_rel=1.52;
    %T_i = 2*pi*(AN_i*(a_rel-a_rel*B*log(a_rel))+a_rel*B^2*log(a_rel));
    T_i = 2*pi*(AN_i*1.35036912+0.02855872);
    OM_i = 2*pi*(0.27941185*AN_i-0.01523959);


    %НЕНУЛЕВЫЕ D
    [pr_0,pv_0] = get_initial_adjoint(AN_i,a_rel,d_coef);
    tspan = linspace(t0,t0+T_i, 10);
    %y0 = cat(2,start_pos,start_vel,zeros([1, 3]),zeros([1, 3]))';

    [a,eMag,iMag,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(start_pos',start_vel',mug);
    [p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,iMag,O,lonPer,nu);
    X = [p;ex;ey;ix;iy;L];
    dxdX = equitoctial2decart_jacobian(X,mug);
    P = dxdX'*[pr_0;pv_0];
    y0 = [X;P];
    tspan = linspace(t0,t0+T_i, AN_i*4);
    options = odeset('AbsTol',1e-10);
    options = odeset(options,'RelTol',1e-10);   

    %options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopL(s, y, AN_i*2*pi));
    [t,y] = ode113(@(t,y) equitoctial_adjoint_for_integration(y,mug), tspan,y0,options);

    traj = arrayfun(@(p,ex,ey,ix,iy,L) equitoctial2decart([p;ex;ey;ix;iy;L], mug),...
        y(:, 1),y(:, 2),y(:, 3),y(:, 4),y(:, 5),y(:, 6),'UniformOutput',false);
    traj = cell2mat(traj')';

    st.t = [t_start, T_i*T_unit/(3600*24)];
    st.planet = planet_end;
    st.mode = orbits;
    st.delta_omega = OM_i;

    [mars_r_f, mars_v_f]=planetModel(st);
    mars_r_f=mars_r_f'*1e+03/r_unit;
    mars_v_f=mars_v_f'*1e+03/V_unit;
    res_r_f = norm(traj(end, 1:3)-mars_r_f');
    res_v_f = norm(traj(end, 4:6)-mars_v_f');
    %res_v_f = norm(y(end, 4:6)-mars_v_f);
    r_error_approx(i) = res_r_f;
    v_error_approx(i) = res_v_f;

    %НУЛЕВЫЕ D
    [pr_0,pv_0] = get_initial_adjoint(AN_i,a_rel,0);
    tspan = linspace(t0,t0+T_i, 10);
    %y0 = cat(2,start_pos,start_vel,zeros([1, 3]),zeros([1, 3]))';

    [a,eMag,iMag,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(start_pos',start_vel',mug);
    [p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,iMag,O,lonPer,nu);
    X = [p;ex;ey;ix;iy;L];
    dxdX = equitoctial2decart_jacobian(X,mug);
    P = dxdX'*[pr_0;pv_0];
    y0 = [X;P];
    tspan = linspace(t0,t0+T_i, AN_i*4);
    options = odeset('AbsTol',1e-10);
    options = odeset(options,'RelTol',1e-10);   

    %options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopL(s, y, AN_i*2*pi));
    [t,y] = ode113(@(t,y) equitoctial_adjoint_for_integration(y,mug), tspan,y0,options);

    traj = arrayfun(@(p,ex,ey,ix,iy,L) equitoctial2decart([p;ex;ey;ix;iy;L], mug),...
        y(:, 1),y(:, 2),y(:, 3),y(:, 4),y(:, 5),y(:, 6),'UniformOutput',false);
    traj = cell2mat(traj')';

    st.t = [t_start, T_i*T_unit/(3600*24)];
    st.planet = planet_end;
    st.mode = orbits;
    st.delta_omega = OM_i;

    [mars_r_f, mars_v_f]=planetModel(st);
    mars_r_f=mars_r_f'*1e+03/r_unit;
    mars_v_f=mars_v_f'*1e+03/V_unit;
    res_r_f = norm(traj(end, 1:3)-mars_r_f');
    res_v_f = norm(traj(end, 4:6)-mars_v_f');
    %res_v_f = norm(y(end, 4:6)-mars_v_f);
    r_error_approx_d(i) = res_r_f;
    v_error_approx_d(i) = res_v_f;

    %ТОЧНОЕ РЕШЕНИЕ

    % yf = [mars_r_f;mars_v_f]';
    % 
    % opt_type = 'fsolve';
    % if strcmp(opt_type,'fsolve')
    %     fsolve_traj_fun=@(z)fsolve_traj(z,y0,yf,t0,T_i);
    %     options = optimoptions('fsolve','Display','off');
    %     %options = optimoptions(options,'UseParallel', true);
    %     %options = optimoptions(options, 'OptimalityTolerance', 1e-10);
    %     %options = optimoptions(options, 'MaxFunctionEvaluations', 1e+10);
    %     %options = optimoptions(options, 'StepTolerance', 1e-10);
    %     %options = optimoptions(options, 'Algorithm', 'levenberg-marquardt');
    %     [P0, fval] = fsolve(fsolve_traj_fun, P, options);
    %     r_opt_error(i) = norm(fval(1:3));
    %     v_opt_error(i) = norm(fval(4:6));
    % end
end
%%
figure(1)
plot(AN_MARS_test, v_error_approx_d,'-', 'DisplayName', 'velocity residual', 'LineWidth',1.0)
hold on;
%plot(AN_MARS_test, v_error_approx,'--', 'DisplayName', 'approximation with D coefficients', 'LineWidth',1.0)
%plot(AN_MARS_test, v_opt_error, 'DisplayName', 'optimization', 'LineWidth',1.0)
plot(AN_MARS_test, r_error_approx_d,'-', 'DisplayName', 'position residual', 'LineWidth',1.0)
hold off
legend('Location','southeast');
grid;
xlabel('Angular distance, revolutions')
ylabel('Residual of velocity, dimensionless')
set(gca, 'YScale', 'log')
%%
figure(2)
plot(AN_MARS_test, r_error_approx_d,'-', 'DisplayName', 'approximation without D coefficients', 'LineWidth',1.0)
hold on;
%plot(AN_MARS_test, r_error_approx,'--', 'DisplayName', 'approximation with D coefficients', 'LineWidth',1.0)
%plot(AN_MARS_test, r_opt_error, 'DisplayName', 'optimization', 'LineWidth',1.0)
hold off
%legend('Location','southeast');
grid;
xlabel('Angular distance, revolutions')
ylabel('Residual of position, dimensionless')
set(gca, 'YScale', 'log')

%% перебираем разное время и дату старта
a_rel = 1.52;
d_coef =1;
orbits='Ephemeris';
savefilename = 'ThesisMarsCheck_naive_apply_equinoctial_1.mat';
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
%%
load(savefilename)
for i = 100:N_count% 
    disp([num2str(i), '/', num2str(N_count)])
    for j = 1:N_count_2 
        if AN_MARS_2_test(i,j)>0
            continue
        end
        disp(['    ',num2str(j), '/', num2str(N_count_2)])
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
        [pr_0,pv_0] = get_initial_adjoint(AN_mars,a_rel,d_coef); 
        pr_0 = R^(-1)*pr_0;
        pv_0 = R^(-1)*pv_0; 
        tspan = linspace(0,t_end*2*pi, round(AN_mars*1000));
        
        [a,eMag,iMag,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(start_pos',start_vel',mug);
        [p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,iMag,O,lonPer,nu);
        X = [p;ex;ey;ix;iy;L];
        dxdX = equitoctial2decart_jacobian(X,mug);
        P = dxdX'*[pr_0;pv_0];
        y0 = [X;P];
        options = odeset('AbsTol',1e-10);
        options = odeset(options,'RelTol',1e-10);   
        
        %options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopL(s, y, AN_i*2*pi));
        [t,y] = ode113(@(t,y) equitoctial_adjoint_for_integration(y,mug), tspan,y0,options);
        
        traj = arrayfun(@(p,ex,ey,ix,iy,L) equitoctial2decart([p;ex;ey;ix;iy;L], mug),...
            y(:, 1),y(:, 2),y(:, 3),y(:, 4),y(:, 5),y(:, 6),'UniformOutput',false);
        traj = cell2mat(traj')';
        
        st.mode = orbits;
        st.t = t_start_MARS_test(j)+t_end*T_earth_days;
        %st.t = [t_start, t(end)*T_unit/(3600*24)];
        st.planet = planet_end;
        st.mode = orbits;
        
        [mars_r_f, mars_v_f]=planetModel(st);
        mars_r_f=mars_r_f'*1e+03/r_unit;
        mars_v_f=mars_v_f'*1e+03/V_unit;
        
        res_r_f = norm(traj(end, 1:3)-mars_r_f');
        res_v_f = norm(traj(end, 4:6)-mars_v_f');
        yf = [mars_r_f;mars_v_f]';
        r_error_approx(i,j) = res_r_f;
        v_error_approx(i,j) = res_v_f;
        
        mars_r_f_new = quatrotate(q_eps',mars_r_f');
        mars_an_f = atan2(mars_r_f_new(2), mars_r_f_new(1));
        traj_new = quatrotate(q_eps',traj(:,1:3));
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
        [pr_0,pv_0] = get_initial_adjoint(AN_mars,a_rel,d_coef); 
        pr_0 = R^(-1)*pr_0;
        pv_0 = R^(-1)*pv_0; 
        P = dxdX'*[pr_0;pv_0];
        y0 = [X;P];
        opt_type = 'fsolve';
        if strcmp(opt_type,'fsolve')
            fsolve_traj_fun=@(z)fsolve_traj(z,y0,yf,0,t_end*2*pi);
            options = optimoptions('fsolve','Display','off');
            %options = optimoptions(options,'UseParallel', true);
            %options = optimoptions(options, 'OptimalityTolerance', 1e-10);
            options = optimoptions(options, 'MaxFunctionEvaluations', 1e+4);
            %options = optimoptions(options, 'StepTolerance', 1e-10);
            %options = optimoptions(options, 'Algorithm', 'levenberg-marquardt');
            [P0, fval] = fsolve(fsolve_traj_fun, P, options);
            r_error_opt(i,j) = norm(fval(1:3));
            v_error_opt(i,j) = norm(fval(4:6));
            y0(7:12) = P0;
            %options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopL(s, y, AN_i*2*pi));
            options = odeset('AbsTol',1e-12);
            options = odeset(options,'RelTol',1e-12);   
            [t,y] = ode113(@(t,y) equitoctial_adjoint_for_integration(y,mug), tspan,y0,options);
            
            traj = arrayfun(@(p,ex,ey,ix,iy,L) equitoctial2decart([p;ex;ey;ix;iy;L], mug),...
                y(:, 1),y(:, 2),y(:, 3),y(:, 4),y(:, 5),y(:, 6),'UniformOutput',false);
            traj = cell2mat(traj')';
            traj_new = quatrotate(q_eps',traj(:,1:3));
            traj_an_f = atan2(traj_new(end,2), traj_new(end,1));
            AN_MARS_2_test(i,j) = calculate_angular_distance(traj_new(:,1:3))/(2*pi);
        end
        save(savefilename);
    end
end
%%
%угловая дальность
AN_i = 40;
%интегрируем
t0 = t_start;
d_coef = 1;
a_rel = 1.52;
T_i = 2*pi*T_a(AN_i, a_rel);
OM_i = 2*pi*OM_a(AN_i, a_rel);

%конечная орбита
ix = 0;
iy = 0;
L  = 0;
ey = 0;
ex = 0.0;
p = 1-ex^2;
X_2 = equitoctial2decart([p;ex;ey;ix;iy;L], mug);
start_pos_2 = X_2(1:3)';
start_vel_2 = X_2(4:6)';
R_2 = calculateRotMatrix(start_pos_2,start_vel_2);

%конечная орбита
ix = 0;
iy = 0;
L  = 0;
ey = 0;
ex = 0;
p = 1-ex^2;
X_3 = equitoctial2decart([p;ex;ey;ix;iy;L], mug);
start_pos_3 = X_3(1:3)';
start_vel_3 = X_3(4:6)';
R_3 = calculateRotMatrix(start_pos_3,start_vel_3);

% start_vel_2(1) = 0.0;
%start_vel_2(2) = 1;
%start_pos_2(1) = 1;

% [a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(start_pos_2',start_vel_2',mug);
% [p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
% disp(['ex=', num2str(ex)])
% disp(['ey=', num2str(ey)])
% disp(['p=', num2str(p)])
% disp(['rx=', num2str(start_pos_2(1))])
% disp(['--------'])

%R_2(1:2,1:2)=R_2(1:2,1:2)*0.3;

% n_rot = [0; 0; 1];


%поворот матрицы R_2
% delta_phi = 0*pi/24;
% R_rot = [cos(delta_phi),-sin(delta_phi),0;
%         sin(delta_phi),cos(delta_phi),0;
%         0,0,1];
% R_2=R_rot*R_2*R_rot';



% q_rot = [cos(delta_phi/2); n_rot*sin(delta_phi/2)];
% R_rot = quat2rotm(q_rot');
% R_rot = R_rot./sqrt(abs(R_rot));
% R_rot(isnan(R_rot)) = 0;
% R_2
% R_2(1,2) = 0;
% R_2(2,1) = 0;
%R_2(1,1) = 0.535;
%R_2(2,2) = 0.45;
%начальная орбита
%start_vel_3 = start_vel;
%start_pos_3 = start_pos;
%start_vel_3(1) = 0.;
%R_3 = calculateRotMatrix(start_pos_2,start_vel_2);

[pr_0,pv_0] = get_initial_adjoint(AN_i,a_rel,d_coef);
% pr_0_new = R_3^(-1)*R_2^(-1)*pr_0;
% pv_0_new  = R_3^(-1)*R_2^(-1)*pv_0;
pr_0_new = R_2^(-1)*R_3^(-1)*pr_0;
pv_0_new  = R_2^(-1)*R_3^(-1)*pv_0;
% pr_0_new = [-pr_0(1);pv_0(1);0];
% pv_0_new = [pr_0(2);-pv_0(2);0];
%y0 = cat(2,start_pos,start_vel,pr_0,pv_0)';
%подменяем на равноденственные
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(start_pos_3',start_vel_3',mug);

[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
X = [p;ex;ey;ix;iy;L];
dxdX = equitoctial2decart_jacobian(X,mug);
P = dxdX'*[pr_0_new;pv_0_new];
y0 = [X;P];
tspan = linspace(t0,t0+T_i, AN_i*400);
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);   

%options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopL(s, y, AN_i*2*pi));
[t,y] = ode113(@(t,y) equitoctial_adjoint_for_integration(y,mug), tspan,y0,options);

traj = arrayfun(@(p,ex,ey,ix,iy,L) equitoctial2decart([p;ex;ey;ix;iy;L], mug),...
    y(:, 1),y(:, 2),y(:, 3),y(:, 4),y(:, 5),y(:, 6),'UniformOutput',false);
traj = cell2mat(traj')';

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
y0 = cat(2,traj(1, 1:3), traj(1, 4:6),[0,0,0],[0,0,0])';
[t,orbit_initial] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
y0 = cat(2,traj(end, 1:3),traj(end, 4:6),[0,0,0],[0,0,0])';
[t,orbit_final] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);


%выводим график
figure(5);
plot3(0, 0, 0, 'k--o');
set(gca,'FontSize',14);
hold on;

plot3(traj(1, 1), traj(1, 2), traj(1, 3), 'O', 'LineWidth', 1);
plot3(traj(:, 1), traj(:, 2), traj(:, 3), 'cyan', 'LineWidth', 1);
plot3(orbit_initial(:, 1), orbit_initial(:, 2), orbit_initial(:, 3), 'k');
plot3(orbit_final(:, 1), orbit_final(:, 2), orbit_final(:, 3), 'r');
%plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k')
%plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'k')
%plot3(1.52*[0,cos(delta_phi)],1.52*[0,sin(delta_phi)],[0,0], 'g', 'LineWidth',2.0)
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(traj(end, 1:3)',traj(end, 4:6)',1);
[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
%plot3(10*[0,-ex],10*[0,-ey],[0,0], 'b--', 'LineWidth',2.0)
disp(['ex=', num2str(ex)])
disp(['ey=', num2str(ey)])
% disp(['e=', num2str(eMag)])
disp(['a=', num2str(a)])
disp(['p=', num2str(p)])
%disp(['i=', num2str(i)])
plot3(traj(end, 1), traj(end, 2), traj(end, 3),'bO')

%plot3(mars_r_f(1), mars_r_f(2), mars_r_f(3),'rO')

hold off;
axis equal

%title('Траектория КА')
xlabel('x, a.e.')
ylabel('y, a.e.')
zlabel('z, a.e.')
view(0,90)
%view(90,0)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
grid on;
box off;
function [value, isterminal, direction] = eventIntegrationTrajStopL(~, y, L_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

L=y(6);

value = L - L_end;
isterminal = 1;
direction = 0;
end
function R = calculateRotMatrix(r0,v0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
e1 = r0;
e2 = v0;
e3 = cross(e1,e2);
R = [e1;e2;e3];
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
res = traj(end,1:6)-yf(1:6);
end