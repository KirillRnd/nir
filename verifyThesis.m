%внешние константы
ae = 149597870700;
mug_0 = 132712.43994*(10^6)*(10^(3*3));
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;
T_earth_days = 365.256363004;
T_mars_days = 365.256363004*1.8808476;

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

% pr_x = @(a)(0.20074211987250948)/(7.780825627926492*a-0.34159993346898454+sin(2*pi*a+0.1231174475945867+0.5853586584957071/(a^2-1.5690074915306857*a+1.9819223815154108)) + (-41.26470185582318+86.35986725318767*a-51.08371071757919*a^2)*exp(-3.424188046763629*a) ); %OK
% pr_y = @(a)(0.0015463527047873246 - 0.0031678016585249486 * cos(2*pi * a+(-2.7056930690628103*a^2+5.176300671261158*a-2.336378597152512)*exp(3.2003639491958467*a-3.2205106566099353*a^2)))/(a^2-0.18328217249773662*a+0.11514749689845788);%OK
% pv_x = @(a)(0.0018661818762522732 - 0.0031537940448213170 * cos(2*pi * a+(-7.277086750486235*a^2+13.74021831427378*a-6.262552219285274)*exp(1.381195073394106*a-2.278956114228774*a^2)))/(a^2-0.19152227381595716*a+0.0993237695309777);
% pv_y = @(a)(0.10456877009148134)/(4.050922375438426*a-0.15902042687671913+sin(2*pi*a+0.11429217276196703+0.012218372217935064/(a^2-2.8347986183365705*a+2.06664561825753)-0.0019877228343783706*a^2) +(-540.5981059117084*a^2+958.0032106019531*a-415.15730979365355)*exp(-5.877538308252977*a));


pr_x = @(a,d)(0.20074211987250948)/(7.780825627926492*a-0.34159993346898454+sin(2*pi*a+0.1231174475945867+d*0.5853586584957071/(a^2-1.5690074915306857*a+1.9819223815154108)) + d*(-41.26470185582318+86.35986725318767*a-51.08371071757919*a^2)*exp(-3.424188046763629*a) ); %OK
pr_y = @(a,d)(0.0015463527047873246 - 0.0031678016585249486 * cos(2*pi * a+d*(-2.7056930690628103*a^2+5.176300671261158*a-2.336378597152512)*exp(3.2003639491958467*a-3.2205106566099353*a^2)))/(a^2-0.18328217249773662*a+0.11514749689845788);%OK
pv_x = @(a,d)(0.0018661818762522732 - 0.0031537940448213170 * cos(2*pi * a+d*(-7.277086750486235*a^2+13.74021831427378*a-6.262552219285274)*exp(1.381195073394106*a-2.278956114228774*a^2)))/(a^2-0.19152227381595716*a+0.0993237695309777);
pv_y = @(a,d)(0.10456877009148134)/(4.050922375438426*a-0.15902042687671913+sin(2*pi*a+0.11429217276196703+d*0.012218372217935064/(a^2-2.8347986183365705*a+2.06664561825753)) +d*(-540.5981059117084*a^2+958.0032106019531*a-415.15730979365355)*exp(-5.877538308252977*a));

B=0.2721831;
A=0.5433279;
T_a  = @(a, a_rel)(a_rel-a_rel*B*log(a_rel))*a+(a_rel*B^2*log(a_rel));
a_OM = @(OM, a_rel,k)(OM+k+A*exp(-a_rel)*log(a_rel)*log(a_rel))/(A*(1+exp(-a_rel))*log(a_rel));
OM_a = @(a, a_rel)A*(a*(exp(a_rel)+1)-log(a_rel))*exp(-a_rel)*log(a_rel);
%%
% PVevery_a_MARS_true = PVevery_a(:,113,:); %a = 1.52, Mars
% PRevery_a_MARS_true = PRevery_a(:,113,:); %a = 1.52, Mars
Tevery_a_MARS_true = Tevery_a(:,113);
Jevery_a_MARS_true = Jevery_a(:,113);
ANevery_a_MARS_true = ANevery_a(:,113);
% OMevery_a_MARS_true = OMevery_a(:,113);
% pr_x = @(a)(0.169249523878027)/(6.62725658086578*a-0.5425599+sin(6.40040102352653*a+0.380350602447347/a)); %OK
% pr_y = @(a)(0.00154886836922945-0.00345254159469051*cos(6.27309592045743*(a-0.9919472999999996)))/(a^2);%OK
% pv_x = @(a)(0.00193992356709378-0.0034504951248164*cos(6.26844568498378*a))/(a^2);%OK
% pv_y = @(a)(0.02592134)/(a+0.2475811*sin(6.2468360277011*a+0.2475811)-0.03856234084404);


%AN_MARS_test = linspace(0.75,6, 169);
PVevery_a_MARS_approx = zeros([169,1,3]); %a = 1.52, Mars
PRevery_a_MARS_approx = zeros([169,1,3]); %a = 1.52, Mars

for i = 1:169
    %AN_i = ANevery_a_MARS_true(i)/(2*pi);
    B=0.2721831;
    a_rel=1.52;
    T_i = 3600*24*Tevery_a_MARS_true(i)/(T_unit*2*pi);
    %AN_i = (T_i-0.03939635)/1.34784702;
    AN_i = (T_i-a_rel*B^2*log(a_rel))/(a_rel-a_rel*B*log(a_rel));
    %AN_i = AN_MARS_test(i);
    AN_i = ANevery_a_MARS_true(i)/(2*pi);
    d_coef = 1;
    PRevery_a_MARS_approx(i,1,:) = [pr_x(AN_i,d_coef),pr_y(AN_i,d_coef),0];
    PVevery_a_MARS_approx(i,1,:) = [pv_x(AN_i,d_coef),pv_y(AN_i,d_coef),0];
end

PVevery_a_MARS_approx_noD = zeros([169,1,3]); %a = 1.52, Mars
PRevery_a_MARS_approx_noD = zeros([169,1,3]); %a = 1.52, Mars

for i = 1:169
    %AN_i = ANevery_a_MARS_true(i)/(2*pi);
    B=0.2721831;
    a_rel=1.52;
    T_i = 3600*24*Tevery_a_MARS_true(i)/(T_unit*2*pi);
    %AN_i = (T_i-0.03939635)/1.34784702;
    AN_i = (T_i-a_rel*B^2*log(a_rel))/(a_rel-a_rel*B*log(a_rel));
    %AN_i = AN_MARS_test(i);
    AN_i = ANevery_a_MARS_true(i)/(2*pi);
    d_coef = 0;
    PRevery_a_MARS_approx_noD(i,1,:) = [pr_x(AN_i,d_coef),pr_y(AN_i,d_coef),0];
    PVevery_a_MARS_approx_noD(i,1,:) = [pv_x(AN_i,d_coef),pv_y(AN_i,d_coef),0];
end
%%
figure(1)
J_unit = 176.62545106129272;
J_ak = @(a,k)(log(k)*(4.747336E-3+1.044802E-2./k-1.517055E-2/(k^2)))./(a+2.020210E-1*(sin(2*pi*a+2.236760-1./a)-a).*exp(-a));
Jevery_a_MARS_pred = J_ak(ANevery_a_MARS_true/(2*pi),1.52);
plot(ANevery_a_MARS_true/(2*pi), Jevery_a_MARS_true/J_unit, 'LineWidth', 1.0, 'DisplayName', 'target functional')
hold on;
plot(ANevery_a_MARS_true/(2*pi), Jevery_a_MARS_pred, '--', 'LineWidth', 1.0, 'DisplayName', 'approximation')
hold off;

grid;
legend;

xlabel('Angular distance, revolutions')
ylabel('Functional, dimensionless')

figure(4)
plot(ANevery_a_MARS_true/(2*pi), Jevery_a_MARS_pred-Jevery_a_MARS_true/J_unit, 'LineWidth', 1.0, 'DisplayName', 'target functional')

xlabel('Angular distance, revolutions')
ylabel('Approximation error')
grid;
%%
figure(2)
plot(ANevery_a_MARS_true/(2*pi), PRevery_a_MARS_true(:,1,1),'LineWidth', 1.0, 'DisplayName', 'p_r_x numerical optimization')
hold on;
plot(ANevery_a_MARS_true/(2*pi), PRevery_a_MARS_true(:,1,2),'LineWidth', 1.0, 'DisplayName', 'p_r_y numerical optimization')
plot(ANevery_a_MARS_true/(2*pi), PVevery_a_MARS_true(:,1,1),'LineWidth', 1.0, 'DisplayName', 'p_v_x numerical optimization')
plot(ANevery_a_MARS_true/(2*pi), PVevery_a_MARS_true(:,1,2),'LineWidth', 1.0, 'DisplayName', 'p_v_y numerical optimization')
plot(ANevery_a_MARS_true/(2*pi), PRevery_a_MARS_approx(:,1,1),'--','LineWidth', 1.0, 'DisplayName', 'p_r_x approximation')

plot(ANevery_a_MARS_true/(2*pi), PRevery_a_MARS_approx(:,1,2),'--','LineWidth', 1.0, 'DisplayName', 'p_r_y approximation')

plot(ANevery_a_MARS_true/(2*pi), PVevery_a_MARS_approx(:,1,1),'--','LineWidth', 1.0, 'DisplayName', 'p_v_x approximation')

plot(ANevery_a_MARS_true/(2*pi), PVevery_a_MARS_approx(:,1,2),'--','LineWidth', 1.0, 'DisplayName', 'p_v_y approximation')
hold off;

legend;
grid;
xlabel('Angular distance, revolutions')
ylabel('Initial value of adjoint variable, dimensionless')

% figure(3)
% plot(ANevery_a_MARS_true/(2*pi), PRevery_a_MARS_true(:,1,1), 'DisplayName', 'p_r_x orig')
% hold on;
% plot(ANevery_a_MARS_true/(2*pi), PRevery_a_MARS_approx(:,1,1),'--', 'DisplayName', 'p_r_x approx')
% 
% plot(ANevery_a_MARS_true/(2*pi), PRevery_a_MARS_true(:,1,2), 'DisplayName', 'p_r_y orig')
% plot(ANevery_a_MARS_true/(2*pi), PRevery_a_MARS_approx(:,1,2),'--', 'DisplayName', 'p_r_y approx')
% hold off;
% 
% legend;
% grid;
% xlabel('Угловая дальность, витки.')

figure(5)
plot(ANevery_a_MARS_true/(2*pi), PRevery_a_MARS_approx(:,1,1)-PRevery_a_MARS_true(:,1,1),'LineWidth', 1.0, 'DisplayName', 'p_r_x error')
hold on;

plot(ANevery_a_MARS_true/(2*pi), PRevery_a_MARS_approx(:,1,2)-PRevery_a_MARS_true(:,1,2),'LineWidth', 1.0, 'DisplayName', 'p_r_y error')
plot(ANevery_a_MARS_true/(2*pi), PVevery_a_MARS_approx(:,1,1)-PVevery_a_MARS_true(:,1,1),'LineWidth', 1.0, 'DisplayName', 'p_v_x error')
plot(ANevery_a_MARS_true/(2*pi), PVevery_a_MARS_approx(:,1,2)-PVevery_a_MARS_true(:,1,2),'LineWidth', 1.0, 'DisplayName', 'p_v_y error')
hold off;

legend;
grid;
xlabel('Angular distance, revolutions')
ylabel('Approximation error, dimensionless ')
disp('Errors');
disp(sqrt(mse(PRevery_a_MARS_approx(:,1,1),PRevery_a_MARS_true(:,1,1))))
disp(sqrt(mse(PRevery_a_MARS_approx(:,1,2),PRevery_a_MARS_true(:,1,2))))
disp(sqrt(mse(PVevery_a_MARS_approx(:,1,1),PVevery_a_MARS_true(:,1,1))))
disp(sqrt(mse(PVevery_a_MARS_approx(:,1,2),PVevery_a_MARS_true(:,1,2))))
%% перебираем разные дальности
r_error_orig = zeros(L2,1);
v_error_orig = zeros(L2,1);
for i = 1:169 
    %интегрируем
    options = odeset('AbsTol',1e-10);
    options = odeset(options,'RelTol',1e-10);   
    
    t0 = t_start;
    tspan = linspace(t0,t0+3600*24*Tevery_a_MARS_true(i)/T_unit, 1000);
    %y0 = cat(2,start_pos,start_vel,zeros([1, 3]),zeros([1, 3]))';
    pr_0 = reshape(PRevery_a_MARS_true(i,1,:), [1, 3]);
    pv_0 = reshape(PVevery_a_MARS_true(i,1,:), [1, 3]);
    y0 = cat(2,start_pos,start_vel,pr_0,pv_0)';
    
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    st.t = [t_start, Tevery_a_MARS_true(i)];
    st.planet = planet_end;
    st.mode = orbits;
    st.delta_omega = OMevery_a_MARS_true(i);
    
    [mars_r_f, mars_v_f]=planetModel(st);
    mars_r_f=mars_r_f'*1e+03/r_unit;
    mars_v_f=mars_v_f'*1e+03/V_unit;
    res_r_f = norm(y(end, 1:3)-mars_r_f');
    res_v_f = norm(y(end, 4:6)-mars_v_f');

    r_error_orig(i) = res_r_f;
    v_error_orig(i) = res_v_f;
end
%% перебираем разные дальности с фиксироваными T и OM
r_error_approx = zeros(L2,1);
v_error_approx = zeros(L2,1);
r_error_approx_noD = zeros(L2,1);
v_error_approx_noD = zeros(L2,1);
r_opt_error = zeros(L2,1);
v_opt_error = zeros(L2,1);
for i = 1:169 
    %интегрируем
    options = odeset('AbsTol',1e-10);
    options = odeset(options,'RelTol',1e-10);   
    
    t0 = t_start;
    
    t_end = 3600*24*Tevery_a_MARS_true(i)/T_unit;
    tspan = linspace(t0,t0+t_end, 1000);
    %y0 = cat(2,start_pos,start_vel,zeros([1, 3]),zeros([1, 3]))';

    pr_0 = reshape(PRevery_a_MARS_approx(i,1,:), [1, 3]);
    %pr_0(1) = PRevery_a_MARS_true(i,1,1);
    %pr_0(2) = PRevery_a_MARS_true(i,1,2);
    pv_0 = reshape(PVevery_a_MARS_approx(i,1,:), [1, 3]);
    %pv_0(1) = PVevery_a_MARS_true(i,1,1);
    %pv_0(2) = PVevery_a_MARS_true(i,1,2);
    y0 = cat(2,start_pos,start_vel,pr_0,pv_0)';
    
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    st.t = [t_start, Tevery_a_MARS_true(i)];
    st.planet = planet_end;
    st.mode = orbits;
    st.delta_omega = OMevery_a_MARS_true(i);
    
    [mars_r_f, mars_v_f]=planetModel(st);
    mars_r_f=mars_r_f'*1e+03/r_unit;
    mars_v_f=mars_v_f'*1e+03/V_unit;
    res_r_f = norm(y(end, 1:3)-mars_r_f');
    res_v_f = norm(y(end, 4:6)-mars_v_f');

    r_error_approx(i) = res_r_f;
    v_error_approx(i) = res_v_f;

    pr_0 = reshape(PRevery_a_MARS_approx_noD(i,1,:), [1, 3]);
    %pr_0(1) = PRevery_a_MARS_true(i,1,1);
    %pr_0(2) = PRevery_a_MARS_true(i,1,2);
    pv_0 = reshape(PVevery_a_MARS_approx_noD(i,1,:), [1, 3]);
    %pv_0(1) = PVevery_a_MARS_true(i,1,1);
    %pv_0(2) = PVevery_a_MARS_true(i,1,2);
    y0 = cat(2,start_pos,start_vel,pr_0,pv_0)';
    
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    st.t = [t_start, Tevery_a_MARS_true(i)];
    st.planet = planet_end;
    st.mode = orbits;
    st.delta_omega = OMevery_a_MARS_true(i);
    
    [mars_r_f, mars_v_f]=planetModel(st);
    mars_r_f=mars_r_f'*1e+03/r_unit;
    mars_v_f=mars_v_f'*1e+03/V_unit;
    res_r_f = norm(y(end, 1:3)-mars_r_f');
    res_v_f = norm(y(end, 4:6)-mars_v_f');

    r_error_approx_noD(i) = res_r_f;
    v_error_approx_noD(i) = res_v_f;

    % yf = [mars_r_f;mars_v_f]';
    % opt_type = 'fsolve';
    % p0 = [pr_0,pv_0]';
    % if strcmp(opt_type,'fsolve')
    % 
    %     fsolve_traj_fun=@(z)fsolve_traj(z,y0,yf,t0,t_end);
    %     options = optimoptions('fsolve','Display','off');
    %     options = optimoptions(options,'UseParallel', true);
    %     options = optimoptions(options, 'OptimalityTolerance', 1e-10);
    %     options = optimoptions(options, 'MaxFunctionEvaluations', 1e+10);
    %     options = optimoptions(options, 'StepTolerance', 1e-10);
    %     options = optimoptions(options, 'Algorithm', 'levenberg-marquardt');
    % 
    %     p0 = fsolve(fsolve_traj_fun, p0, options);
    % end
    % 
    % y0 = cat(2,start_pos,start_vel,p0(1:3)',p0(4:6)')';
    % 
    % options = odeset('AbsTol',1e-10);
    % options = odeset(options,'RelTol',1e-10);   
    % %options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 1.52));
    % [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    % 
    % res_r_f = norm(y(end, 1:3)-mars_r_f');
    % res_v_f = norm(y(end, 4:6)-mars_v_f');
    % r_opt_error(i) = res_r_f;
    % v_opt_error(i) = res_v_f;
end
%%
figure(1)
plot(ANevery_a_MARS_true/(2*pi), r_error_orig, 'k', 'LineWidth', 1.0,'DisplayName', 'Numerical optimization')
hold on;
plot(ANevery_a_MARS_true/(2*pi), r_error_approx,'b--', 'LineWidth', 1.0,'DisplayName', 'Approximation with D coefficients')
plot(ANevery_a_MARS_true/(2*pi), r_error_approx_noD,'r-.','LineWidth', 1.0, 'DisplayName', 'Approximation without D coefficients')
%plot(ANevery_a_MARS_true/(2*pi), r_opt_error,'--', 'DisplayName', 'Невязка положения новое решение')
%plot(ANevery_a_MARS_true/(2*pi), v_opt_error,'--', 'DisplayName', 'Невязка скорости новое решение')
hold off
legend('Location','best');
grid;
xlabel('Angular distance, revolutions')
ylabel('Residual of position, dimensionless')
set(gca, 'YScale', 'log')
set(gca, 'FontSize', 11)

figure(2)
plot(ANevery_a_MARS_true/(2*pi), v_error_orig,'k', 'LineWidth', 1.0, 'DisplayName', 'Numerical optimization')
hold on;
plot(ANevery_a_MARS_true/(2*pi), v_error_approx,'b--','LineWidth', 1.0, 'DisplayName', 'Approximation with D coefficients')
plot(ANevery_a_MARS_true/(2*pi), v_error_approx_noD,'r-.','LineWidth', 1.0, 'DisplayName', 'Approximation without D coefficients')
%plot(ANevery_a_MARS_true/(2*pi), r_opt_error,'--', 'DisplayName', 'Невязка положения новое решение')
%plot(ANevery_a_MARS_true/(2*pi), v_opt_error,'--', 'DisplayName', 'Невязка скорости новое решение')
hold off
legend('Location','best');
grid;
xlabel('Angular distance, revolutions')
ylabel('Residual of velocity, dimensionless')
set(gca, 'YScale', 'log')
set(gca, 'FontSize', 11)
%% перебираем разные дальности
N_count = 197;
r_error_approx = zeros(N_count,1);
v_error_approx = zeros(N_count,1);

r_opt_error = zeros(N_count,1);
v_opt_error = zeros(N_count,1);
AN_MARS_test = linspace(1,50, N_count);
J_MARS_opt = zeros(N_count,1);
d_coef = 0;
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

    pr_0 = [pr_x(AN_i,d_coef),pr_y(AN_i,d_coef),0];
    pv_0 = [pv_x(AN_i,d_coef),pv_y(AN_i,d_coef),0];
    %t_end = 3600*24*Tevery_a_MARS_true(i)/T_unit;
    tspan = linspace(t0,t0+T_i, 10);
    %y0 = cat(2,start_pos,start_vel,zeros([1, 3]),zeros([1, 3]))';

    y0 = cat(2,start_pos,start_vel,pr_0,pv_0)';
    
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    %options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 1.52));
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    an_end = atan2(y(end,2),y(end,1));
    %mars_v_f = [cos(an_end + pi/2),sin(an_end + pi/2),0];

    st.t = [t_start, T_i*T_unit/(3600*24)];
    st.planet = planet_end;
    st.mode = orbits;
    st.delta_omega = OM_i;
    
    [mars_r_f, mars_v_f]=planetModel(st);
    mars_r_f=mars_r_f'*1e+03/r_unit;
    mars_v_f=mars_v_f'*1e+03/V_unit;
    res_r_f = norm(y(end, 1:3)-mars_r_f');
    res_v_f = norm(y(end, 4:6)-mars_v_f');
    %res_v_f = norm(y(end, 4:6)-mars_v_f);
    r_error_approx(i) = res_r_f;
    v_error_approx(i) = res_v_f;
    yf = [mars_r_f;mars_v_f]';
    opt_type = 'fsolve';
    p0 = [pr_0,pv_0]';
    if strcmp(opt_type,'fsolve')

        fsolve_traj_fun=@(z)fsolve_traj(z,y0,yf,t0,T_i);
        options = optimoptions('fsolve','Display','off');
        %options = optimoptions(options,'UseParallel', true);
        options = optimoptions(options, 'OptimalityTolerance', 1e-10);
        options = optimoptions(options, 'MaxFunctionEvaluations', 1e+10);
        options = optimoptions(options, 'StepTolerance', 1e-10);
        %options = optimoptions(options, 'Algorithm', 'levenberg-marquardt');
        p0 = fsolve(fsolve_traj_fun, p0, options);
    elseif strcmp(opt_type,'fmincon')
        fmincon_traj_fun=@(z)fmincon_traj(z,y0,yf,t0,T_i);
        options = optimoptions('fmincon','Display','off');
        options = optimoptions(options,'UseParallel', true);
        options = optimoptions(options, 'OptimalityTolerance', 1e-10);
        options = optimoptions(options, 'MaxFunctionEvaluations', 1e+10);
        options = optimoptions(options, 'StepTolerance', 1e-10);
        options = optimoptions(options, 'ConstraintTolerance', 1e-10);
        options = optimoptions(options, 'MaxIterations', 250);
        %options = optimoptions(options, 'FiniteDifferenceType', 'central');
        %options = optimoptions(options, 'EnableFeasibilityMode', true);
        p_limit = 10^(-2);
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        lb = p0-[1, 1, 1, 1, 1, 1]'*p_limit;
        ub = p0+[1, 1, 1, 1, 1, 1]'*p_limit;
        p0 = fmincon(@(x)1, p0, A, b, Aeq, beq, lb, ub, fmincon_traj_fun, options);
    end
    
    y0 = cat(2,start_pos,start_vel,p0(1:3)',p0(4:6)')';
    
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    %options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 1.52));
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    a_reactive = y(:,10:12);
    a_vec=vecnorm(a_reactive, 2, 2).^2;
    Jt = cumtrapz(t, a_vec)/(2);
    J_MARS_opt(i) = Jt(end);
    res_r_f = norm(y(end, 1:3)-mars_r_f');
    res_v_f = norm(y(end, 4:6)-mars_v_f');
    r_opt_error(i) = res_r_f;
    v_opt_error(i) = res_v_f;
end
figure(1)
plot(AN_MARS_test, r_error_approx,'--', 'DisplayName', 'Невязка положения аппроксимация')
hold on;
plot(AN_MARS_test, v_error_approx,'--', 'DisplayName', 'Невязка скорости аппроксимация')
plot(AN_MARS_test, r_opt_error, 'DisplayName', 'Невязка положения решение')
plot(AN_MARS_test, v_opt_error, 'DisplayName', 'Невязка скорости решение')
hold off
legend('Location','best');
grid;
xlabel('Угловая дальность, витки.')
ylabel('Невязка на правом конце')
set(gca, 'YScale', 'log')

%%
J_a  = @(a)(0.002117*exp(a))/(a*exp(a)+0.17856*(sin(7.44443*a)-a));
J_MARS_approx= zeros(N_count,1);
for i = 1:N_count 
    J_MARS_approx(i) = J_a(AN_MARS_test(i));
end
figure(2)
plot(AN_MARS_test, J_MARS_opt, 'DisplayName', 'Целевой функционал')
hold on;
plot(AN_MARS_test, J_MARS_approx,'--', 'DisplayName', 'Аппроксимированный функционал')
hold off
grid;
xlabel('Угловая дальность, витки.')
ylabel('Функционал')

%%

vecPV_cartesian =[

  -0.017828836330895
  -0.069754179573323
  -0.000000000000005];


vecPR_cartesian = [
  -0.083400848307302
  -0.014320246933389
  -0.000000000000004];
T_point =   3.138140048076477;
T_point_2 = 1.146426831543537e+03;


% T_point =   7.668293435690914;
% 
% vecPV_cartesian =[
%    0.000120968372462
%    0.004782528163732
%   -0.000000000000000];
% 
% vecPR_cartesian =[
%    0.004694925114471
%    0.000109525174757
%   -0.000000000000000];

%интегрируем
%options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 1.52));
t0 = t_start;
i = 5;
%y0 = cat(2,start_pos,start_vel,zeros([1, 3]),zeros([1, 3]))';
%pr_0 = reshape(PRevery_a_MARS_approx(i,1,:), [1, 3]);
%pv_0 = reshape(PVevery_a_MARS_approx(i,1,:), [1, 3]);
%AN_i = ANevery_a_MARS_true(i)/(2*pi);
AN_i = 10;

B=0.2721831;
a_rel=1.52;
%T_i = 2*pi*(AN_i*(a_rel-a_rel*B*log(a_rel))+a_rel*B^2*log(a_rel));
start_vel_2 = start_vel;
start_pos_2 = start_pos;
%start_vel_2(1) = 0.4;
%start_pos_2(1) = 1.2;
%start_vel_2(2) = sqrt(1/start_pos_2(1));
%a_rel = 0.387;
T_i = 2*pi*T_a(AN_i, a_rel);
%T_i = 2*pi*(AN_i*1.35036912+0.02855872);
%OM_i = 2*pi*(0.27941185*AN_i-0.01523959);
OM_i = 2*pi*OM_a(AN_i, a_rel);
%[Kr,Kv] = calculateKMatrix(a_rel);
K = calculateKMatrix_ver2(a_rel,AN_i,1);
%выбираем вектор поворота
n_rot = [1; 0; 0];
delta_i = 0*pi/180;
q_rot = [cos(delta_i/2); n_rot*sin(delta_i/2)];
R_2 = quat2rotm(q_rot');
R_2(1:3,3)=R_2(1:3,3)*R_2(1,1)*R_2(2,2);
R_2(1:3,2)=R_2(1:3,2)*R_2(2,2);
R_2(1:3,1)=R_2(1:3,1)*R_2(1,1);
% R_2(1,1)=1.013;
% R_2(2,2)=R_2(2,2)*1.013;
% R_2(3,3)=R_2(3,3)*1;
n_rot = [0; 0; 1];
delta_i = 0*pi/6;
q_rot = [cos(delta_i/2); n_rot*sin(delta_i/2)];
R_3 = quat2rotm(q_rot');
R_2=R_3*R_2*R_3';
R = calculateRotMatrix(start_pos_2,start_vel_2);
d_coef = 1;
pr_0 = K*R_2^(-1)*R^(-1)*[pr_x(AN_i,d_coef);pr_y(AN_i,d_coef);0.00];
pr_0 = pr_0';
pv_0 = K*R_2^(-1)*R^(-1)*[pv_x(AN_i,d_coef);pv_y(AN_i,d_coef);0.00]; 
pv_0 = pv_0';

pr_0 = vecPR_cartesian';
pv_0 = vecPV_cartesian';
T_i = T_point*2*pi;
y0 = cat(2,start_pos_2,start_vel_2,pr_0,pv_0)';
M = cross(pv_0, start_vel_2)+cross(pr_0, start_pos_2);
tspan = linspace(t0,t0+T_i, AN_i*400);
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
%options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopE0(s, y));

%%
ddeltady0=zeros([6,6]);
step_h=sqrt(eps);
for i=7:12
    y0_delta=zeros([12,1]);
    y0_delta(i)=step_h;
    
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0+y0_delta,options);
    p_plus=y(end,1:6);

    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0-y0_delta,options);
    p_minus=y(end,1:6);
    partial=(p_plus-p_minus)/(2*step_h);
    ddeltady0(:,i-6)=partial;
end

C=cond(ddeltady0);

%%
[t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);

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
%earth_traj=1e+3*earth_traj/ae;
earth_traj=earth_traj/ae;

t_orbit = linspace(t0,t0+T_mars/(24*3600), 1000);

st.t = t_orbit';
st.planet = planet_end;
st.mode = orbits;
st.delta_omega = 0;

mars_traj = planetModel(st);
%mars_traj=1e+3*mars_traj/ae;
mars_traj=mars_traj/ae;
y0 = cat(2,start_pos_2,start_vel_2,[0,0,0],[0,0,0])';
[t,orbit_initial] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
y0 = cat(2,y(end, 1:3),y(end, 4:6),[0,0,0],[0,0,0])';
[t,orbit_final] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);


%выводим график
figure(5);
plot3(0, 0, 0, 'k--o');
set(gca,'FontSize',14);
hold on;

plot3(y(:, 1), y(:, 2), y(:, 3), 'cyan', 'LineWidth', 2);
%plot3(orbit_initial(:, 1), orbit_initial(:, 2), orbit_initial(:, 3), 'k');
%plot3(orbit_final(:, 1), orbit_final(:, 2), orbit_final(:, 3), 'r');
plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k')
plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'k')
plot3(1.52*[cos(delta_i+pi),cos(delta_i)],1.52*[sin(delta_i+pi),sin(delta_i)],[0,0])
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(y(end, 1:3)',y(end, 4:6)',1);
disp(['e=', num2str(eMag)])
disp(['a=', num2str(a)])
disp(['i=', num2str(i)])
%plot3(y(end, 1), y(end, 2), y(end, 3),'bO')

%plot3(mars_r_f(1), mars_r_f(2), mars_r_f(3),'rO')

hold off;
axis equal

title('Траектория КА')
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
%%
function res = fsolve_traj(z,y0,yf,t0,t_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tspan = linspace(t0,t0+t_end, 10);
y0(7:12)=z;
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
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
function [value, isterminal, direction] = eventIntegrationTrajStopR(~, y, r_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

r=y(1:3);

value = norm(r) - r_end;
isterminal = 1;
direction = 0;
end
function [value, isterminal, direction] = eventIntegrationTrajStopE0(~, y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

r=y(1:3);
v=y(4:6);
e = getEccentricity(r,v,1);
value = e;
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
function [Kr,Kv] = calculateKMatrix(k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
a1 = (k+2.087217245873957)*log(k)/k;
a2 = -k*exp(-8.563932009313163/k)+69.18759116365857*exp(-8.563932009313163*k);
a3 = -0.07470468492737094*log(k)+(0.6602914816446814*log(k)*exp(-k))/(k^3/2);
a4 = (-0.8348043664990472*k+2.922882836448715*log(k)+3.6786032392739227)*log(k)/k;
a2=0;
a3=0;
a4=a1;
Kr = [[a1,a2,0];[a3,a4,0];[0,0,1]];
Kv = [[a4,a3,0];[a2,a1,0];[0,0,1]];
end
function K = calculateKMatrix_ver2(k_target,a,d)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%коррекция kappa
c0 = 1.568311392349845;
c1 = (1.52-1)/(exp(-c0/1.52)-exp(-c0));
c2 = 1 - c1*exp(-c0);
k = c2+c1*exp(-c0/k_target);
g = (-8.234450871689381-1.835384888679592*cos(2*pi*a+4.568045799135197)+12.140923790903319*a)/...
    (6.521910293172844*a-4.627546590315015+sin(2*pi*a+2.963335050753793+d*(0.010167241207198103)/(a^2-1.4200382540034495*a+0.5258798683545889))...
    + d*(-22.253725602121175*a^2+47.5730432773046*a-24.334701239930087)*exp(-3.0200712580477367*a) );
f = -g/1.52+1/log(1.52);
gamma = (f*k+g)*log(k)/k;

K = [[gamma,0,0];[0,gamma,0];[0,0,1]];
end