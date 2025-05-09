%внешние константы
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
t_start=juliandate(2026,1,1);
orbits='Simple';
N=1350;
m0=367;
eta=0.45;
planet_start = 'Earth';
planet_end = 'Mars';
omega = -pi;

st.t = t_start;
st.planet = planet_end;
st.mode = orbits;
st.delta_omega = 0;

[mars_pos, mars_vel] = planetModel(st);
mars_pos=1e+3*mars_pos/ae;
mars_vel=1e+3*mars_vel/V_unit;
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(mars_pos',mars_vel',1);
Omega_Mars = O*180/pi;
omega_Mars = lonPer*180/pi;
i_Mars = i*180/pi;
e_Mars = eMag;
a_Mars = a;

% Omega_Mars = 0;
% omega_Mars = 0;
% i_Mars = 0.01;
% e_Mars = 0;
% a_Mars = 1.3;

%определяем стартовое положение
st.t = t_start;
st.planet = planet_start;
st.mode = orbits;
st.delta_omega = omega;

[start_pos, start_vel] = planetModel(st);
start_pos=start_pos*1e+03/ae;
start_vel=start_vel*1e+03/V_unit;
%start_pos=[1 0 0];
%start_vel=[0 1 0];
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(start_pos',start_vel',1);

start_an = truLon*180/pi;
B=0.2721831;
A=0.5433279;
T_a  = @(a, a_rel)(a_rel-a_rel*B*log(a_rel))*a+(a_rel*B^2*log(a_rel));
a_OM = @(OM, a_rel,k)(OM+k+A*exp(-a_rel)*log(a_rel)*log(a_rel))/(A*(1+exp(-a_rel))*log(a_rel));
OM_a = @(a, a_rel)A*(a*(exp(a_rel)+1)-log(a_rel))*exp(-a_rel)*log(a_rel);

%%
%интегрируем
t0 = t_start;                 %дата старта, безразм
d_coef = 1;                   %учёт D-коэффициентов, 0 или 1
%OM_i = 2*pi*OM_a(AN_i, a_rel);%оптимальная фаза

%изменение начальной орбиты
start_pos_2 = start_pos;
start_vel_2 = start_vel;
%start_vel_2(1)=0.1;

%новая аппроксимация времени
%T_i = 2*pi*(AN_i*(0.03772064567906565*a_rel + 1.405055577366462)*log(a_rel + 1)+(0.08297934358955021*a_rel + 0.009838514316752341).*log(a_rel).^2);
%i_target_array = [i_Mars*0.1,i_Mars*0.2,i_Mars*0.3,i_Mars*0.4,i_Mars*0.5,i_Mars*0.6,i_Mars*0.7,i_Mars*0.8,i_Mars*0.9,i_Mars];%0:2:4;

O_initial = Omega_Mars-start_an;
o_initial = omega_Mars-start_an;

i_target_array = i_Mars;
%i_target_array = [0.95*i_Mars, i_Mars, 1.05*i_Mars];
O_target_array = O_initial;
e_target_array =  e_Mars;
%e_target_array =  [0.8*e_Mars, e_Mars, 1.2*e_Mars];
%o_target_array = [o_initial+10, o_initial, o_initial-10];
o_target_array = o_initial;

AN_array = 0.75:0.01:6;
%AN_array = 1:1:6;
k_array = 1.2:0.01:2.2;

N_count_1 = length(i_target_array);
N_count_2 = length(O_target_array);
N_count_3 = length(e_target_array);
N_count_4 = length(o_target_array);

N_count_5 = length(AN_array);
N_count_6 = length(k_array);

%rf_massive = zeros(N_count_1,N_count_2,N_count_3,N_count_4,N_count_5,N_count_6,3);
%vf_massive = zeros(N_count_1,N_count_2,N_count_3,N_count_4,N_count_5,N_count_6,3);
%J_massive = zeros(N_count_1,N_count_2,N_count_3,N_count_4,N_count_5,N_count_6,1);
total_iter = N_count_1 * N_count_2 * N_count_3 * N_count_4 * N_count_5 * N_count_6;


rf_all = zeros(total_iter, 3);
vf_all = zeros(total_iter, 3);
J_all  = zeros(total_iter, 1);
T_all  = zeros(total_iter, 1);

i_all = zeros(total_iter, 1);
O_all = zeros(total_iter, 1);
e_all  = zeros(total_iter, 1);
o_all  = zeros(total_iter, 1);
a_all  = zeros(total_iter, 1);
AN_all  = zeros(total_iter, 1);
parfor idx = 1:total_iter
    [j_1,j_2,j_3,j_4,j_5,j_6] = ind2sub([N_count_1,N_count_2,N_count_3,N_count_4,N_count_5,N_count_6], idx);
    i_target = pi*i_target_array(j_1)/180;
    O_target = pi*O_target_array(j_2)/180;
    e_target = e_target_array(j_3);
    o_target = pi*o_target_array(j_4)/180;

    AN_i = AN_array(j_5);         %угловая дальность, витки
    a_rel = k_array(j_6);         %соотношение радиусов
    T_i = 2*pi*T_a(AN_i, a_rel);  %время перелёта

    R = calculateRotMatrix(start_pos_2,start_vel_2);
    
    [pr_0,pv_0] = get_initial_adjoint(AN_i,a_rel,d_coef);
    R_total = R;
    p0 = [pr_0, pv_0];
    
    %[rf,vf,Jf,r_all] = integrate_Rtotal_equinoctial(R_total, start_pos_2,start_vel_2,p0,t0,T_i);
    [rf,vf,~,~] = integrate_Rtotal(R_total, start_pos_2,start_vel_2,p0,t0,T_i);
    [a_fix,e_fix,~,~,~,~,~,~,o_fix,~] = rv2orb(rf,vf,1);

    R_eo_2 = calculateReoMatrix(a_fix, e_fix,o_fix+pi-start_an/180*pi);
    R_total = R_eo_2*R;
    
    [rf,vf,~,~] = integrate_Rtotal(R_total, start_pos_2,start_vel_2,p0,t0,T_i);
    [a_fix,~,~,~,~,~,~,~,~,~] = rv2orb(rf,vf,1);
    R_eo = calculateReoMatrix(a_fix, e_target,o_target);
    R_iO = calculateRiOMatrix(a_fix, i_target,O_target);

    R_total = R_iO*R_eo*R_eo_2*R;
    %R_total = R_eo_2*R;
    [rf,vf,Jf,~] = integrate_Rtotal(R_total, start_pos_2,start_vel_2,p0,t0,T_i);

    %rf_massive(j_1,j_2,j_3,j_4,j_5,j_6,:) = rf;
    %vf_massive(j_1,j_2,j_3,j_4,j_5,j_6,:) = vf;
    %J_massive(j_1,j_2,j_3,j_4,j_5,j_6,:) = Jf;
    rf_all(idx,:) = rf;    % rf - это вектор размера 1x3
    vf_all(idx,:) = vf;    % vf - это вектор размера 1x3
    J_all(idx)    = Jf;    % Jf - скаляр    
    T_all(idx)    = T_i;    % Jf - скаляр    


    [a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(rf,vf,1);
    i_all(idx,:) = i;  
    O_all(idx,:) = O;    
    e_all(idx,:) = eMag;   
    o_all(idx,:) = lonPer;    
    a_all(idx,:) = a;   
    AN_all(idx,:) = AN_i;   
end
%%
T_earth_days = 365.256363004;
T_earth = T_earth_days*3600*24;
%коэффициенты обезразмеривания
r_unit=ae;
V_unit=sqrt(mug_0/ae);
T_unit_days = T_earth_days/(2*pi);
T_all_days = T_unit_days*T_all;
%%
diff_a = abs(a_all-a_Mars);
diff_i = abs(180*i_all/pi-i_Mars);
diff_e = abs(e_all-e_Mars);
ok_condition = diff_a < 0.01;
disp(sum(ok_condition));

i_all_ok = 180*i_all(ok_condition)/pi;
e_all_ok = e_all(ok_condition);
O_all_ok = 180*unwrap(O_all(ok_condition))/pi;
o_all_ok = 180*unwrap(o_all(ok_condition))/pi;

figure(1)
plot(AN_all(ok_condition),i_all_ok,'.')
hold on;
plot([AN_all(1), AN_all(end)],[i_Mars, i_Mars], 'LineWidth',2)
hold off;
title('i')
grid()
figure(2)
plot(AN_all(ok_condition),O_all_ok,'.')
hold on;
plot([AN_all(1), AN_all(end)],[Omega_Mars, Omega_Mars], 'LineWidth',2)
hold off;
title('O')
grid()
figure(3)
plot(AN_all(ok_condition),e_all_ok,'.')
hold on;
%plot(AN_all(ok_condition),e_target-0.25*e_target*sin(AN_all(ok_condition)*2*pi),'.')
plot([AN_all(1), AN_all(end)],[e_Mars, e_Mars], 'LineWidth',2)
hold off;
title('e')
grid()
figure(4)
plot(AN_all(ok_condition),o_all_ok,'.')
hold on;
plot([AN_all(1), AN_all(end)],[omega_Mars, omega_Mars],'LineWidth',2)
hold off;
title('o')
grid()
figure(5)
plot(AN_all(ok_condition),a_all(ok_condition),'.')
hold on;
plot([AN_all(1), AN_all(end)],[a_Mars, a_Mars],'LineWidth',2)
hold off;
title('a')
grid()
%%

i_all_ok_err = sqrt( mean( (i_all_ok - i_Mars).^2 ,'all') );
e_all_ok_err = sqrt( mean( (e_all_ok - e_Mars).^2 ,'all') );
O_all_ok_err = sqrt( mean( (O_all_ok - Omega_Mars).^2 ,'all') );
o_all_ok_err = sqrt( mean( (o_all_ok - omega_Mars).^2 ,'all') );

disp(['e RMSE=', num2str(e_all_ok_err)])
disp(['i RMSE=', num2str(i_all_ok_err)])
disp(['O RMSE=', num2str(O_all_ok_err)])
disp(['o RMSE=', num2str(o_all_ok_err)])
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

% Марс
t_orbit = linspace(t0,t0+T_mars/(24*3600), 1000);

st.t = t_start+T_all_days;
st.planet = planet_end;
st.mode = orbits;
st.delta_omega = 0;

[mars_traj, mars_vel] = planetModel(st);
mars_traj=1e+3*mars_traj/ae;
mars_vel=1e+3*mars_vel/V_unit;

diff_rf = vecnorm(rf_all-mars_traj, 2, 2);
diff_vf = vecnorm(vf_all-mars_vel, 2, 2);
good_condition = diff_rf < 0.1 & diff_vf < 0.1;
disp(sum(good_condition));

options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
y0 = cat(2,start_pos_2,start_vel_2,[0,0,0],[0,0,0])';
tspan_ini = linspace(t0,t0+2*pi, 400);
[t,orbit_initial] = ode113(@(t,y) internalIntegration3D(t,y), tspan_ini,y0,options);
%y0 = cat(2,rf',vf',[0,0,0],[0,0,0])';
%tspan_fin = linspace(t0,t0+2*pi*a_rel^(3/2), AN_i*400);
%[t,orbit_final] = ode113(@(t,y) internalIntegration3D(t,y), tspan_fin,y0,options);


figure(6)
scatter(T_all_days(good_condition),J_all(good_condition),'.')
title('Парето-Фронт')
grid()

% figure(7);
% scatter3(rf_all(good_condition,1),rf_all(good_condition,2),180*i_all(good_condition)/pi, '.')
% figure(8);
% scatter3(rf_all(good_condition,1),rf_all(good_condition,2),e_all(good_condition), '.')
% figure(9);
% %scatter3(rf_all(good_condition,1),rf_all(good_condition,2),diff_vf(good_condition), '.')
% scatter3(rf_all(:,1),rf_all(:,2),diff_vf(:), '.')
% figure(10);
% %scatter3(rf_all(good_condition,1),rf_all(good_condition,2),diff_vf(good_condition), '.')
% plot(T_all_days,diff_vf, '.')
%выводим график
figure(11);
plot3(0, 0, 0, 'k--o');
set(gca,'FontSize',14);
set(gcf,'Color','white')
hold on;

plot3(orbit_initial(:, 1), orbit_initial(:, 2), orbit_initial(:, 3), 'k--');
plot3(start_pos(1),start_pos(2), start_pos(3), 'bO', 'LineWidth', 1);

scatter3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3)-1, 'k.');
%scatter3(mars_traj(180:200, 1), mars_traj(180:200, 2), mars_traj(180:200, 3)-1, 'bO');
%plot3(rf(1), rf(2), rf(3),'bO')
%rf_cloud = reshape(rf_massive, [numel(J_massive),3]);

%scatter3(rf_all(180:200,1),rf_all(180:200,2),rf_all(180:200,3)-1, '*', 'MarkerEdgeColor',"red")
scatter3(rf_all(:,1),rf_all(:,2),rf_all(:,3)-1, '.', 'MarkerEdgeColor',"#77AC30")
scatter3(rf_all(ok_condition,1),rf_all(ok_condition,2),rf_all(ok_condition,3), '.b')
scatter3(rf_all(good_condition,1),rf_all(good_condition,2),rf_all(good_condition,3)+2, '*r')
hold off;
axis equal

%title('Траектория КА')
xlabel('x, DU')
ylabel('y, DU')
zlabel('z, DU')
view(0,90)
%view(90,0)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
grid on;
box off;

function R = calculateRotMatrix(r0,v0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
e1 = r0;
e2 = v0;
e3 = cross(e1,e2);
R = [e1;e2;e3];
end
function R_eo = calculateReoMatrix(k, e_target,omega_target)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%настройка эксцентриситета
theta = 0.642;
%коэффициенты для линейной аппроксимации без весов
% C_ecc_0 = (k-1)*3.565068035759543;
% C_ecc_1 = (k-1)*-2.3261252072805805;
% C_ecc_2 = (k-1)*-0.6264065291171595;

%коэффициенты для линейной аппроксимации с весами
C_ecc_0 = (k-1)*3.272584669562773;
C_ecc_1 = (k-1)*-1.9315126322573586;
C_ecc_2 = (k-1)*-0.4543650380418002;

C_ecc = C_ecc_0 + sin(omega_target)*(C_ecc_1*sin(omega_target)+C_ecc_2*sin(3*omega_target));
dm = -e_target/C_ecc;
m = -sin(theta)+dm;%X_rotated
n = cos(theta); %Y_rotated 

e_param = m*cos(theta)+n*sin(theta);
a_param = -m*sin(theta)+n*cos(theta);

R_e = [a_param*(1-e_param), 0, 0;
    0, (1+e_param)/sqrt(a_param*(1-e_param^2)), 0;
    0, 0, (1+e_param)/sqrt(a_param*(1-e_param^2))*a_param*(1-e_param)];

o_input = omega_target/2+(1/6)*sin(2*omega_target)+(1/32)*sin(4*omega_target);

R_rot = [cos(o_input),-sin(o_input),0;
        sin(o_input),cos(o_input),0;
        0,0,1];
R_eo=R_rot*R_e*R_rot';
end
function R_iO = calculateRiOMatrix(k, i_target,Omega_target)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Ci0 = (0.38159616342652647-0.0243395353632856*k-0.2767778010385592*exp(-k^2))*log(k);
%Ci1 = -k*exp(-7.8328443*sinh(k))+5.5964956*exp(-sinh((1.38935843796911*k^2+1.09321695795059)/(0.7103882*k^2+0.04272333)));

Ci0 = (-(k*(0.039389372*k - 0.2516674) - 0.044592466)*log(k));
Ci1 = 0.0012084354*(k-1.4255633)*(k^3-0.99985284);
delta_i_coef = Ci0+Ci1*cos(2*Omega_target);
if i_target == 0
    j=0;
else
    j = -i_target/delta_i_coef;
end
R_i = [
    1, 0, 0;
    0, 1/(1+j^2), -j/(1+j^2);
    0, j/(1+j^2), 1/(1+j^2);
    ];

%восходящий узел
R_rot = [cos(Omega_target),-sin(Omega_target),0;
        sin(Omega_target),cos(Omega_target),0;
        0,0,1];
R_iO=R_rot*R_i*R_rot';
end
function [pr_0,pv_0] = get_initial_adjoint(a,k,d)
%Adjoint variables approximation
% a - angular distance, revolutions
% k - radii ratio
% d - multiplier of D coefficients, 0 or 1
%MARS
pr_x_M = @(a,d)(0.20074211987250948)/(7.780825627926492*a-0.34159993346898454...
    +sin(2*pi*a+0.1231174475945867+d*0.5853586584957071/(a^2-1.5690074915306857*a+1.9819223815154108))...
    + d*(-41.26470185582318+86.35986725318767*a-51.08371071757919*a^2)*exp(-3.424188046763629*a) ); %OK
pr_y_M = @(a,d)(0.0015463527047873246 - 0.0031678016585249486 * cos(2*pi * a+...
    d*(-2.7056930690628103*a^2+5.176300671261158*a-2.336378597152512)*exp(3.2003639491958467*a-3.2205106566099353*a^2)))/(a^2-0.18328217249773662*a+0.11514749689845788);%OK
pv_x_M = @(a,d)(0.0018661818762522732 - 0.0031537940448213170 * cos(2*pi * a+...
    d*(-7.277086750486235*a^2+13.74021831427378*a-6.262552219285274)*exp(1.381195073394106*a-2.278956114228774*a^2)))/(a^2-0.19152227381595716*a+0.0993237695309777);
pv_y_M = @(a,d)(0.10456877009148134)/(4.050922375438426*a-0.15902042687671913...
    +sin(2*pi*a+0.11429217276196703+d*0.012218372217935064/(a^2-2.8347986183365705*a+2.06664561825753)) ...
    +d*(-540.5981059117084*a^2+958.0032106019531*a-415.15730979365355)*exp(-5.877538308252977*a));
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

k1=0.72;
k2=1.52;

%prx
p1 = pr_x_V(a,d);
p2 = pr_x_M(a,d);
% f = (p1/log(k1)-(k2/k1)*p2/log(k2))/(1-k2/k1);
% g = p2*k2/log(k2)-f*k2;
% pr_x = @(k)(log(k)*(f*k+g)/k);
g = (p2*power(k2,1/4)/log(k2) - p1*k1*power(k2,-3/4)/log(k1))/(1-power(k1/k2,3/4));
f = p1*k1/log(k1)-g*power(k1,3/4);
pr_x = @(k)(log(k)*(f/k+g*power(k,-1/4)));
%pry
p1 = pr_y_V(a,d);
p2 = pr_y_M(a,d);
g = (p2*k2^(3/2)/log(k2) - p1*k1*sqrt(k2)/log(k1))/(1-sqrt(k2/k1));
f = p1*k1/log(k1)-g/sqrt(k1);
pr_y = @(k)(log(k)*(f/k+g/(k^(3/2))));
%pvx
p1 = pv_x_V(a,d);
p2 = pv_x_M(a,d);
g = (p2*k2^(3/2)/log(k2) - p1*k1*sqrt(k2)/log(k1))/(1-sqrt(k2/k1));
f = p1*k1/log(k1)-g/sqrt(k1);
pv_x = @(k)(log(k)*(f/k+g/(k^(3/2))));
%pvy
p1 = pv_y_V(a,d);
p2 = pv_y_M(a,d);
% f = (p1/log(k1)-(k2/k1)*p2/log(k2))/(1-k2/k1);
% g = p2*k2/log(k2)-f*k2;
% pv_y = @(k)(log(k)*(f*k+g)/k);
g = (p2*power(k2,1/4)/log(k2) - p1*k1*power(k2,-3/4)/log(k1))/(1-power(k1/k2,3/4));
f = p1*k1/log(k1)-g*power(k1,3/4);
pv_y = @(k)(log(k)*(f/k+g*power(k,-1/4)));

pr_0 = [pr_x(k);pr_y(k);0.];
pv_0 = [pv_x(k);pv_y(k);0.]; 
end
function [rf,vf,Jf,r_all] = integrate_Rtotal_equinoctial(R_total, start_pos,start_vel,p0,t0,t_end)
p0 = reshape(p0,[6,1]);
pr = R_total^(-1)*p0(1:3);
pv = R_total^(-1)*p0(4:6);

mug = 1;

[~,eMag,i,O,~,nu,~,~,lonPer,p] = rv2orb(start_pos',start_vel',mug);

[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
X0 = [p;ex;ey;ix;iy;L];
dxdX = equitoctial2decart_jacobian(X0,mug);

P0 = dxdX'*[pr;pv];
y0 = [X0;P0;0];
tspan = linspace(t0,t0+t_end, 2);
options = odeset('AbsTol',1e-8);
options = odeset(options,'RelTol',1e-8);   
[~,y] = ode113(@(t,y) equitoctial_adjoint_for_integration_withJ(y,mug), tspan,y0,options);

traj = arrayfun(@(p,ex,ey,ix,iy,L) equitoctial2decart([p;ex;ey;ix;iy;L], mug),...
    y(:, 1),y(:, 2),y(:, 3),y(:, 4),y(:, 5),y(:, 6),'UniformOutput',false);
traj = cell2mat(traj')';
rf = traj(end, 1:3)';
vf = traj(end, 4:6)';
Jf = y(end, 13);
r_all = traj(:, 1:3);
%[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(rf,vf,1);
%[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
end
function [rf,vf,Jf,r_all] = integrate_Rtotal(R_total, start_pos,start_vel,p0,t0,t_end)
p0 = reshape(p0,[6,1]);
pr = R_total^(-1)*p0(1:3);
pv = R_total^(-1)*p0(4:6);

mug = 1;


y0 = [start_pos,start_vel,pr',pv',0,0];
tspan = linspace(t0,t0+t_end, 2);
options = odeset('AbsTol',1e-8);
options = odeset(options,'RelTol',1e-8);   
[~,traj] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);

rf = traj(end, 1:3)';
vf = traj(end, 4:6)';
Jf = traj(end, 13);
r_all = traj(:, 1:3);
%[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(rf,vf,1);
%[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
end