%clearvars -except symF
%clc;
%5if exist('symF','var') ~= 1
%    symbolic_Jacob
%end
%t_start = juliandate(2022,1,1);
t_start=0;
orbits='Flat';
N=1350;
m0=367;
eta=0.45;
%условия на fmincon
%ЗАДАЧА ПРОЛЁТА case_traj=1; ЗАДАЧА сопровождения case_traj=2;
case_traj=2;
%Выбор сходимости по физическим координатам ('r') или по параметрическим ('u')
decreaseNonPsysical = 0;
%Начальные условия
%x0(1:8)=PX(:,41)';
%x0(1:8)=[-0.081312208804168  -0.165241434409968  -0.137910347360105  -0.029456685044885  -0.151663938192582  -0.306939705835856  -0.241150551359049  -0.080582966706315];
%x0=[0.0103    0.0048    0.0046    0.0048   -0.0036    0.0065    0.0162    0.0104    1.7770         0];
A = [];
b = [];
Aeq = [];
beq = [];
ae = 149597870700;
mug_0 = 132712.43994*(10^6)*(10^(3*3));
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;
T_mars_days = 365.256363004*1.8808476;

eul = [pi/12 pi/4 pi/12];
rotmZYX = eul2rotm(eul);

r_unit=ae;
V_unit=sqrt(mug_0/ae);
T_unit = T_earth/(2*pi);
planet_start = 'Earth';
planet_end = 'Mars';

mug=1;

n=0;
angle=0.7496;
delta_s=(n+angle)*2*pi;
j=169;
%delta_s=ds(j)*2*pi;
x0=zeros([1, 10]);
%x0(12)=x0(11)/2;
%modifier_p=1e-08;
modifier_p=10^(-4-sqrt(delta_s));
modifier_f=1e+10;
integration_acc=1e-10;

%Одиночный запуск метода и получение всех необходимых для графиков
%переменных
display = 1;
terminal_state = 's';
UorR = 'u_hat';
rad=0;
omega = 1.2257;
%x0(9)=n+angle+rad/2;
x0(1:8)=PXevery_a(1,113,:);
x0(9)=delta_s/(2*pi);
calculate_condition=1;

checkMethod_params.t_start = t_start;
checkMethod_params.delta_s = delta_s;
checkMethod_params.rad = rad;
checkMethod_params.UorR = 'u_hat';
checkMethod_params.decreaseNonPsysical = decreaseNonPsysical;
checkMethod_params.modifier_p = modifier_p;
checkMethod_params.modifier_f = modifier_f;
checkMethod_params.x0_sec = x0;
checkMethod_params.eta = eta;
checkMethod_params.case_traj = case_traj;
checkMethod_params.planet_start = planet_start;
checkMethod_params.planet_end = planet_end;
checkMethod_params.display = display;
checkMethod_params.terminal_state = terminal_state;
checkMethod_params.integration_acc = integration_acc;
checkMethod_params.calculate_condition = calculate_condition;
checkMethod_params.orbits = orbits;
checkMethod_params.omega = omega;
checkMethod_params.a_rel = 1.52;

[dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] = checkMethod(checkMethod_params);
% 
if terminal_state == 's'
    x0_sec = [px s_f/(2*pi) phi/(2*pi)];
elseif terminal_state == 't'
    x0_sec = [px t_end/365.256363004 phi/(2*pi)];
end
%modifier_p=1e-08;
% terminal_state = 's';
% UorR = 'u_hat';
% integration_acc=1e-16;
rad=0;
% decreaseNonPsysical=0;
calculate_condition=0;
% [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time_2] = checkMethod(t_start,n+angle,rad,UorR,decreaseNonPsysical,modifier_p,modifier_f,x0_sec,eta, case_traj,planet_end, display,terminal_state,integration_acc,calculate_condition, orbits, omega);
% evaluation_time=evaluation_time+evaluation_time_2;

%Коррекция фи для графика
phi=atan2(uu(end,4),uu(end,1));
phi0=atan2(uu(1,4),uu(1,1));
functional = Jt(end);
z0=zeros([1, 6]);
[rr_cont, Jt_cont, C_cont, evaluation_time_cont, dr_cont, dV_cont, pr0, pv0] =...
    checkContinuation(t_start, t_end, t,z0, case_traj,planet_end,eta, n, orbits, omega);
functional_cont = Jt_cont(end);
m_cont=massLP(Jt_cont, m0, N);
%t_end=t(end);
figure(2);
plot(t/T_earth, vecnorm(a_ks, 2, 2)*1e+03, 'LineWidth', 3);
%title('Зависимость ускорения силы тяги от времени')
xlabel('Физическое время, годы','FontSize',14)
ylabel('Реактивное ускорение, мм/c^2','FontSize',14)
box off;
set(gca,'FontSize',14)
set(gca,'TickLabelInterpreter', 'latex')
grid on;

figure(3);
plot(t/T_earth, s/(2*pi), 'LineWidth', 3);
title('Зависимость мнимого времени от физического')
xlabel('Физическое время, годы')
ylabel('Фиктивное время, 2\pi радиан')
box off;
set(gca,'FontSize',14)
set(gca,'TickLabelInterpreter', 'latex')
grid on;

figure(4);
m=massLP(Jt, m0, N);
plot(t/T_earth, m, 'LineWidth', 3);
title('Зависимость массы от времени')
xlabel('Физическое время, годы')
ylabel('m, масса, кг')
box off;
set(gca,'FontSize',14)
set(gca,'TickLabelInterpreter', 'latex')
grid on;

%Проверка "на глаз"
figure(1);
plot3(0, 0, 0, 'k--o')
set(gca,'FontSize',14)
hold on;

t0 = t_start;
t_orbit = linspace(t0,t0+T_earth/(24*3600), 1000);

st.t = t_orbit';
st.planet = 'Earth';
st.mode = orbits;
st.delta_omega = omega;

earth_traj = planetModel(st);
earth_traj=earth_traj*1e+03/ae;
%Для KS
earth_traj_New = arrayfun(@(x,y,z)rotmZYX*[x, y, z]', earth_traj(:, 1),earth_traj(:, 2),earth_traj(:, 3),'UniformOutput',false);
earth_traj_New = cell2mat(earth_traj_New')';

t_orbit = linspace(t0,t0+T_mars/(24*3600), 1000);

st.t = t_orbit';
st.planet = planet_end;
st.mode = orbits;
st.delta_omega = omega;

mars_traj = planetModel(st);
mars_traj=mars_traj*1e+03/ae;

%Для KS
mars_traj_New = arrayfun(@(x,y,z)rotmZYX*[x, y, z]', mars_traj(:, 1),mars_traj(:, 2),mars_traj(:, 3),'UniformOutput',false);
mars_traj_New = cell2mat(mars_traj_New')';

plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k')
plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'k')

st.t = [t_start, t_end];
st.planet = planet_end;
st.mode = orbits;
st.delta_omega = omega;

[mars_r_f, mars_v_f]=planetModel(st);
mars_r_f=mars_r_f'*1e+03;
mars_v_f=mars_v_f'*1e+03;

rr_old = arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', rr(:, 1),rr(:, 2),rr(:, 3),'UniformOutput',false);
rr_old = cell2mat(rr_old')';

VV_old = arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', VV(:, 1),VV(:, 2),VV(:, 3),'UniformOutput',false);
VV_old = cell2mat(VV_old')';


plot3(rr_old(:, 1), rr_old(:, 2), rr_old(:, 3), 'cyan', 'LineWidth', 1);


a_ks_old= arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', a_ks(:, 1),a_ks(:, 2),a_ks(:, 3),'UniformOutput',false);
a_ks_old = cell2mat(a_ks_old')';

a_scale=3e-01/mean(vecnorm(a_ks, 2, 2));
%a_scale=0;
d = 24*3600;
idxes=1;
for i=1:ceil(t(end)/d)
    ix = find(t>d*i*10, 1);
    idxes=[idxes, ix];
end   
%тяга
for i = idxes
    plot3([rr_old(i, 1), rr_old(i, 1)+a_scale*a_ks_old(i, 1)], [rr_old(i, 2), rr_old(i, 2)+...
        a_scale*a_ks_old(i, 2)],[rr_old(i, 3), rr_old(i, 3)+a_scale*a_ks_old(i, 3)],'k')
end

plot3(rr_old(end, 1), rr_old(end, 2), rr_old(end, 3),'bO')
plot3(mars_r_f(1)/ae, mars_r_f(2)/ae,mars_r_f(3)/ae,'rO')

plot3(rr_cont(:, 1)/ae, rr_cont(:, 2)/ae, rr_cont(:, 3)/ae, 'g', 'LineWidth', 2.5);
%эти две точки должны находиться рядом
%plot3(rr_cont(500, 1)/ae, rr_cont(500, 2)/ae, rr_cont(500, 3)/ae, 'gO', 'LineWidth', 2.5);
%plot3(rr_old(500, 1), rr_old(500, 2), rr_old(500, 3), 'bO', 'LineWidth', 2.5);
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
hold off;

%Выводим траекторию в параметрических переменных"
figure(5);
plot3(0, 0, 0, 'k--o')
set(gca,'FontSize',14)
hold on;

th = linspace(0 ,4*pi,1000)';

mars_traj_ks = arrayfun(@(r1, r2, r3) rToU([r1,r2,r3], phi), mars_traj_New(:, 1),mars_traj_New(:, 2),mars_traj_New(:, 3),'UniformOutput',false);
mars_traj_ks = cell2mat(mars_traj_ks')';

mars_traj_ks_zero = arrayfun(@(r1, r2, r3) rToU([r1,r2,r3], phi0), mars_traj_New(:, 1),mars_traj_New(:, 2),mars_traj_New(:, 3),'UniformOutput',false);
mars_traj_ks_zero = cell2mat(mars_traj_ks_zero')';

%phi === 0
earth_traj_ks = arrayfun(@(r1, r2, r3) rToU([r1,r2,r3], phi0), earth_traj_New(:, 1),earth_traj_New(:, 2),earth_traj_New(:, 3),'UniformOutput',false);
earth_traj_ks = cell2mat(earth_traj_ks')';
plot3(earth_traj_ks(:, 1), earth_traj_ks(:, 2), earth_traj_ks(:, 3), 'k')
plot3(mars_traj_ks(:, 1), mars_traj_ks(:, 2), mars_traj_ks(:, 3), 'r')
%plot3(-mars_traj_ks_zero(:, 1), -mars_traj_ks_zero(:, 2), -mars_traj_ks_zero(:, 3), 'r--')
%plot3(mars_traj_ks_zero(:, 1), mars_traj_ks_zero(:, 2), mars_traj_ks_zero(:, 3), 'r--')
plot3(uu(:, 1), uu(:, 2), uu(:, 3), 'b', 'LineWidth', 2.5);
%a_scale=3e-01/mean(vecnorm(a_ks, 2, 2));
a_scale=0;
d = 24*3600;
idxes=1;
for i=1:ceil(t(end)/d)
    ix = find(t>d*i*10, 1);
    idxes=[idxes, ix];
end    
for i = idxes
    %plot3([uu(i, 1), uu(i, 1)+a_scale*a_ks(i, 1)], [uu(i, 2), uu(i, 2)+a_scale*a_ks(i, 2)],[uu(i, 3), uu(i, 3)+a_scale*a_ks(i, 3)],'k')
end

%plot3(uu(end, 1), uu(end, 2), uu(end, 3),'bO')
%plot3(mars_r_f(1), mars_r_f(2),mars_r_f(3),'rO')
axis equal

title('Траектория КА KS')
xlabel('u1')
ylabel('u2')
zlabel('u3')
grid on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
hold off;
d = cmp2Trajectories(rr_old(:, 1:3)*ae,rr_cont);

disp(['Максимальная разница в  координатах ', num2str(d/ae,'%10.2e\n'), 'а.е.'])
disp(['Расход массы в KS-координатах ', num2str(m(1)-m(end),'%10.5e\n'), 'кг'])
disp(['Расход массы методом продолжения ', num2str(m_cont(1)-m_cont(end),'%10.6e\n'), 'кг'])
disp(['Невязка координаты ', num2str(norm(ae*rr_old(end, 1:3)-mars_r_f(1:3)'),'%10.2e\n'),',м'])
disp(['Невязка скорости ', num2str((norm(V_unit*VV_old(end, 1:3)-mars_v_f(1:3)')),'%10.2e\n'),',м/с'])
% относительное число обусловленности
disp(['Число обусловленности в KS-переменных ', num2str(C,'%10.2e\n')])
disp(['Число обусловленности в методе продолжения ', num2str(C_cont,'%10.2e\n')])
% абсолютное число обусловленности
%disp(['Абсолютное число обусловленности ', num2str(1/norm(grad),'%10.2e\n')])