%Эот скрипт перебирает угловые дальности с заданным радиусом поиска
t_start = juliandate(2022,1,1);
display=1;
UorR='u_hat';
N=1350;
m0=367;
eta=0.45;
case_traj = 2;
step = 1/32;
ds = 3/4:step:12/2;
rad = step/2;
L=length(ds);
T_NONLINEAR=zeros([1e+03,L]);
RR=zeros([1e+03,3,L]);
T=zeros([1,L]);
T_cont=zeros([1,L]);
T_END=zeros([1,L]);
DR=zeros([1,L]);
DV=zeros([1,L]);
CONV=zeros([1,L]);
CONV_CONT=zeros([1,L]);
SF=zeros([1,L]);
PHI=zeros([1,L]);
%Затраты массы
M=zeros([1,L]);
M_cont=zeros([1,L]);
%Разница по координатам вдоль траектории
D=zeros([1,L]);
PX=zeros([8,L]);
S=zeros([1,L]);
modifier_p=1e-04;
modifier_f=1e+10;
%Начальное приближение
x0=zeros([1, 10]);
warning('off');
planet_start = 'Earth';
planet_end = 'Mars';
decreaseNonPsysical = 0;
skipPoint=zeros([1,L]);

px0=zeros([1, 8]);
phi0=0;

rad=1/8;

for i=1:L
    i
    ds(i)
    %x0=zeros([1, 10]);
    
    x0=[px0 0 phi0/(2*pi)];
    x0(9)=ds(i);
    delta_s=ds(i);

    integration_acc=1e-10;
    %Одиночный запуск метода и получение всех необходимых для графиков
    %переменных
    n = floor(delta_s);
    modifier_p=10^(-4-sqrt(delta_s));
    display = 0;
    terminal_state = 's';
    UorR = 'u_hat';
    
    %delta_s=1.23*(n+angle)-0.24;
    %delta_s=1.2*(n+angle)-0.2;
    %delta_s=n+angle;
    evaluation_time=0;
    if rad > 0
        calculate_condition=0;
        [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] = checkMethod(t_start,delta_s,rad,UorR,decreaseNonPsysical,modifier_p,modifier_f,x0,eta, case_traj,planet_end, display,terminal_state,integration_acc,calculate_condition);
        % 
        if terminal_state == 's'
            x0_sec = [px s_f/(2*pi) phi/(2*pi)];
        elseif terminal_state == 't'
            x0_sec = [px t_end/365.256363004 phi/(2*pi)];
        end
    end
    rad=0;
    % decreaseNonPsysical=0;
    calculate_condition=1;
    [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time_2] = checkMethod(t_start,delta_s,rad,UorR,decreaseNonPsysical,modifier_p,modifier_f,x0_sec,eta, case_traj,planet_end, display,terminal_state,integration_acc,calculate_condition);
    evaluation_time=evaluation_time+evaluation_time_2;
    %Неудачное интегрирование
    if length(t) < 1000
        skipPoint(i)=1;
    end
    %Метод не сошёлся
    if dr > 10000
        skipPoint(i)=2;
    end
    
    if length(t) < 1000 ||  dr > 10000
        px0=zeros([1, 8]);
        phi0=0;
        rad=1/8;
        %сразу же пересчитываем неудачную точку с раширенной областью
        x0=[px0 0 phi0/(2*pi)];
        x0(9)=ds(i);

        calculate_condition=0;
        [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] = checkMethod(t_start,delta_s,rad,UorR,decreaseNonPsysical,modifier_p,modifier_f,x0,eta, case_traj,planet_end, display,terminal_state,integration_acc,calculate_condition);
        x0_sec = [px s_f/(2*pi) phi/(2*pi)];

        rad=0;
        calculate_condition=1;
        [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time_2] = checkMethod(t_start,delta_s,rad,UorR,decreaseNonPsysical,modifier_p,modifier_f,x0_sec,eta, case_traj,planet_end, display,terminal_state,integration_acc,calculate_condition);
        evaluation_time=evaluation_time+evaluation_time_2;
        
        if length(t) == 1000 && dr < 10000
            skipPoint(i)=0;
        else
            continue
        end

        
    end
    px0=px;
    phi0=phi;


    T_NONLINEAR(:,i) = t;
    T(i)=evaluation_time;
    %S(i)=s(end);
    S(i)=delta_s;
    DR(i)=dr;
    DV(i)=dV;
    CONV(i)=C;
    PX(:,i)=px;
    SF(i)=s_f;
    PHI(i)=phi;
    m=massLP(Jt, m0, N);
    T_END(i)=t_end;
    M(i)=m(1)-m(end);
    RR(:,:,i)=rr(:,1:3);
end

px_new=zeros([1, 8]);
phi_new=0;

for i=L:-1:1
    if skipPoint(i) < 1
        px_new=PX(:,i)';
        phi_new=PHI(i);
        continue
    end
    i
    ds(i)
    delta_s=ds(i);
    modifier_p=10^(-4-sqrt(delta_s));
    integration_acc=1e-12;
    %Одиночный запуск метода и получение всех необходимых для графиков
    %переменных

    display = 0;
    terminal_state = 's';
    UorR = 'u_hat';
    
    x0_sec = [px_new delta_s/(2*pi) phi_new/(2*pi)];
    rad=0;
    calculate_condition=1;
    [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] = checkMethod(t_start,delta_s,rad,UorR,decreaseNonPsysical,modifier_p,modifier_f,x0_sec,eta, case_traj,planet_end, display,terminal_state,integration_acc,calculate_condition);
    %evaluation_time=evaluation_time+evaluation_time_2;

    if dr < 10000
        skipPoint(i)=-1;
        px_new=px;
        phi_new=phi;
        s_f_new=s_f;
        t_end_new=t_end;
    else
        continue
    end
    T_NONLINEAR(:,i) = t;
    T(i)=evaluation_time;
    %S(i)=s(end);
    S(i)=delta_s;
    DR(i)=dr;
    DV(i)=dV;
    CONV(i)=C;
    PX(:,i)=px;
    SF(i)=s_f;
    PHI(i)=phi;
    m=massLP(Jt, m0, N);
    T_END(i)=t_end;
    M(i)=m(1)-m(end);
    RR(:,:,i)=rr(:,1:3);
end

%нанесём получивщиеся траектории


%t_end=t(end);
eul = [0 pi/4 0];
rotmZYX = eul2rotm(eul);
%Проверка "на глаз"
figure(1);
plot3(0, 0, 0, 'y--o')
set(gca,'FontSize',14)
hold on;

% mars_traj = 1.52*[cos(th), sin(th), zeros(100,1)];
% earth_traj  = [cos(th), sin(th), zeros(100,1)];
t0 = t_start;
ae = 149597870700;
T_earth = 365.256363004*3600*24;
T_mars_days = 365.256363004*1.8808476;

mug_0 = 132712.43994*(10^6)*(10^(3*3));
V_unit=sqrt(mug_0/ae);
T_unit = T_earth/(2*pi);
ae = 149597870700;
T_mars=T_earth*1.8808476;
T_earth = 365.256363004*3600*24;
t_orbit = linspace(t0,t0+T_earth/(24*3600), 1000);
earth_traj = planetEphemeris(t_orbit','SolarSystem','Earth','430');
earth_traj=earth_traj*1e+03/ae;

earth_traj_New = arrayfun(@(x,y,z)rotmZYX*[x, y, z]', earth_traj(:, 1),earth_traj(:, 2),earth_traj(:, 3),'UniformOutput',false);
earth_traj_New = cell2mat(earth_traj_New')';

t_orbit = linspace(t0,t0+T_mars/(24*3600), 1000);
mars_traj = planetEphemeris(t_orbit','SolarSystem','Mars','430');
mars_traj=mars_traj*1e+03/ae;

mars_traj_New = arrayfun(@(x,y,z)rotmZYX*[x, y, z]', mars_traj(:, 1),mars_traj(:, 2),mars_traj(:, 3),'UniformOutput',false);
mars_traj_New = cell2mat(mars_traj_New')';

plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k')
plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'r')

axis equal

%title('Траектория КА')
xlabel('x, a.e.')
ylabel('y, a.e.')
zlabel('z, a.e.')

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
view(0,90);
box off;
hold off;

% %Выводим траекторию в параметрических переменных"
% figure(5);
% plot3(0, 0, 0, 'y--o')
% set(gca,'FontSize',14)
% hold on;
% 
% th = linspace(0 ,4*pi,1000)';
% 
% mars_traj_ks = arrayfun(@(r1, r2, r3) rToU([r1,r2,r3], 0), mars_traj_New(:, 1),mars_traj_New(:, 2),mars_traj_New(:, 3),'UniformOutput',false);
% mars_traj_ks = cell2mat(mars_traj_ks')';
% earth_traj_ks = arrayfun(@(r1, r2, r3) rToU([r1,r2,r3], 0), earth_traj_New(:, 1),earth_traj_New(:, 2),earth_traj_New(:, 3),'UniformOutput',false);
% earth_traj_ks = cell2mat(earth_traj_ks')';
% plot3(earth_traj_ks(:, 1), earth_traj_ks(:, 2), earth_traj_ks(:, 3), 'k')
% plot3(mars_traj_ks(:, 1), mars_traj_ks(:, 2), mars_traj_ks(:, 3), 'r--')
% plot3(-mars_traj_ks(:, 1), -mars_traj_ks(:, 2), -mars_traj_ks(:, 3), 'r--')
% 
% axis equal
% 
% title('Траектория КА KS')
% xlabel('u1')
% ylabel('u2')
% zlabel('u3')
% 
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% box off;
% hold off;

[r0, V0] = planetEphemeris(t_start,'SolarSystem',planet_start,'430');

eul = [0 pi/4 0];
rotmZYX = eul2rotm(eul);

disp('--------------------')
for i=1:L
    if skipPoint(i)>0
        continue
    end
    rr=RR(:,:,i);
    t = T_NONLINEAR(:,i);
    t_end=T_END(i);
    [rr_cont, Jt_cont, C_cont, evaluation_time_cont, dr_cont, dV_cont, pr0, pv0] =...
        checkContinuation(t_start, t_end, t, case_traj,planet_end,eta, floor(ds(i)));
    
    
    functional_cont = Jt_cont(end);
    m_cont=massLP(Jt_cont, m0, N);
    M_cont(i)=m_cont(1)-m_cont(end);
    %Неправильное число витков влияет на затраты массы
    if M_cont(i)-M(i) > 10
        [rr_cont, Jt_cont, C_cont, evaluation_time_cont, dr_cont, dV_cont, pr0, pv0] =...
            checkContinuation(t_start, t_end, t, case_traj,planet_end,eta, floor(ds(i))+1);
        functional_cont = Jt_cont(end);
        m_cont=massLP(Jt_cont, m0, N);
        M_cont(i)=m_cont(1)-m_cont(end);
    end
    CONV_CONT(i)=C_cont;
    T_cont(i) = evaluation_time_cont;
    figure(1);
    hold on;
    [mars_r_f, mars_v_f]=planetEphemeris([t_start, t_end],'SolarSystem',planet_end,'430');
    mars_r_f=mars_r_f'*1e+03;
    mars_v_f=mars_v_f'*1e+03;
    
    rr_old = arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', rr(:, 1),rr(:, 2),rr(:, 3),'UniformOutput',false);
    rr_old = cell2mat(rr_old')';
    plot3(rr_old(:, 1), rr_old(:, 2), rr_old(:, 3), 'r', 'LineWidth', 2.5);

    plot3(rr_old(end, 1), rr_old(end, 2), rr_old(end, 3),'bO')
    plot3(mars_r_f(1)/ae, mars_r_f(2)/ae,mars_r_f(3)/ae,'rO')
    
    hold off;
    
%     figure(5)
%     hold on;
%     plot3(uu(:, 1), uu(:, 2), uu(:, 3), 'b', 'LineWidth', 2.5);
%     plot3(uu(end, 1), uu(end, 2), uu(end, 3), 'bO');
%     hold off;
    d = cmp2Trajectories(rr_old(:, 1:3)*ae,rr_cont)/ae;
    
    D(i)=d;
    disp(['Величина мнимого времени ', num2str(ds(i)*2*pi,'%10.2f\n'), ' рад.'])
    disp(['Невязка по координате предложенного метода ', num2str(DR(i),'%10.2e\n'), 'м.'])
    disp(['Невязка по скорости предложенного метода ', num2str(DV(i),'%10.2e\n'), 'м/с'])
    disp(['Время работы предложенного метода ', num2str(T(i),'%10.2f\n'), 'сек.'])
    disp(['Время работы метода продолжения по параметру ', num2str(T_cont(i),'%10.2f\n'), 'сек.'])
    disp(['Затраты массы предложенного метода ', num2str(M(i),'%10.2f\n'), 'кг.'])
    disp(['Затраты массы метода продолжения по параметру ', num2str(M_cont(i),'%10.2f\n'), 'кг.'])
    disp(['Максимальная разница в координатах ', num2str(D(i),'%10.2e\n'), 'а.е.'])
    disp(['Число обусловленности в KS-переменных ', num2str(CONV(i),'%10.2e\n')])
    disp(['Число обусловленности в методе продолжения ', num2str(CONV_CONT(i),'%10.2e\n')])
    disp('--------------------')
end


THETA=zeros([1,L]);
for i=1:L
    rr=RR(:,:,i);
    rr_old = arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', rr(:, 1),rr(:, 2),rr(:, 3),'UniformOutput',false);
    rr_old = cell2mat(rr_old')';
    THETA(i)=traj2TrueAnomaly(rr_old);
end

figure(4);
plot(S,THETA,'.')
xlabel("Фиктивное время, 2pi радиан")
ylabel("Угловая дальность, 2pi радиан")
xlim([0 6.5])
ylim([0 6.5])

figure(6);
plot(S,M,'O')
hold on;
plot(S,M_cont,'+')
plot([0, 6.5],[m0, m0],'k')

xlim([0 6.5])

xlabel("Фиктивное время, 2pi радиан")
%ylabel("Угловая дальность, 2pi радиан")
ylabel("Расход топлива, dm, кг")
set(gca,'FontSize',14)
text(3,380,'Начальная масса КА','FontSize',14)
hold off;

figure(10);
plot(T_END/365,M,'O')
hold on;
plot(T_END/365,M_cont,'+')
plot([0, 9],[m0, m0],'k')

xlim([0 9])

xlabel("Физическое время, годы")
ylabel("Расход топлива, dm, кг")
set(gca,'FontSize',14)
text(3,380,'Начальная масса КА','FontSize',14)
hold off;

%Получаем точки парето-фронта
pareto=[1];
M(skipPoint>0)=m0;
for i=2:L
    if skipPoint(i)>0
        continue
    end

    if sum(M(1:i-1) < M(i)) == 0
        pareto(end+1)=i;
    end
end
%Рисуем парето-фронт
figure(7);
plot(S(pareto),M(pareto),'O')
hold on;
plot(S(pareto),M_cont(pareto),'+')
plot([0, 6.5],[m0, m0],'k')

xlim([0 6.5])

xlabel("Фиктивное время, число витков")
ylabel("Расход топлива, dm, кг")
set(gca,'FontSize',14)
text(3,380,'Начальная масса КА','FontSize',14)
hold off;

%рисуем числа обусловленности
figure(8);
semilogy(S(pareto),CONV(pareto),'O')
hold on;
semilogy(S(pareto),CONV_CONT(pareto),'+')

xlim([0 6.5])

xlabel("Фиктивное время, число витков")
ylabel("Числа обусловленности")
set(gca,'FontSize',14)
hold off;

%рисуем время работы
figure(9);
plot(S(pareto),T(pareto),'O')
hold on;
plot(S(pareto),T_cont(pareto),'+')

xlim([0 6.5])

xlabel("Фиктивное время, число витков")
ylabel("Время работы с")
set(gca,'FontSize',14)
hold off;

% figure(8);
% pareto=T<100;
% pareto_cont=T_cont<100;
% plot(S(pareto),T(pareto),'O')
% hold on;
% plot(S(pareto_cont),T_cont(pareto_cont),'+')
% 
% xlim([0 5.5])
% 
% xlabel("Фиктивное время, число витков")
% ylabel("Время вычислений, секунды")
% set(gca,'FontSize',14)
% hold off;