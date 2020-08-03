% Обходим витки и ищем минимум функционала
clc;
clear;
symbolic_Jacob
warning('off');

angle = 6*pi/6;
rad=0.01;
case_traj=2;
[functional, dis, s, y] = trajectorySearch(0, angle, rad, case_traj, symF);
n = 0;
for i = 1:10
    [functional_tmp, dis_tmp, s_tmp, y_tmp] = trajectorySearch(i, angle, rad, case_traj, symF);
    i
    if (functional_tmp < functional) && (dis_tmp <= dis*100)
        functional = functional_tmp
        dis = dis_tmp;
        s = s_tmp;
        y = y_tmp;
        n = i
    end
end

u = y(:, 1:4);
r=zeros(length(u),4);
for i = 1:length(u)
    rr = u(i,:);
    L = [[rr(1) -rr(2) -rr(3) rr(4)];
    [rr(2) rr(1) -rr(4) -rr(3)];
    [rr(3) rr(4) rr(1) rr(2)];
    [rr(4) -rr(3) rr(2) -rr(1)]];
    r(i,:)=L*rr';
end

ae = 149597870700;
mug = 132712.43994*(10^6)*(10^(3*3));

T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;

modifier=1e-8;
tf_a = T_earth*(n + angle-rad)*modifier;
tf_b = T_earth*(n + angle+rad)*modifier;

x_tmp=T_earth*(n + angle)*modifier;

n_M = floor((x_tmp/modifier)/T_mars);
angle_M = (x_tmp/modifier)/T_mars-n_M;
t_Mars_0 = (angle-angle_M-0.03)*T_mars;

u_f=y(end, 1:4);
v_f=y(end, 5:8);
h_f=y(end, 9);
tau_f=y(end, 10);
tf = tau_f-2*(u_f'*v_f)/(-2*h_f);

n_M = floor((tf+t_Mars_0)/T_mars);
angle_M = ((tf+t_Mars_0)/T_mars-n_M)*2*pi;

%Проверка "на глаз"
plot(0, 0,'y--o')
hold on;
th = 0:pi/50:2*pi;
plot(ae*cos(th),ae*sin(th),'k');
plot(1.52*ae*cos(th),1.52*ae*sin(th),'r');
plot(r(:, 1), r(:, 2),'b')
plot(r(end, 1), r(end, 2),'bO')
plot(1.52*ae*cos(angle_M), 1.52*ae*sin(angle_M),'rO')
axis equal
hold off;