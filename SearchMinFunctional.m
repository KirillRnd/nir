% Обходим витки и ищем минимум функционала
clc;
clear;
symbolic_Jacob
warning('off');

angle = 6*pi/6;
rad=0.01;
case_traj=1;
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
a=zeros(length(u),4);
t=zeros(length(u),1);
for i = 1:length(u)
    rr = u(i,:)';
    L = [[rr(1) -rr(2) -rr(3) rr(4)];
    [rr(2) rr(1) -rr(4) -rr(3)];
    [rr(3) rr(4) rr(1) rr(2)];
    [rr(4) -rr(3) rr(2) -rr(1)]];
    r(i,:)=L*rr;
    u2=norm(u)^2;
    v=y(i, 5:8)';
    h=y(i, 9)';
    tau=y(i ,10)';
    pu=y(i, 11:14)';
    pv=y(i, 15:18)';
    ph=y(i, 19)';
    ptau=y(i, 20)';
    aa=L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*(rr'*rr)*rr/(-2*h)^(3/2));
    
    La = [[aa(1) -aa(2) -aa(3) aa(4)];
    [aa(2) aa(1) -aa(4) -aa(3)];
    [aa(3) aa(4) aa(1) aa(2)];
    [aa(4) -aa(3) aa(2) -aa(1)]];
    a(i, :)=La*aa;
    t(i) = tau-2*(rr'*v)/(-2*h);
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

figure(2);
plot(t, vecnorm(a, 2, 2))
figure(3);
plot(s, t)
%Проверка "на глаз"
figure(1);
plot(0, 0,'y--o')
hold on;
th = 0:pi/50:2*pi;
plot(ae*cos(th),ae*sin(th),'k');
plot(1.52*ae*cos(th),1.52*ae*sin(th),'r');
plot(r(:, 1), r(:, 2),'b')
a_scale=1e+10;

d = 24*3600;
idxes=1;
for i=1:ceil(t(end)/d)
    ix = find(t>d*i*10, 1);
    idxes=[idxes, ix];
end    
for i = idxes
    plot([r(i, 1), r(i, 1)+a_scale*a(i, 1)], [r(i, 2), r(i, 2)+a_scale*a(i, 2)],'k')
end
plot(r(end, 1), r(end, 2),'bO')
plot(1.52*ae*cos(angle_M), 1.52*ae*sin(angle_M),'rO')
axis equal
hold off;
