% Обходим витки и ищем минимум функционала
clearvars -except symF
clc;
if exist('symF','var') ~= 1
    symbolic_Jacob
end
warning('off');

N=1350;
m0=367;
eta=0.45;

n_max=4;

angle = 0.5;
rad=0.1;
case_traj=2;
n = 0
[Jt, dis, s, y] = trajectorySearch(n, angle, rad, case_traj, symF, eta);
functional = Jt(end)

dm=zeros(n_max+1, 1);
m=massLP(Jt, m0, N);
m(1)-m(end)
dm(1)=m(1)-m(end);
dis
for i = 1:n_max
    [Jt_tmp, dis_tmp, s_tmp, y_tmp] = trajectorySearch(i, angle, rad, case_traj, symF, eta);
    i
    functional_tmp = Jt_tmp(end)
    m=massLP(Jt_tmp, m0, N);
    m(1)-m(end)
    dm(i+1)=m(1)-m(end);
    dis_tmp
    if functional_tmp < functional
        functional = functional_tmp;
        Jt = Jt_tmp;
        dis = dis_tmp;
        s = s_tmp;
        y = y_tmp;
        n = i;
    end
end

uu = y(:, 1:4);
rr=zeros(length(uu),4);
a=zeros(length(uu),4);
t=zeros(length(uu),1);
VV=zeros(length(uu),4);
for i = 1:length(uu)
    u = uu(i,:)';
    r=KS(u);
    rr(i,:)=r;
    L=L_KS(u);
    u2=norm(u)^2;
    v=y(i, 5:8)';
    h=y(i, 9)';
    tau=y(i ,10)';
    pu=y(i, 11:14)';
    pv=y(i, 15:18)';
    ph=y(i, 19)';
    ptau=y(i, 20)';
    %aa=L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*(rr'*rr)*rr/(-2*h)^(3/2));
    res=symF(h,ph,ptau,pu(1),pu(2),pu(3),pu(4),pv(1),pv(2),pv(3),pv(4),u(1),u(2),u(3),u(4),v(1),v(2),v(3),v(4));
    dvds=res(5:8);
    dhds=res(9);
    V = 2*sqrt(-2*h)*L*v/(u2);
    VV(i, :)=V;
    a(i, :)=(-2*h/(norm(r)^2))*(2*(L_KS(v)*v+L_KS(u)*dvds)-(2*u'*v/(sqrt(-2*h)) + norm(r)*dhds/((-2*h)^(3/2)))*V)+mug*r/(norm(r)^3);
    %a(i, :)=KS(aa);
    t(i) = tau-2*(u'*v)/(-2*h);
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
t_Mars_0 = (angle-angle_M)*T_mars;

u_f=y(end, 1:4);
v_f=y(end, 5:8);
h_f=y(end, 9);
tau_f=y(end, 10);
tf = tau_f-2*(u_f'*v_f)/(-2*h_f);

n_M = floor((tf+t_Mars_0)/T_mars);
angle_M = ((tf+t_Mars_0)/T_mars-n_M)*2*pi;

figure(2);
plot(t/(24*3600), vecnorm(a, 2, 2)*1e+03);
title('Зависимость ускорения силы тяги от времени')
xlabel('t, время в днях')
ylabel('a, ускорение силы тяги. мм/с^2')

figure(3);
plot(t/(24*3600), s);
title('Зависимость мнимого времени от обычного')
xlabel('t, время в днях')
ylabel('s, мнимое время')
figure(4);
m=massLP(Jt, m0, N);
plot(t/(24*3600), m);
title('Зависимость массы от времени')
xlabel('t, время в днях')
ylabel('m, масса, кг')

figure(5);
plot(0:n_max, dm, '*');
title('Зависимость расхода топлива от количества витков')
xlabel('n, витки')
ylabel('dm, масса, кг')

%Проверка "на глаз"
figure(1);
plot(0, 0,'y--o')
hold on;
th = 0:pi/50:2*pi;
plot(ae*cos(th),ae*sin(th),'k');
plot(1.52*ae*cos(th),1.52*ae*sin(th),'r');

a_scale=3e+10/mean(vecnorm(a, 2, 2));

d = 24*3600;
idxes=1;
for i=1:ceil(t(end)/d)
    ix = find(t>d*i*10, 1);
    idxes=[idxes, ix];
end    
for i = idxes
    plot([rr(i, 1), rr(i, 1)+a_scale*a(i, 1)], [rr(i, 2), rr(i, 2)+a_scale*a(i, 2)],'k')
end
plot(rr(:, 1), rr(:, 2),'b')
plot(rr(end, 1), rr(end, 2),'bO')
plot(1.52*ae*cos(angle_M), 1.52*ae*sin(angle_M),'rO')
axis equal

title('Траектория КА')
xlabel('x, м')
ylabel('y, м')

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

hold off;
