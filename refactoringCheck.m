clc;
clear;
x0=zeros([1, 10]);
t_start = juliandate(2022,1,1);

n=0;
angle=0.75;
x0(9)=n+angle;

case_traj = 2;

N=1350;
m0=367;
eta=0.45;

terminal_state = 't';
integration_acc=1e-16;
calculate_condition = 0;

ae = 149597870700;
mug_0 = 132712.43994*(10^6)*(10^(3*3));
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;
T_mars_days = 365.256363004*1.8808476;

r_unit=ae;
V_unit=sqrt(mug_0/ae);
T_unit = T_earth/(2*pi);

planet_start = 'Earth';
planet_end = 'Mars';

[r0, V0] = planetEphemeris(t_start,'SolarSystem',planet_start,'430');

eul = [0 pi/4 0];
rotmZYX = eul2rotm(eul);

r0 = [rotmZYX*r0'/ae; 0]*1e+03;
V0 = [rotmZYX*V0'/V_unit; 0]*1e+03;



mug=1;

if terminal_state == 's'
    s_f=x0(9)*2*pi;
elseif terminal_state == 't'
    s_f=1.5*x0(9)*2*pi;
end
t_end_0=x0(9)*365.256363004;
days2sec=24*3600;
tspan=linspace(0, t_end_0*days2sec, 1000);
[rr_cont, Jt_cont, C_cont, evaluation_time_cont, dr_cont, dV_cont, PR, PV, t_cont] =...
    checkContinuation(t_start, t_end_0, tspan, case_traj,planet_end,eta, n);


%%
%задаем начальные условия
phi0=0;
u0=rToU(r0,phi0);
w0=vFromV(V0, r0, mug, phi0);
h0=-mug/(u0'*u0+4*(w0'*w0));

pr0=PR(1,:)';
pv0=PV(1,:)';

px0 = prpV2pupw(u0,w0,pr0,pv0)';
pu0=px0(1:4);
pw0=px0(5:8);
%pu0 = [0; 0; 0; 0];
%pw0 = [0; 0; 0; 0];
a_check = a_reactive(u0,w0,pu0,pw0);
if norm(a_check(1:3)-pv0) > 1e-14
    disp('a=/=pv problem');
end
tau0=2*(u0'*w0)/sqrt(-2*h0);
y0 = [u0; w0; pu0; pw0; tau0];

%t_start_fix=T_unit*(y0(10)-2*(y0(1:4)*y0(5:8)')/sqrt(-2*(y0(9)')))/(24*60*60);
int_s0sf = linspace(0, s_f, 1e+3);
time0 = tic;
%options = odeset('Events', @(s, y) eventIntegrationTraj(s, y, tf));
acc=integration_acc;
options = odeset('AbsTol',acc);
options = odeset(options,'RelTol',acc);
%максимальное время интегрирования
maxtime=1000;
if terminal_state == 's'
    options = odeset(options, 'Events',@(s, y) eventIntegrationTraj(s, y, time0, maxtime));
elseif terminal_state == 't'
    options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopTime(s, y, time0,maxtime, t_end_0));
end

dydy0=reshape(eye(20),[1 400]);
if calculate_condition == 1
    [s,Y] = ode113(@(s,y) integrateTraectoryWithVariations(s,y),int_s0sf,[y0, dydy0], options);
    y=Y(:,1:20);
else
    [s,y] = ode113(@(s,y) integrateTraectory(s,y),int_s0sf,y0, options);
end
uu = y(:, 1:4);
rr=zeros(length(uu),4);

t=zeros(length(uu),1);
HH=zeros(length(uu),1);
P_TR=zeros(length(uu),2);
VV=zeros(length(uu),4);
a_ks=zeros(length(uu),4);
for i = 1:length(uu)
    u = uu(i,:)';
    r=KS(u);
    rr(i,:)=r;
    L=L_KS(u);
    u2=norm(u)^2;
    w=y(i, 5:8)';
    h=-mug/(u'*u+4*w'*w);
    tau=y(i ,17)';
    pu=y(i, 9:12)';
    pw=y(i, 13:16)';
    f_ortdgduv=get_ortdgduv(u,w);
    p_tr=[[pu;pw]'*f_ortdgduv]';
    %ph=y(i, 19)';
    %ptau=y(i, 20)';
    dtds=u2/sqrt(-2*h);
    aa_ks=a_reactive(u,w,pu,pw);
    a_ks(i, :)=aa_ks/(ae/sqrt(mug_0)).^2;

    V = 2*sqrt(-2*h)*L*w/(u2);
    VV(i, :)=V;
    t(i) = T_unit*(tau-2*(u'*w)/sqrt(-2*h));
    H=calculateHamiltonian(u,w,pu,pw);
    HH(i)=H;
    P_TR(i,:)=p_tr;
end
t = t - t(1);
Jt = integrateFunctional(t, y, eta);

if terminal_state == 's'
    t_end = T_unit*(tau-2*(u'*w)/sqrt(-2*h))/(24*60*60);
elseif terminal_state == 't'
    t_end = t_end_0;
end




%Проверка "на глаз"
figure(1);
plot3(0, 0, 0, 'y--o')
set(gca,'FontSize',14)
hold on;

t0 = t_start;
t_orbit = linspace(t0,t0+T_earth/(24*3600), 1000);
earth_traj = planetEphemeris(t_orbit','SolarSystem','Earth','430');
earth_traj=earth_traj*1e+03/ae;
%Для KS
earth_traj_New = arrayfun(@(x,y,z)rotmZYX*[x, y, z]', earth_traj(:, 1),earth_traj(:, 2),earth_traj(:, 3),'UniformOutput',false);
earth_traj_New = cell2mat(earth_traj_New')';

t_orbit = linspace(t0,t0+T_mars/(24*3600), 1000);
mars_traj = planetEphemeris(t_orbit','SolarSystem',planet_end,'430');
mars_traj=mars_traj*1e+03/ae;

%Для KS
mars_traj_New = arrayfun(@(x,y,z)rotmZYX*[x, y, z]', mars_traj(:, 1),mars_traj(:, 2),mars_traj(:, 3),'UniformOutput',false);
mars_traj_New = cell2mat(mars_traj_New')';

plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k')
plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'r')

[mars_r_f, mars_v_f]=planetEphemeris([t_start, t_end],'SolarSystem',planet_end,'430');
mars_r_f=mars_r_f'*1e+03;
mars_v_f=mars_v_f'*1e+03;

rr_old = arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', rr(:, 1),rr(:, 2),rr(:, 3),'UniformOutput',false);
rr_old = cell2mat(rr_old')';

VV_old = arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', VV(:, 1),VV(:, 2),VV(:, 3),'UniformOutput',false);
VV_old = cell2mat(VV_old')';


plot3(rr_old(:, 1), rr_old(:, 2), rr_old(:, 3), 'b', 'LineWidth', 2.5);


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
for i = idxes
    plot3([rr_old(i, 1), rr_old(i, 1)+a_scale*a_ks_old(i, 1)], [rr_old(i, 2), rr_old(i, 2)+...
        a_scale*a_ks_old(i, 2)],[rr_old(i, 3), rr_old(i, 3)+a_scale*a_ks_old(i, 3)],'k')
end

plot3(rr_old(end, 1), rr_old(end, 2), rr_old(end, 3),'bO')
plot3(mars_r_f(1)/ae, mars_r_f(2)/ae,mars_r_f(3)/ae,'rO')

a_cont_old= arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', PV(:, 1),PV(:, 2),PV(:, 3),'UniformOutput',false);
a_cont_old = cell2mat(a_cont_old')';

a_scale_cont=3e-01/mean(vecnorm(PV, 2, 2));
%a_scale=0;
d = 24*3600;
idxes=1;
for i=1:ceil(t_cont(end)/d)
    ix = find(t_cont>d*i*10, 1);
    idxes=[idxes, ix];
end    
for i = idxes
    plot3([rr_cont(i, 1)/ae, rr_cont(i, 1)/ae+a_scale_cont*a_cont_old(i, 1)], [rr_cont(i, 2)/ae, rr_cont(i, 2)/ae+...
        a_scale_cont*a_cont_old(i, 2)],[rr_cont(i, 3)/ae, rr_cont(i, 3)/ae+a_scale_cont*a_cont_old(i, 3)],'m')
end

plot3(rr_cont(:, 1)/ae, rr_cont(:, 2)/ae, rr_cont(:, 3)/ae, 'g', 'LineWidth', 2.5);
axis equal

title('Траектория КА')
xlabel('x, a.e.')
ylabel('y, a.e.')
zlabel('z, a.e.')
view(0,90)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
hold off;