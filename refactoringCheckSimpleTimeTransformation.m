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

r0 = [rotmZYX*r0'/ae;]*1e+03;
V0 = [rotmZYX*V0'/V_unit;]*1e+03;



mug=1;

if terminal_state == 's'
    s_f=x0(9)*2*pi;
elseif terminal_state == 't'
    s_f=1.5*x0(9)*2*pi;
end
t_end_0=x0(9)*365.256363004;
phi = x0(10)*2*pi;
days2sec=24*3600;
tspan=linspace(0, t_end_0*days2sec, 1000);
[rr_cont, Jt_cont, C_cont, evaluation_time_cont, dr_cont, dV_cont, PR, PV, t_cont] =...
    checkContinuation(t_start, t_end_0, tspan, case_traj,planet_end,eta, n);


%%
%задаем начальные условия
dtds0=norm(r0);
k   = 1;

rs0 = r0;
Vs0 = V0*dtds0;
%tau0=(r0'*Vs0)/norm(r0)/sqrt(-2*h0);
t0  = 0;
pr0 = PR(1,:)';
pv0 = PV(1,:)';
%Определение начальных параметров

%%
px=prpV2pXpVXSimple(rs0,Vs0,pr0,pv0);
%px=[0, 0, 0, 0, 0, 0];
%prs0=px(1:3)';
%pVs0=px(4:6)';

prs0=k*(pr0-(pv0'*V0)*r0/norm(r0)^2);
pVs0=k*(pv0/norm(r0));
y0 = [rs0; Vs0; prs0; pVs0;t0;];

%t_start_fix=T_unit*(y0(10)-2*(y0(1:4)*y0(5:8)')/sqrt(-2*(y0(9)')))/(24*60*60);
acc=integration_acc;
options = odeset('AbsTol',acc);
options = odeset(options,'RelTol',acc);
time0 = tic;
maxtime=1000;
if terminal_state == 's'
    options = odeset(options, 'Events',@(s, y) eventIntegrationTraj(s, y, time0, maxtime));
elseif terminal_state == 't'
    options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopTimeSimple(s, y, time0,maxtime, t_end_0));
end


%максимальное время интегрирования

[t,y] = ode113(@(t,y) integrateTraectorySimpleTimeTransformation(t,y),tspan/T_unit,y0, options);

XX = y(:, 1:3);
PVX = y(:, 10:12);
rr=zeros(length(XX),3);

aa=zeros(length(XX),3);
for i = 1:length(XX)
    r=XX(i,:)';
    pVX=PVX(i,:);
    Vs=y(i, 4:6)';
    rr(i,:)=r;
    dtds=norm(r);
    a=pVX*dtds;
    aa(i,:)=a;
end
t = t - t(1);


t_end = t_end_0;





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

plot3(rr_old(:, 1), rr_old(:, 2), rr_old(:, 3), 'b', 'LineWidth', 2.5);


aa_old= arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', aa(:, 1),aa(:, 2),aa(:, 3),'UniformOutput',false);
aa_old = cell2mat(aa_old')';

a_scale=3e-01/mean(vecnorm(aa, 2, 2));
%a_scale=0;
d = 24*3600;
idxes=1;
for i=1:ceil(t(end)/d)
    ix = find(t>d*i*10, 1);
    idxes=[idxes, ix];
end    
for i = idxes
    plot3([rr_old(i, 1), rr_old(i, 1)+a_scale*aa_old(i, 1)], [rr_old(i, 2), rr_old(i, 2)+...
        a_scale*aa_old(i, 2)],[rr_old(i, 3), rr_old(i, 3)+a_scale*aa_old(i, 3)],'k')
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
plot3(rr_cont(end, 1)/ae, rr_cont(end, 2)/ae, rr_cont(end, 3)/ae, 'gO', 'LineWidth', 2.5);
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