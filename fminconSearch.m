%clearvars -except symF
%clc;
%5if exist('symF','var') ~= 1
%    symbolic_Jacob
%end
n=1;
t0=juliandate(2001,12,1);
tf=t0+365*n;

N=1350;
m0=367;
eta=0.45;
%������� �� fmincon
%������ ���˨�� case_traj=1; ������ ������������� case_traj=2;
case_traj=2;
%��������� �������
ae = 149597870700;
mug_0 = 132712.43994*(10^6)*(10^(3*3));
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;

days2sec = 24*3600;

r_norm=ae;
V_unit=sqrt(mug_0/ae);
T_unit = T_earth/(2*pi);

[r0, V0] = planetEphemeris(t0,'SolarSystem','Earth','430');
[rf, Vf] = planetEphemeris(tf,'SolarSystem','Mars','430');
r0=r0*1e+03;
V0=V0*1e+03;
rf=rf*1e+03;
Vf=Vf*1e+03;
z0 = [0 0 0 0 0 0];
y0 = cat(2,r0,V0,z0,...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(eye(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(eye(3),[1,9])...
    )';
tspan = linspace(0, (tf-t0)*days2sec, (n+1)*1e+4);
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);
[t, y] = ode45(@(t,y)integrateTraectory(t,y),tspan, y0,options);

r_end =y(end, 1:3);
V_end =y(end, 4:6);
pV_end=y(end, 7:9);

if case_traj == 1
    b = cat(2,r_end-rf, pV_end);
elseif  case_traj == 2
    b = cat(2,r_end-rf, V_end-Vf);
end

tic;
[tau, z] = ode45(@(tau, z)integrateExternal(tau, z, b, case_traj, tspan, r0,V0), [0 1], z0);
toc
y_found = cat(2,r0,V0,z(end,:),...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(eye(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(eye(3),[1,9])...
    )';

%options = odeset('Events', @(s, y) eventIntegrationTraj(s, y, tf));
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);
%�����������, ��������� ����������� ���������� �� fmincon

[s,y] = ode113(@(s,y) integrateTraectory(s,y),tspan, y0, options);

rr=y(:,1:3);
%Jt = integrateFunctional(s, y, eta, h0);
%functional = Jt(end);

% figure(2);
% plot(t/(24*3600), vecnorm(a, 2, 2)*1e+03, 'LineWidth', 3);
% %title('����������� ��������� ���� ���� �� �������')
% xlabel('t, �����, ���','FontSize',14)
% ylabel('���������� ���������, ��/c^2','FontSize',14)
% box off;
% set(gca,'FontSize',14)
% 
% figure(4);
% m=massLP(Jt, m0, N);
% plot(t/(24*3600), m, 'LineWidth', 3);
% title('����������� ����� �� �������')
% xlabel('t, �����, ���')
% ylabel('m, �����, ��')
% box off;
% set(gca,'FontSize',14)

%�������� "�� ����"
figure(1);
%plot3(0, 0, 0, 'y--o')
set(gca,'FontSize',14)
hold on;
th = 0:pi/50:2*pi;

t_orbit = linspace(t0,t0+T_earth/(24*3600), 1000);
earth_traj = planetEphemeris(t_orbit','SolarSystem','Earth','430');
earth_traj=earth_traj*1e+03;

t_orbit = linspace(t0,t0+T_mars/(24*3600), 1000);
mars_traj = planetEphemeris(t_orbit','SolarSystem','Mars','430');
mars_traj=mars_traj*1e+03;

%plot(cos(th),sin(th),'k');
%plot(1.52*cos(th),1.52*sin(th),'r');
plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k')
plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'r')

mars_r_f=planetEphemeris(tf,'SolarSystem','Mars','430','AU');

plot3(rr(:, 1), rr(:, 2), rr(:, 3), 'b', 'LineWidth', 2.5);
%a_scale=3e-01/mean(vecnorm(a, 2, 2));
%a_scale=0;
d = 24*3600;
idxes=1;
for i=1:ceil(t(end)/d)
    ix = find(t>d*i*10, 1);
    idxes=[idxes, ix];
end    
for i = idxes
    %plot3([rr(i, 1), rr(i, 1)+a_scale*a(i, 1)], [rr(i, 2), rr(i, 2)+a_scale*a(i, 2)],[rr(i, 3), rr(i, 3)+a_scale*a(i, 3)],'k')
end

%plot3(rr(end, 1), rr(end, 2), rr(end, 3),'bO')
%plot3(mars_r_f(1), mars_r_f(2),mars_r_f(3),'rO')
axis equal

%title('���������� ��')
xlabel('x, a.e.')
ylabel('y, a.e.')

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
hold off;

[mars_r_f, mars_v_f]=planetEphemeris(tf,'SolarSystem','Mars','430');

%disp(['������ ����� ', num2str(m(1)-m(end)), '��'])
disp(['������� ���������� ', num2str(norm(rr(end, 1:3)*r_norm-mars_r_f*1e+03),'%10.2e\n'),',�'])
%disp(['������� �������� ', num2str(norm(VV(end, 1:3)*V_unit-mars_v_f*1e+03),'%10.2e\n'),',�/�'])
% ������������� ����� ���������������
%disp(['������������� ����� ��������������� ', num2str(norm(x)*norm(grad)/fval,'%10.2e\n')])
% ���������� ����� ���������������
%disp(['���������� ����� ��������������� ', num2str(1/norm(grad),'%10.2e\n')])