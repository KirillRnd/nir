function [functional, dis, s, y] = trajectorySearch(n, angle, rad, case_traj, symF)
%trajectorySearch Summary of this function goes here
%������ ��� fmincon, ����� ��������� ���������� � ���� ����� ����������

%��������� �������
x0=[0 0 0 0 0 0 0 0 0 0 0];
A = [];
b = [];
Aeq = [];
beq = [];
ae = 149597870700;
mug = 132712.43994*(10^6)*(10^(3*3));
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;

modifier=1e-8;
modifier_p=1e-15; 
tf_a = T_earth*(n + angle-rad)*modifier;
tf_b = T_earth*(n + angle+rad)*modifier;

x0(11)=T_earth*(n + angle)*modifier;

n_M = floor((x0(11)/modifier)/T_mars);
angle_M = (x0(11)/modifier)/T_mars-n_M;
t_Mars_0 = (angle-angle_M-0.03)*T_mars;
lb = -[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]*10000;
ub = -lb;

lb(11) = tf_a;
ub(11) = tf_b;
%��������� �� ����������� 1�-12, ����� fmincon ������� � ����� ��������
%���������� � �� ������� ������ ���������

fun=@(x)fun2min([x(1:10)*modifier_p x(11)/modifier], case_traj, symF, t_Mars_0);
x = fmincon(fun, x0, A, b, Aeq, beq, lb, ub);
px = x(1:10)*modifier_p;
tf = x(11)/modifier;
%������ ��������� �������
r0 = [1*ae 0 0 0]';
V0 = [0 (mug/(1*ae))^(1/2) 0 0]';

t0=0;
u0 = [0 0 0 0]';
h0=(norm(V0)^2)/2-mug/norm(r0);

u0(4) = 0;
u0(1) = sqrt((norm(r0)+r0(1))/2);
u0(2) = r0(2)/(2*u0(1));
u0(3) = r0(3)/(2*u0(1));

L = [[u0(1) -u0(2) -u0(3) u0(4)];
    [u0(2) u0(1) -u0(4) -u0(3)];
    [u0(3) u0(4) u0(1) u0(2)];
    [u0(4) -u0(3) u0(2) -u0(1)]]; 
v0 = L'*V0/(2*sqrt(-2*h0));
tau0=0;
y0 = cat(1, u0, v0, h0, tau0,  px')';


n = floor(tf/T_earth);
angle = (tf/T_earth-n)*2*pi;

sf = (n*2*pi+angle)*1.5;
int_s0sf = linspace(0, sf, (n+1)*1e+4);
options = odeset('Events', @(s, y) eventIntegrationTraj(s, y, tf));
options = odeset(options,'AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);
%�����������, ��������� ����������� ���������� �� fmincon

[s,y] = ode113(@(s,y) integrateTraectory(s,y,symF),int_s0sf,y0, options);
functional = integrateFunctional(s, y)

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


n_M = floor((tf+t_Mars_0)/T_mars);
angle_M = ((tf+t_Mars_0)/T_mars-n_M)*2*pi;

rf = 1.52*ae*[cos(angle_M) sin(angle_M) 0 0];
Vf = ((mug/(1.52*ae))^(1/2))*[cos(angle_M+pi/2) sin(angle_M+pi/2) 0 0];

%����� ��������� ������� ��� ������ ������
if case_traj == 1
    pv=y(end, 14:17);
    r_end=r(end,:);
    dis = norm((rf-r_end)/norm(r0))^2 + (norm(pv)^2)*1e+20;
elseif case_traj == 2
    v = y(end, 5:8)';
    h = y(end, 9);
    Lend = [[u(end, 1) -u(end, 2) -u(end, 3) u(end, 4)];
    [u(end, 2) u(end, 1) -u(end, 4) -u(end, 3)];
    [u(end, 3) u(end, 4) u(end, 1) u(end, 2)];
    [u(end, 4) -u(end, 3) u(end, 2) -u(end, 1)]];
    V = 2*sqrt(-2*h)*Lend*v/(norm(u(end,:))^2);
    r_end=r(end,:);
    dis = norm((rf-r_end)/norm(r0))^2 + (norm((V'-Vf)/norm(V0))^2);
end
end

