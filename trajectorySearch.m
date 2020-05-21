function [functional, dis, s, y] = trajectorySearch(n,angle,case_traj)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%clear;
%clc;
%������� �� fmincon
%������ ���˨�� case_traj=1; ������ ������������� case_traj=2;
%case_traj=2;
%���������� ������
%n = 4;
%angle = 6*pi/6;
%��������� �������
x0=[0 0 0 0 0];
A = [];
b = [];
Aeq = [];
beq = [];
ae = 149597870700;
mug = 132712.43994*(10^6)*(10^(3*3));

rf = 1.52*ae*[cos(angle) sin(angle)];
Vf = ((mug/(1.52*ae))^(1/2))*[cos(angle+pi/2) sin(angle+pi/2)];
lb = -[1, 1, 1, 1, 1e-4]*1000;
ub = [1, 1, 1, 1, 1e-4]*1000;
%��������� �� ����������� 1�-12, ����� fmincon ������� � ����� ��������
%���������� � �� ������� ������ ���������
options = optimoptions ('fmincon','Display', 'none');
fun=@(x)fun2min(x*1e-12, rf, Vf, case_traj, n, angle);
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[],options)*1e-12;
%������ ��������� �������
r0 = [1*ae 0 ]';
V0 = [0 (mug/(1*ae))^(1/2) ]';

t0=0;
u0 = [0 0]';
h0=(norm(V0)^2)/2-mug/norm(r0);

u0(2) = sqrt((norm(r0)-r0(1))/2);
u0(1) = sqrt(norm(r0)-u0(2)^2);

L = [[u0(1) -u0(2)];
    [u0(2) u0(1)]];  
v0 = L'*V0/(2*sqrt(-2*h0));
y0 = cat(1, u0, v0, h0, x', t0)';

sf = (n*2*pi+angle)*1.5;
int_s0sf = linspace(0, sf, (n+1)*1e+4);
options = odeset('Events', @(s, y) eventIntegrationTraj(s, y, angle, n));
options = odeset(options,'AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);
%�����������, ��������� ����������� ���������� �� fmincon

[s,y] = ode113(@(s,y) integrateTraectory(s,y,mug),int_s0sf,y0, options);
functional = integrateFunctional(s, y);

u = y(:, 1:2);
r=zeros(length(u),2);
for i = 1:length(u)
    rr = u(i,:);
    L = [[rr(1) -rr(2)];
    [rr(2) rr(1)]];
    r(i,:)=L*rr';
end


%����� ��������� ������� ��� ������ ������
if case_traj == 1
    pv=y(end, 8:9);
    r_end=r(end,:);
    dis = norm((rf-r_end)/norm(r0))^2 + (norm(pv)^2)*1e+20;
elseif case_traj == 2
    v = y(end, 3:4)';
    h = y(end, 5);
    Lend = [[u(end, 1) -u(end, 2)];[u(end, 2) u(end, 1)]];
    V = 2*sqrt(-2*h)*Lend*v/(norm(u(end,:))^2);
    r_end=r(end,:);
    dis = norm((rf-r_end)/norm(r0))^2 + (norm((V'-Vf)/norm(V0))^2);
end
end

