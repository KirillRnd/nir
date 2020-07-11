function dis = fun2min(x, rf, Vf, case_traj, n, angle, symF)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

pu0=x(1:4)';
pv0=x(5:8)';
ph0=x(9);

mug = 132712.43994*(10^6)*(10^(3*3));
ae = 149597870700;

r0 = [1*ae 0 0 0]';
V0 = [0 (mug/(1*ae))^(1/2) 0 0]';

u0 = [0 0 0 0]';
h0 = (norm(V0)^2)/2-mug/norm(r0);

u0(4) = 0;
u0(1) = sqrt((norm(r0)+r0(1))/2);
u0(2) = r0(2)/(2*u0(1));
u0(3) = r0(3)/(2*u0(1));

L = [[u0(1) -u0(2) -u0(3) u0(4)];
    [u0(2) u0(1) -u0(4) -u0(3)];
    [u0(3) u0(4) u0(1) u0(2)];
    [u0(4) -u0(3) u0(2) -u0(1)]]; 

v0 = L'*V0/(2*sqrt(-2*h0));
%pu0=[1 0]'*1e-12;
%pv0=[0 0]'*1e-08;
%ph0=0';
t0 = 0;
pt0=0;
y0 = cat(1, u0, v0, h0, pu0, pv0, ph0, t0, pt0)';
%Îïðåäåëÿåì tf
T=2*pi*sqrt((1*ae)^3/mug);
%tf=3*T/12;
sf = (n*2*pi+angle)*1.5;
%angle = 3*pi/2;

options = odeset('Events', @(s, y) eventIntegrationTraj(s, y, angle, n));
options = odeset(options,'AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);

[s,y] = ode113(@(s,y) integrateTraectory(s,y, symF),[0 sf],y0, options);
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
%ÇÀÄÀ×À ÏÐÎË¨ÒÀ
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

