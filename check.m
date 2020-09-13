clear;
clc;
symbolic_Jacob
n = 0;
angle = 6*pi/6;

x=[0 1e-20 0 0]

ae = 149597870700;
mug = 132712.43994*(10^6)*(10^(3*3));

r0 = [1*ae 0 ]';
V0 = [0 (mug/(1*ae))^(1/2) ]';

t0=0;
u0 = [0 0]';
h0=(norm(V0)^2)/2-mug/norm(r0);

u0(1) = sqrt((norm(r0)+r0(1))/2);
u0(2) = r0(2)/(2*u0(1));

L = [[u0(1) -u0(2)];
    [u0(2) u0(1)]];  
v0 = L'*V0/(2*sqrt(-2*h0));
y0 = cat(1, u0, v0, x', t0)';

sf = (n*2*pi+angle)*1.5;
int_s0sf = linspace(0, sf, (n+1)*1e+4);
options = odeset('Events', @(s, y) eventIntegrationTraj(s, y, angle, n));
options = odeset(options,'AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);
%Интегрируем, используя сопряженные переменные из fmincon

[s,y] = ode113(@(s,y) integrateTraectory(s, y, symF),int_s0sf,y0, options);

u = y(:, 1:2);
r=zeros(length(u),2);
for i = 1:length(u)
    rr = u(i,:);
    L = [[rr(1) -rr(2)];
    [rr(2) rr(1)]];
    r(i,:)=L*rr';
end

%Проверка "на глаз"
plot(0, 0,'y--o')
hold on;
th = 0:pi/50:2*pi;
plot(ae*cos(th),ae*sin(th),'k');
plot(1.52*ae*cos(th),1.52*ae*sin(th),'r');
plot(r(:, 1), r(:, 2),'b')
plot(r(end, 1), r(end, 2),'bO')
plot(1.52*ae*cos(angle), 1.52*ae*sin(angle),'rO')
axis equal
hold off;