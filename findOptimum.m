
mug = 132712.43994*(10^6)*(10^(3*3));
ae = 149597870700;

r0 = [1*ae 0 ]';
V0 = [0 (mug/(1*ae))^(1/2) ]';

u0 = [0 0]';
h0=(norm(V0)^2)/2-mug/norm(r0);

u0(2) = sqrt((norm(r0)-r0(1))/2);
u0(1) = sqrt(norm(r0)-u0(2)^2);

L = [[u0(1) -u0(2)];
    [u0(2) u0(1)]];  
v0 = L'*V0/(2*sqrt(-2*h0));
pu0=[0 0]'*1e-13;
pv0=[0 0]'*1e-12;
ph0=0'*1e-15;
t0=0;
y0 = cat(1, u0, v0, h0, pu0, pv0, ph0, t0)';
%Определяем tf
T=2*pi*sqrt((1*ae)^3/mug);
%tf=3*T/12;
sf = 4*pi;
angle = 3*pi/2;
rf = 1.52*ae*[cos(angle) sin(angle)];
n = 1;
options = odeset('Events', @(s, y) eventIntegrationTraj(s, y, angle, n));
options = odeset(options,'AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);
int_s0sf = linspace(0, sf, n*1e+4);
[s,y] = ode113(@(s,y) integrateTraectory(s,y,mug),int_s0sf,y0, options);
u = y(:, 1:2);
r=zeros(length(u),2);
for i = 1:length(u)
    x = u(i,:);
    L = [[x(1) -x(2)];
    [x(2) x(1)]];
    r(i,:)=L*x';
end
plot(0, 0,'y--o')
hold on;
th = 0:pi/50:2*pi;
plot(ae*cos(th),ae*sin(th),'k');
plot(r(:, 1), r(:, 2),'b')
plot(r(end, 1), r(end, 2),'bO')
axis equal
hold off;