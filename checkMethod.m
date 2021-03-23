function [dr,dv, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks] = checkMethod(t_start,psi,rad, UorR,direction,modifier_p,modifier_f, x0, eta)
%UNTITLED9 Summary of this function goes here
%   Вычисляет невязку в зависимости от входных параметров
%условия на fmincon
%ЗАДАЧА ПРОЛЁТА case_traj=1; ЗАДАЧА сопровождения case_traj=2;
case_traj=2;
%Начальные условия
%x0=zeros([1, 11]);

A = [];
b = [];
Aeq = [];
beq = [];
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

% modifier_p=1e-04;
% modifier_f=1e+04;
modifier_b=1e+13;

s_a = psi-rad;
s_b = psi+rad;

x0(11)=psi;


lb = -[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0]*modifier_b;
ub = -lb;

lb(11) = s_a;
ub(11) = s_b;

lb(12) = 0;
ub(12) = 1;
%домножаем на коэффициент 1е-12, чтобы fmincon работал с более крупными
%величинами и не выдавал лишних ворнингов
tic;
fun=@(x)fun2min([x(1:10)*modifier_p x(11), x(12)], case_traj, t_start, r0, V0, planet_end, modifier_f, UorR,direction);

options = optimoptions('fmincon','UseParallel', true);
options = optimoptions(options, 'Display', 'iter');
options = optimoptions(options, 'OptimalityTolerance', 1e-15);
options = optimoptions(options, 'MaxFunctionEvaluations', 1e+10);
options = optimoptions(options, 'StepTolerance', 1e-12);
options = optimoptions(options, 'ConstraintTolerance', 1e-12);
options = optimoptions(options, 'MaxIterations', 1500);
options = optimoptions(options,'OutputFcn',@myoutput);

[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub,[], options);
toc
px = x(1:10)*modifier_p;
s_f = x(11)*2*pi;
phi = x(12)*2*pi;
%задаем начальные условия
%options = optimoptions(options,'OutputFcn',@myoutput);
%options = optimoptions(options, 'Algorithm', 'sqp');

t0=0;
u0 = [0 0 0 0]';
h0=(norm(V0)^2)/2-mug/norm(r0);

u0 = rToU(r0, 0);
L = L_KS(u0); 
v0 = vFromV(V0,r0,mug,0);

tau0=getEccentricAnomaly(r0(1:3),V0(1:3),mug);
y0 = cat(1, u0, v0, 0, tau0,  px')';
n = 1;
int_s0sf = linspace(0, s_f, (n+1)*1e+4);
time0 = tic;
%options = odeset('Events', @(s, y) eventIntegrationTraj(s, y, tf));
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);
options = odeset(options, 'Events',@(s, y) eventIntegrationTraj(s, y, time0));
%Интегрируем, используя сопряженные переменные из fmincon

[s,y] = ode113(@(s,y) integrateTraectory(s,y,h0),int_s0sf, y0, options);
%Jt = integrateFunctional(s, y, eta, h0);
%functional = Jt(end);
Jt = integrateFunctional(s, y, eta, h0);
t_start_fix=T_unit*(y(1, 10)-2*(y(1, 1:4)*y(1, 5:8)')/sqrt(-2*(y(1, 9)'+h0)))/(24*60*60);
uu = y(:, 1:4);
rr=zeros(length(uu),4);

t=zeros(length(uu),1);
VV=zeros(length(uu),4);
a_ks=zeros(length(uu),4);
for i = 1:length(uu)
    u = uu(i,:)';
    r=KS(u);
    rr(i,:)=r;
    L=L_KS(u);
    u2=norm(u)^2;
    v=y(i, 5:8)';
    h=y(i, 9)'+h0;
    tau=y(i ,10)';
    pv=y(i, 15:18)';
    ph=y(i, 19)';
    ptau=y(i, 20)';
    dtds=u2/sqrt(-2*h);
    aa_ks=L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*(u2)*u/((-2*h)^(3/2)))/dtds;
    a_ks(i, :)=aa_ks/(ae/sqrt(mug_0)).^2;

    V = 2*sqrt(-2*h)*L*v/(u2);
    VV(i, :)=V;
    t(i) = T_unit*(tau-2*(u'*v)/sqrt(-2*h));
end
t = t - t(1);
t_end = T_unit*(tau-2*(u'*v)/sqrt(-2*h))/(24*60*60)-t_start_fix;


[mars_r_f, mars_v_f]=planetEphemeris([t_start, t_end],'SolarSystem',planet_end,'430');
mars_r_f=rotmZYX*mars_r_f'*1e+03;
mars_v_f=rotmZYX*mars_v_f'*1e+03;

dr=norm(ae*rr(end, 1:3)-mars_r_f(1:3)');
dv=norm(V_unit*VV(end, 1:3)-mars_v_f(1:3)');
C=norm(x)*norm(grad)/fval;
end

