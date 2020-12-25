function [dr,dv,C] = checkMethod(t_start,phi)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
%clearvars -except symF
%clc;
%5if exist('symF','var') ~= 1
%    symbolic_Jacob
%end
% t_start = juliandate(2001,12,1);
% N=1350;
% m0=367;
% eta=0.45;
%условия на fmincon
%ЗАДАЧА ПРОЛЁТА case_traj=1; ЗАДАЧА сопровождения case_traj=2;
case_traj=2;
%Начальные условия
x0=[0 0 0 0 0 0 0 0 0 0 0];
A = [];
b = [];
Aeq = [];
beq = [];
ae = 149597870700;
mug_0 = 132712.43994*(10^6)*(10^(3*3));
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;

r_norm=ae;
V_unit=sqrt(mug_0/ae);
T_unit = T_earth/(2*pi);

[r0, V0] = planetEphemeris(t_start,'SolarSystem','Earth','430');

r0 = [r0/ae, 0]'*1e+03;
V0 = [V0/V_unit, 0]'*1e+03;

mug=1;

n=2;
angle=0.5;
rad=pi/48;
%d_mars=-0.25;

modifier_p=1e-05;


s_a = (phi-rad);
s_b = (phi+rad);

x0(11)=(phi);

lb = -[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]*1e+10;
ub = -lb;

lb(11) = s_a;
ub(11) = s_b;
%домножаем на коэффициент 1е-12, чтобы fmincon работал с более крупными
%величинами и не выдавал лишних ворнингов
tic;
fun=@(x)fun2min([x(1:10)*modifier_p x(11)], case_traj, t_start, r0, V0);

options = optimoptions('fmincon','UseParallel', true);
%options = optimoptions(options, 'Display', 'iter');
options = optimoptions(options, 'OptimalityTolerance', 1e-10);
options = optimoptions(options, 'MaxFunctionEvaluations', 1e+10);
options = optimoptions(options, 'StepTolerance', 1e-10);
options = optimoptions(options, 'ConstraintTolerance', 1e-10);
options = optimoptions(options,'OutputFcn',@myoutput);
%options = optimoptions(options, 'Algorithm', 'sqp');

[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub,[], options);
toc
px = x(1:10)*modifier_p;
s_f = x(11);
%задаем начальные условия

t0=0;
u0 = [0 0 0 0]';
h0=(norm(V0)^2)/2-mug/norm(r0);

u0(4) = 0;
u0(1) = sqrt((norm(r0)+r0(1))/2);
u0(2) = r0(2)/(2*u0(1));
u0(3) = r0(3)/(2*u0(1));

L = L_KS(u0); 
v0 = L'*V0/(2*sqrt(-2*h0));
tau0=0;
y0 = cat(1, u0, v0, 0, tau0,  px')';

int_s0sf = linspace(0, s_f, (n+1)*1e+4);
%options = odeset('Events', @(s, y) eventIntegrationTraj(s, y, tf));
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);
%Интегрируем, используя сопряженные переменные из fmincon

[s,y] = ode113(@(s,y) integrateTraectory(s,y,h0),int_s0sf, y0, options);
%Jt = integrateFunctional(s, y, eta, h0);
%functional = Jt(end);

uu = y(:, 1:4);
rr=zeros(length(uu),4);
%a=zeros(length(uu),4);

t=zeros(length(uu),1);
VV=zeros(length(uu),4);
for i = 1:length(uu)
    u = uu(i,:)';
    r=KS(u);
    rr(i,:)=r;
    L=L_KS(u);
    u2=norm(u)^2;
    v=y(i, 5:8)';
    h=y(i, 9)'+h0;
    tau=y(i ,10)';
%     pu=y(i, 11:14)';
%     pv=y(i, 15:18)';
%     ph=y(i, 19)';
%     ptau=y(i, 20)';
    %aa=L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*(rr'*rr)*rr/(-2*h)^(3/2));
    %res=symF(u,v,h,pu,pv,ph,ptau);
    %dvds=res(5:8);
    %dhds=res(9);
    V = 2*sqrt(-2*h)*L*v/(u2);
    VV(i, :)=V;
    %a(i, :)=((-2*h/(norm(r)^2))*(2*(L_KS(v)*v+L_KS(u)*dvds)-(2*u'*v/(sqrt(-2*h)) + norm(r)*dhds/((-2*h)^(3/2)))*V)+mug*r/(norm(r)^3))/(ae/sqrt(mug_0)).^2;
    
    %a(i, :)=KS(aa);
    t(i) = T_unit*(tau-2*(u'*v)/sqrt(-2*h));
end

t_end=t(end);

[mars_r_f, mars_v_f]=planetEphemeris([t_start, t_end/(24*3600)],'SolarSystem','Mars','430');

dr=norm(rr(end, 1:3)*r_norm-mars_r_f*1e+03);
dv=norm(VV(end, 1:3)*V_unit-mars_v_f*1e+03);
C=norm(x)*norm(grad)/fval;
end

