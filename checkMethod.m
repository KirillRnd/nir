function [dr,dv, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] = checkMethod(t_start,psi,rad, UorR,decreaseNonPsysical,modifier_p,modifier_f, x0, eta, case_traj,planet_end,display,terminal_state,integration_acc, calculate_condition, orbits, omega)
%UNTITLED9 Summary of this function goes here
%   Вычисляет невязку в зависимости от входных параметров
%условия на fmincon
%ЗАДАЧА ПРОЛЁТА case_traj=1; ЗАДАЧА сопровождения case_traj=2;

%Начальные условия
%x0=zeros([1, 11]);

A = [];
b = [];
Aeq = [];
beq = [];


ae = 149597870700;
mug_0 = 132712.43994*(10^6)*(10^(3*3));
T_earth = 365.256363004*3600*24;

r_unit=ae;
V_unit=sqrt(mug_0/ae);
T_unit = T_earth/(2*pi);
planet_start = 'Earth';

st.t = t_start;
st.planet = planet_start;
st.mode = orbits;
st.delta_omega = omega;

[r0, V0] = planetModel(st);

eul = [pi/12 pi/4 pi/12];
rotmZYX = eul2rotm(eul);

r0 = [rotmZYX*r0'/ae; 0]*1e+03;
V0 = [rotmZYX*V0'/V_unit; 0]*1e+03;
%V0 = V0 + 2700*V0/norm(V0)/V_unit;%гип избыток
%nonlcon = @(x)ubOrtPv(x, rToU(r0, 0));
nonlcon=[];
mug=1;
x0(1:8)=x0(1:8)/modifier_p;
% modifier_p=1e-04;
% modifier_f=1e+04;
modifier_b=1e+13;

%s_a = psi-rad;
s_a = psi;
s_b = psi+rad*2*pi;

%x0(11)=psi;


lb = -[1, 1, 1, 1, 1, 1, 1, 1, 0, 0]*modifier_b;
ub = -lb;

lb(9) = s_a/(2*pi);
ub(9) = s_b/(2*pi);


if strcmp(UorR,'u') || strcmp(UorR,'r')
    lb(10) = -1.0;
    ub(10) = 1.0;
elseif strcmp(UorR,'u_hat')
    lb(10) = x0(10);
    ub(10) = x0(10);
%     lb(10) = 0.0;
%     ub(10) = 1.0;
end
%домножаем на коэффициент 1е-12, чтобы fmincon работал с более крупными
%величинами и не выдавал лишних ворнингов
st_f2m.case_traj = case_traj;
st_f2m.t_start = t_start;
st_f2m.r0 = r0;
st_f2m.V0 = V0;
st_f2m.planet_end = planet_end;
st_f2m.modifier_f = modifier_f;
st_f2m.UorR = UorR;
st_f2m.decreaseNonPsysical = decreaseNonPsysical;
st_f2m.terminal_state = terminal_state;
st_f2m.integration_acc = integration_acc;
st_f2m.orbits = orbits;
st_f2m.omega = omega;
tic;
fun=@(x)fun2min([x(1:8)*modifier_p, x(9), x(10)], st_f2m);

options = optimoptions('fmincon','UseParallel', true);
if display == 1
    options = optimoptions(options, 'Display', 'iter');
end
options = optimoptions(options, 'OptimalityTolerance', 1e-10);
options = optimoptions(options, 'MaxFunctionEvaluations', 1e+10);
options = optimoptions(options, 'StepTolerance', 1e-10);
options = optimoptions(options, 'ConstraintTolerance', 1e-10);
options = optimoptions(options, 'MaxIterations', 250);
options = optimoptions(options, 'FiniteDifferenceType', 'central');
options = optimoptions(options, 'EnableFeasibilityMode', true);
%options = optimoptions(options, 'Algorithm', 'sqp');

options = optimoptions(options,'OutputFcn',@myoutput);

[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x)1, x0, A, b, Aeq, beq, lb, ub, fun, options);
evaluation_time = toc;
px = x(1:8)*modifier_p;
if terminal_state == 's'
    s_f=x(9)*2*pi;
elseif terminal_state == 't'
    s_f=1.5*x(9)*2*pi;
end
t_end_0=x(9)*365.256363004;
phi = x(10)*2*pi;
%задаем начальные условия
%options = optimoptions(options,'OutputFcn',@myoutput);
%options = optimoptions(options, 'Algorithm', 'sqp');
phi0=phi;
%phi0=0;
h0=(norm(V0)^2)/2-mug/norm(r0);

u0 = rToU(r0, phi0);
w0 = vFromV(V0,r0,mug,phi0);

%tau0=getEccentricAnomaly(r0(1:3),V0(1:3),mug);
%tau0=0;
tau0=2*u0'*w0/sqrt(-2*h0);
y0 = cat(1, u0, w0, px', tau0)';

%t_start_fix=T_unit*(y0(10)-2*(y0(1:4)*y0(5:8)')/sqrt(-2*(y0(9)')))/(24*60*60);
int_s0sf = linspace(0, s_f, 1e+3);
time0 = tic;
%options = odeset('Events', @(s, y) eventIntegrationTraj(s, y, tf));
acc=integration_acc;
options = odeset('AbsTol',acc);
options = odeset(options,'RelTol',acc);
%максимальное время интегрирования
maxtime=100;
if terminal_state == 's'
    options = odeset(options, 'Events',@(s, y) eventIntegrationTraj(s, y, time0, maxtime));
elseif terminal_state == 't'
    options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopTime(s, y, time0,maxtime, t_end_0));
end
%Интегрируем, используя сопряженные переменные из fmincon

ddeltady0=zeros([8,8]);
if calculate_condition == 1
    step_h=1e-6;
    for i=9:16
        y0_delta=zeros([1,17]);
        y0_delta(i)=step_h;
        [s,y] = ode113(@(s,y) integrateTraectory(s,y),int_s0sf,y0+y0_delta, options);
        u_end=y(end,1:4)';
        w_end=y(end,5:8)';
        px_end=y(end,9:16)';
        g=get_target_g(u_end,w_end);
        ortdgduv=get_ortdgduv(u_end,w_end);
        tr=[px_end'*ortdgduv]';
        p_plus=[g;tr];

        [s,y] = ode113(@(s,y) integrateTraectory(s,y),int_s0sf,y0-y0_delta, options);
        u_end=y(end,1:4)';
        w_end=y(end,5:8)';
        px_end=y(end,9:16)';
        g=get_target_g(u_end,w_end);
        ortdgduv=get_ortdgduv(u_end,w_end);
        tr=[px_end'*ortdgduv]';
        p_minus=[g;tr];
        partial=(p_plus-p_minus)/(2*step_h);
        ddeltady0(:,i-8)=partial;
    end

    C=cond(ddeltady0);
else
    C=1;
end

[s,y] = ode113(@(s,y) integrateTraectory(s,y),int_s0sf,y0, options);
%Jt = integrateFunctional(s, y, eta, h0);
%functional = Jt(end);
%на случай, если всё сломается
if length(s)<=100
    uu=zeros(length(s),4);
    rr=zeros(length(s),4);
    VV=zeros(length(s),4);
    t=zeros(length(s),1);
    VV=zeros(length(s),4);
    a_ks=zeros(length(s),4);
    Jt=zeros(length(s),1);
    C=1;
    dr=ae;
    dv=ae;
    t_end=T_earth;
    return
end
%t_start_fix=T_unit*(y(1, 10)-2*(y(1, 1:4)*y(1, 5:8)')/sqrt(-2*(y(1, 9)')))/(24*60*60);
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

st.t = [t_start, t_end];
st.planet = planet_end;
st.mode = orbits;
st.delta_omega = omega;

[mars_r_f, mars_v_f]=planetModel(st);
mars_r_f=rotmZYX*mars_r_f'*1e+03;
mars_v_f=rotmZYX*mars_v_f'*1e+03;

dr=norm(ae*rr(end, 1:3)-mars_r_f(1:3)');
dv=norm(V_unit*VV(end, 1:3)-mars_v_f(1:3)');

% if calculate_condition == 1
%     dfdy0 = reshape(Y(end,18:306),[17 17]);
%     dfdz0 = dfdy0(9:16,1:8);
%     C=cond(dfdz0);
% else
%     C=1;
% end


%C=1;
%C=norm(x)*norm(grad)/fval;
end

