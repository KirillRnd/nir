function checkMethod_results = checkMethod_ver2(checkMethod_params_INPUT)
%UNTITLED9 Summary of this function goes here
%   Вычисляет невязку в зависимости от входных параметров
%условия на fmincon
%ЗАДАЧА ПРОЛЁТА case_traj=1; ЗАДАЧА сопровождения case_traj=2;

%Начальные условия
%x0=zeros([1, 11]);
checkMethod_params = checkMethod_params_INPUT;
t_start = checkMethod_params.t_start;
psi = checkMethod_params.delta_s;
rad = checkMethod_params.rad;
UorR = checkMethod_params.UorR;
decreaseNonPsysical = checkMethod_params.decreaseNonPsysical;
modifier_p = checkMethod_params.modifier_p;
modifier_f = checkMethod_params.modifier_f;
x0 = checkMethod_params.x0_sec;
eta = checkMethod_params.eta;
case_traj = checkMethod_params.case_traj;
planet_end = checkMethod_params.planet_end;
display = checkMethod_params.display;
terminal_state = checkMethod_params.terminal_state;
integration_acc = checkMethod_params.integration_acc;
calculate_condition = checkMethod_params.calculate_condition;
orbits = checkMethod_params.orbits;
omega = checkMethod_params.omega;

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
st.a_rel = checkMethod_params.a_rel;

[r0, V0] = planetModel(st);

eul = [pi/12 pi/4 pi/12];
rotmZYX = eul2rotm(eul);


if isfield(checkMethod_params,'dV_add')
    V0 = V0+checkMethod_params.dV_add*1e-03; %добавляем гип. избыток в км/с
elseif isfield(checkMethod_params,'dV_value')
    hi = checkMethod_params.hi;
    %dV_add = checkMethod_params.dV_value*[cos(hi), sin(hi), 0]*1e-03;
    st_f2m.hi = hi;
%     st_f2m.dV_add = [dV_add'/V_unit]*1e+03;
    st_f2m.dV_value = checkMethod_params.dV_value/V_unit;
    st_f2m.Vp = [rotmZYX*V0'/V_unit]*1e+03;
    st_f2m.rotmZYX = rotmZYX;
    %V0 = V0+dV_add; %добавляем гип. избыток в км/с
elseif isfield(checkMethod_params,'dV_value_end')
    hi = checkMethod_params.hi;
    %dV_add = checkMethod_params.dV_value*[cos(hi), sin(hi), 0]*1e-03;
    st_f2m.hi = hi;
%     st_f2m.dV_add = [dV_add'/V_unit]*1e+03;
    st_f2m.dV_value = checkMethod_params.dV_value_end/V_unit;
    st_f2m.Vp = [rotmZYX*V0'/V_unit]*1e+03;
    st_f2m.rotmZYX = rotmZYX;
    %V0 = V0+dV_add; %добавляем гип. избыток в км/с
end

r0 = [rotmZYX*r0'/ae; 0]*1e+03;
V0 = [rotmZYX*V0'/V_unit; 0]*1e+03;

%V0 = V0 + 2700*V0/norm(V0)/V_unit;%гип избыток
%nonlcon = @(x)ubOrtPv(x, rToU(r0, 0));
%nonlcon=[];
mug=1;
%x0(1:8)=x0(1:8)/modifier_p;
% modifier_p=1e-04;
% modifier_f=1e+04;
modifier_b=1e-1;

%s_a = psi-rad;
s_a = psi;
s_b = psi+rad*2*pi;

%x0(11)=psi;
box_initial_point = [x0(1:8), 0, 0];
box_restrictions = [1, 1, 1, 1, 1, 1, 1, 1, 0, 0]*modifier_b;
lb = box_initial_point-box_restrictions;
ub = box_initial_point+box_restrictions;

lb(9) = s_a/(2*pi);
ub(9) = s_b/(2*pi);


if strcmp(UorR,'u') || strcmp(UorR,'r')
    lb(10) = -1.0;
    ub(10) = 1.0;
elseif strcmp(UorR,'u_hat')
    lb(10) = x0(10);
    ub(10) = x0(10);
    lb(11) = 0;
    ub(11) = 0;
    x0(11) = 0;
%     lb(10) = 0.0;
%     ub(10) = 1.0;
elseif strcmp(UorR,'u_hat_v_hyp')
    lb(10) = x0(10);
    ub(10) = x0(10);
     lb(11) = -pi;
     ub(11) = pi;
     x0(11) = pi/64;
%     lb(10) = 0.0;
%     ub(10) = 1.0;
elseif strcmp(UorR,'u_hat_v_hyp_end')
    lb(10) = x0(10);
    ub(10) = x0(10);
     lb(11) = -pi;
     ub(11) = pi;
     x0(11) = pi/64;
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
st_f2m.a_rel = checkMethod_params.a_rel;
tic;
fun=@(x)fun2min([x(1:8), x(9), x(10), x(11)], st_f2m);

options = optimoptions('fmincon','UseParallel', true);
if display == 1
    options = optimoptions(options, 'Display', 'iter');
else
    options = optimoptions(options, 'Display', 'off');
end
options = optimoptions(options, 'OptimalityTolerance', 1e-10);
options = optimoptions(options, 'MaxFunctionEvaluations', 1e+10);
options = optimoptions(options, 'StepTolerance', 1e-10);
options = optimoptions(options, 'ConstraintTolerance', 1e-10);
options = optimoptions(options, 'MaxIterations', 250);
%options = optimoptions(options, 'FiniteDifferenceType', 'central');
%options = optimoptions(options, 'EnableFeasibilityMode', true);
%options = optimoptions(options, 'Algorithm', 'sqp');

options = optimoptions(options,'OutputFcn',@myoutput);

[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x)1, x0, A, b, Aeq, beq, lb, ub, fun, options);
evaluation_time = toc;
px = x(1:8);
if terminal_state == 's'
    s_f=x(9)*2*pi;
elseif terminal_state == 't'
    s_f=1.5*x(9)*2*pi;
end
t_end_0=x(9)*365.256363004;
phi = x(10)*2*pi;
delta_hi = x(11);
%задаем начальные условия
%options = optimoptions(options,'OutputFcn',@myoutput);
%options = optimoptions(options, 'Algorithm', 'sqp');
phi0=phi;
%phi0=0;
if strcmp(UorR,'u_hat_v_hyp')
    hi = st_f2m.hi + delta_hi;
    dV_add = st_f2m.dV_value*st_f2m.rotmZYX*[cos(hi), sin(hi), 0]';
    V0 = st_f2m.V0+[dV_add; 0];
end
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
        tr=(px_end'*ortdgduv)';
        p_plus=[g;tr];

        [s,y] = ode113(@(s,y) integrateTraectory(s,y),int_s0sf,y0-y0_delta, options);
        u_end=y(end,1:4)';
        w_end=y(end,5:8)';
        px_end=y(end,9:16)';
        g=get_target_g(u_end,w_end);
        ortdgduv=get_ortdgduv(u_end,w_end);
        tr=(px_end'*ortdgduv)';
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
    p_tr=([pu;pw]'*f_ortdgduv)';
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
st.a_rel = checkMethod_params.a_rel;

[mars_r_f, mars_v_f]=planetModel(st);
mars_r_f=rotmZYX*mars_r_f'*1e+03;
mars_v_f=rotmZYX*mars_v_f'*1e+03;

if strcmp(st_f2m.UorR,'u_hat_v_hyp_end')
     hi = st_f2m.hi + delta_hi;
     dV_add = st_f2m.dV_value*st_f2m.rotmZYX*[cos(hi), sin(hi), 0]';
     mars_v_f = mars_v_f+dV_add;
end

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
checkMethod_results = struct();
checkMethod_results.dr = dr;
checkMethod_results.dv = dv;
checkMethod_results.C = C;
checkMethod_results.px = px;
checkMethod_results.s_f = s_f;
checkMethod_results.phi = phi;
checkMethod_results.t_end = t_end;
checkMethod_results.s = s;
checkMethod_results.uu = uu;
checkMethod_results.rr = rr;
checkMethod_results.VV = VV;
checkMethod_results.t = t;
checkMethod_results.Jt = Jt;
checkMethod_results.a_ks = a_ks;
checkMethod_results.evaluation_time = evaluation_time;
checkMethod_results.hi_opt = st_f2m.hi+delta_hi;
end

