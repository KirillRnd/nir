%% В этом скрипте получаем ршение задачи встречи для 40 витков
% делал для препринта весной 2024
%внешние константы
%clc;
clear;
ae = 149597870700;
mug_0 = 132712.43994*(10^6)*(10^(3*3));
mug=1;
T_earth_days = 365.256363004;
T_earth = T_earth_days*3600*24;
T_mars_days = T_earth_days*1.8808476;
T_mars=T_earth*1.8808476;

r_unit=ae;
V_unit=sqrt(mug_0/ae);
T_unit = T_earth/(2*pi);
%гиперпараметры
t_start=0;
orbits='Flat';
N=1350;
m0=367;
eta=0.45;
planet_start = 'Earth';
planet_end = 'Mars';
omega = -pi;

%определяем стартовое положение
st.t = t_start;
st.planet = planet_start;
st.mode = orbits;
st.delta_omega = omega;

[start_pos, start_vel] = planetModel(st);
start_pos=start_pos*1e+03/ae;
start_vel=start_vel*1e+03/V_unit;


%%
k_coef = 0.1:0.01:10;
N_count_1 = length(k_coef);
T_40_array = zeros(N_count_1, 1);
P_40_array = zeros(N_count_1, 6);
T_40_array(91)=2*pi*40;
savefilename = 'mat-files/Preprint2024_1_universal_2.mat';
load(savefilename)
for j1 = 396:1:400%171
    
    a_rel = k_coef(j1);
    if a_rel == 1
        continue
    end
    disp(a_rel)
    %интегрируем
    param = struct();
    param.p_scaler = 1e-1;
    param.t_scaler = 1e+1;
    t0=0;
    d_coef = 0;
    AN_i = 40;
    [pr_0,pv_0] = get_initial_adjoint(AN_i,a_rel,d_coef);
    pr_0 = (pr_0)';
    pv_0 = (pv_0)';
    %z0 = [pr_0,pv_0];
    delta_z0 = P_40_array(j1-1,:)-P_40_array(j1-2,:);
    z0 = P_40_array(j1-1,:)+delta_z0;
    z0 = P_40_array(j1,:);
    [a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(start_pos',start_vel',mug);
        
    [p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
    X0 = [p;ex;ey;ix;iy;L];
    dxdX = equitoctial2decart_jacobian(X0,mug);
    Z0 = dxdX'*z0';
    Y0 = [X0;Z0];
    %t_end_0 = 2*pi*TfromAN(AN_i, a_rel);
    delta_t_end_0 = T_40_array(j1-1)-T_40_array(j1-2);
    t_end_0 = T_40_array(j1-1)+delta_t_end_0;
    t_end_0 = T_40_array(j1);
    Z0_withT = [Z0/param.p_scaler;0];
    
    
    % mars_r_f = a_rel*[cos(AN_i*2*pi);sin(AN_i*2*pi);0];
    % mars_v_f = (1/sqrt(a_rel))*[cos(AN_i*2*pi+pi/2);sin(AN_i*2*pi+pi/2);0];
    % yf = [mars_r_f;mars_v_f]';
    yf = [a_rel,0,0,0,0,AN_i*2*pi];
    
    param.y0 = Y0;
    param.yf = yf;
    param.t0 = t0;
    param.t_end_0 = t_end_0;
    fun2min_border_ex = @(z)fun2min_border(z,param);
    fun2min_J_ex = @(z)fun2min_J(z,param);
    options = optimoptions('fmincon','Display','iter');
    options = optimoptions(options,'UseParallel', true);
    options = optimoptions(options, 'OptimalityTolerance', 1e-10);
    options = optimoptions(options, 'MaxFunctionEvaluations', 1e+10);
    options = optimoptions(options, 'StepTolerance', 1e-10);
    options = optimoptions(options, 'ConstraintTolerance', 1e-10);
    options = optimoptions(options, 'MaxIterations', 250);
    %options = optimoptions(options, 'FiniteDifferenceType', 'central');
    %options = optimoptions(options, 'EnableFeasibilityMode', true);
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = -1e-0*ones(1,7);
    ub = 1e-0*ones(1,7);
    angle_range = 16*pi/param.t_scaler;
    lb(7)=-angle_range;
    ub(7)=+angle_range;
    ZF_withT = fmincon(fun2min_J_ex, Z0_withT, A, b, Aeq, beq, lb, ub, fun2min_border_ex, options);
    zf = inv(dxdX)'*ZF_withT(1:6)*param.p_scaler;
    t_end_delta = ZF_withT(7)*param.t_scaler;
    t_end = t_end_0+t_end_delta;
    T_40_array(j1) = t_end;
    P_40_array(j1,:) = zf;
    j1_last=j1;
    save(savefilename);

end
%%
j1=245;
param = struct();
param.p_scaler = 1e-1;
param.t_scaler = 1e+1;
t0=0;
d_coef = 0;
a_rel = k_coef(j1);
T_orb = 2*pi*a_rel^(3/2);
AN_i = 40;
[pr_0,pv_0] = get_initial_adjoint(AN_i,a_rel,d_coef);
pr_0 = (pr_0)';
pv_0 = (pv_0)';
z0 = [pr_0,pv_0];
z0 = P_40_array(j1,:);
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(start_pos',start_vel',mug);

[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
X0 = [p;ex;ey;ix;iy;L];
dxdX = equitoctial2decart_jacobian(X0,mug);
Z0 = dxdX'*z0';
Y0 = [X0;Z0];
%t_end_0 = 2*pi*TfromAN(AN_i, a_rel);
%t_end_0 = 2*pi*TfromK(a_rel);
t_end_0 = T_40_array(j1);
Z0_withT = [Z0/param.p_scaler;0];


mars_r_f = a_rel*[cos(AN_i*2*pi);sin(AN_i*2*pi);0];
mars_v_f = (1/sqrt(a_rel))*[cos(AN_i*2*pi+pi/2);sin(AN_i*2*pi+pi/2);0];
%yf = [mars_r_f;mars_v_f]';
yf = [a_rel,0,0,0,0,AN_i*2*pi];
param.y0 = Y0;
param.yf = yf;
param.t0 = t0;
param.t_end_0 = t_end_0;
fun2min_border_ex = @(z)fun2min_border(z,param);
fun2min_J_ex = @(z)fun2min_J(z,param);
options = optimoptions('fmincon','Display','iter');
options = optimoptions(options,'UseParallel', true);
options = optimoptions(options, 'OptimalityTolerance', 1e-10);
options = optimoptions(options, 'MaxFunctionEvaluations', 1e+10);
options = optimoptions(options, 'StepTolerance', 1e-10);
options = optimoptions(options, 'ConstraintTolerance', 1e-10);
options = optimoptions(options, 'MaxIterations', 250);
%options = optimoptions(options, 'FiniteDifferenceType', 'central');
%options = optimoptions(options, 'EnableFeasibilityMode', true);
A = [];
b = [];
Aeq = [];
beq = [];
lb = -1e-0*ones(1,7);
ub = 1e-0*ones(1,7);
angle_range = 8*pi/param.t_scaler;
%lb(4:5)=0;
%ub(4:5)=0;
lb(7)=-angle_range;
ub(7)=+angle_range;
ZF_withT = fmincon(fun2min_J_ex, Z0_withT, A, b, Aeq, beq, lb, ub, fun2min_border_ex, options);
zf = inv(dxdX)'*ZF_withT(1:6)*param.p_scaler;
t_end_delta = ZF_withT(7)*param.t_scaler;
t_end = t_end_0+t_end_delta;
% %T_40_array(j1) = t_end;
% P_40_array(j1,:) = zf;
%%
figure(1)
plot(k_coef,T_40_array)
grid()

figure(2)
plot(k_coef,P_40_array(:,1))
hold on;
plot(k_coef,P_40_array(:,5))
hold off
grid()
figure(3)
plot(k_coef,P_40_array(:,2))
hold on;
plot(k_coef,P_40_array(:,4))
hold off
grid()
%%
%j1 = 990;
%AN_i=40;
a_rel = 2.5;
%t_end = 2*pi*TfromK(a_rel);
%t_end = 2*pi*TfromAN(AN_i, a_rel);
%t_end = T_40_array(j1);
%zf = PfromK(a_rel)';
%zf = P_40_array(j1,:)';
zf = spline(k_coef_for_spline,P_40_array',10);
t_end = spline(k_coef_for_spline,T_40_array',a_rel);
%zf=z0';
%t_end = t_end_0;

y0 = cat(2,start_pos,start_vel,zf(1:3)',zf(4:6)',0,0)';
tspan = linspace(t0,t0+t_end, AN_i*400);
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
[t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);

J = y(end,13);       %функционал
AN_final = y(end,14);%угловая дальность
res_r_f = norm(y(end, 1:3)-mars_r_f');
res_v_f = norm(y(end, 4:6)-mars_v_f');
%Z0_Mars = ZF;
t_end_Mars = t_end;

%рисуем
t0 = t_start;
t_orbit = linspace(t0,t0+T_earth/(24*3600), 1000);

st.t = t_orbit';
st.planet = 'Earth';
st.mode = orbits;
st.delta_omega = omega;

earth_traj = planetModel(st);
earth_traj=1e+3*earth_traj/ae;

t_orbit = linspace(t0,t0+T_mars/(24*3600), 1000);

% st.t = t_orbit';
% st.planet = planet_end;
% st.mode = orbits;
% st.delta_omega = 0;
% 
% mars_traj = planetModel(st);
% mars_traj=1e+3*mars_traj/ae;
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
y0 = cat(2,y(1, 1:3),y(1, 4:6),[0,0,0],[0,0,0])';
[t,orbit_initial] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
y0 = cat(2,y(end, 1:3),y(end, 4:6),[0,0,0],[0,0,0])';
[t,orbit_final] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);


%выводим график
figure(5);
plot3(0, 0, 0, 'k--o');
set(gca,'FontSize',14);
hold on;

plot3(y(1, 1), y(1, 2), y(1, 3), 'O', 'LineWidth', 1);
plot3(y(:, 1), y(:, 2), y(:, 3), 'cyan', 'LineWidth', 1);
plot3(orbit_initial(:, 1), orbit_initial(:, 2), orbit_initial(:, 3), 'k');
plot3(orbit_final(:, 1), orbit_final(:, 2), orbit_final(:, 3), 'r');
%plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k')
%plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'k')
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(y(end, 1:3)',y(end, 4:6)',1);
disp(['e=', num2str(eMag)])
disp(['a=', num2str(a)])
disp(['J=', num2str(J)])
%disp(['i=', num2str(i)])
plot3(y(end, 1), y(end, 2), y(end, 3),'bO')

%plot3(mars_r_f(1), mars_r_f(2), mars_r_f(3),'rO')

hold off;
axis equal

%title('Траектория КА')
xlabel('x, AU')
ylabel('y, AU')
zlabel('z, AU')
view(0,90)
%view(90,0)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
grid on;
box off;
%%
function Jt_end = find_close_solution(delta_t_end, t_end_0, y0, Z0, t0,AN_i,kappa)
    
    mars_r_f = kappa*[cos(AN_i*2*pi);sin(AN_i*2*pi);0];
    mars_v_f = (1/sqrt(kappa))*[cos(AN_i*2*pi+pi/2);sin(AN_i*2*pi+pi/2);0];
    t_end = t_end_0+delta_t_end;
    start_pos = y0(1:3);
    start_vel = y0(4:6);
    mug = 1;
    [~,eMag,i,O,~,nu,~,~,lonPer,p] = rv2orb(start_pos',start_vel',mug);

    [p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
    X0 = [p;ex;ey;ix;iy;L];
    dxdX = equitoctial2decart_jacobian(X0,mug);
    Y0 = [X0;Z0];

    yf = [mars_r_f;mars_v_f]';
    % options_fsolve = optimoptions('fsolve','Display','off');
    % fsolve_traj_fun=@(z)fsolve_traj(z,Y0,yf,t0,t_end);
    % ZF = fsolve(fsolve_traj_fun, Z0, options_fsolve);
    p_scaler=1e-0;
    fmincon_fun=@(z)fmincon_traj(z,Y0,yf,t0,t_end,p_scaler);
    options = optimoptions('fmincon','Display','iter');
    options = optimoptions(options,'UseParallel', true);
    options = optimoptions(options, 'OptimalityTolerance', 1e-10);
    options = optimoptions(options, 'MaxFunctionEvaluations', 1e+10);
    options = optimoptions(options, 'StepTolerance', 1e-10);
    options = optimoptions(options, 'ConstraintTolerance', 1e-10);
    options = optimoptions(options, 'MaxIterations', 1000);
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = -1e-0*ones(1,6);
    ub = 1e-0*ones(1,6);
    ZF = fmincon(@(x)1, Z0/p_scaler, A, b, Aeq, beq, lb, ub, fmincon_fun, options);
    Y0(7:12)=ZF*p_scaler;
    Y0(13)=0;
    % zf = inv(dxdX)'*ZF;
    % 
    % pr0=zf(1:3);
    % pv0=zf(4:6);
    % y0(7:9)=pr0;
    % y0(10:12)=pv0;
    % y0(13:14)=0;
    tspan = linspace(t0,t0+t_end, round(t_end*10));
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    % [t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);
    % 
    % Jt_end = y(end, 13);
    [t,y] = ode113(@(t,y) equitoctial_adjoint_for_integration_withJ(y,mug), tspan,Y0,options);

    traj = arrayfun(@(p,ex,ey,ix,iy,L) equitoctial2decart([p;ex;ey;ix;iy;L], mug),...
    y(:, 1),y(:, 2),y(:, 3),y(:, 4),y(:, 5),y(:, 6),'UniformOutput',false);
    traj = cell2mat(traj')';
    res = traj(end,1:6)-yf(1:6);
    if kappa>1
        Jt_end = -y(end, 13);
    else
        Jt_end = y(end, 13);
    end
end
function [c, ceq] = fmincon_traj(z,y0,yf,t0,t_end, p_scaler)
    c =[];
    ceq = fsolve_traj(z*p_scaler,y0,yf,t0,t_end);
end
function res = fsolve_traj(z,y0,yf,t0,t_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
y0(7:12)=z;
tspan = linspace(t0,t0+t_end, 10);
mug = 1;
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   

[t,y] = ode113(@(t,y) equitoctial_adjoint_for_integration(y,mug), tspan,y0,options);

traj = arrayfun(@(p,ex,ey,ix,iy,L) equitoctial2decart([p;ex;ey;ix;iy;L], mug),...
    y(:, 1),y(:, 2),y(:, 3),y(:, 4),y(:, 5),y(:, 6),'UniformOutput',false);
traj = cell2mat(traj')';
res = (traj(end,1:6)-yf(1:6));
end
function [c, ceq] = fun2min_border(z,param)
y0 = param.y0;
yf = param.yf;
t0 = param.t0;
t_end_0 = param.t_end_0;
t_end = t_end_0+z(7)*param.t_scaler;
y0(7:12)=param.p_scaler*z(1:6);
tspan = linspace(t0,t0+t_end, 1000);
mug = 1;
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   

[~,y] = ode113(@(t,y) equitoctial_adjoint_for_integration(y,mug), tspan,y0,options);

% traj = arrayfun(@(p,ex,ey,ix,iy,L) equitoctial2decart([p;ex;ey;ix;iy;L], mug),...
%     y(:, 1),y(:, 2),y(:, 3),y(:, 4),y(:, 5),y(:, 6),'UniformOutput',false);
% traj = cell2mat(traj')';
% c =[];
% ceq = (traj(end,1:6)-yf(1:6));
c =[];
ceq = (y(end,1:6)-yf(1:6));
end
function J  = fun2min_J(z,param)
y0 = param.y0;
t0 = param.t0;
t_end_0 = param.t_end_0;
t_end = t_end_0+z(7)*param.t_scaler;
y0(7:12)=param.p_scaler*z(1:6);
y0(13)=0;
tspan = linspace(t0,t0+t_end, 10);
mug = 1;
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   

[~,y] = ode113(@(t,y) equitoctial_adjoint_for_integration_withJ(y,mug), tspan,y0,options);

J = y(end,13);
end
function out1 = TfromK(k)
  %TFROMK  Autogenerated by SymPy
  %   Code generated with SymPy 1.12
  %
  %   See http://www.sympy.org/ for more information.
  %
  %   This file is part of 'project'

  out1 = 26.116552*k - exp(k) + 5.8698006*log(k) + 16.684269;

end
function out1 = PfromK(k)
  %TFROMK  Autogenerated by SymPy
  %   Code generated with SymPy 1.12
  %
  %   See http://www.sympy.org/ for more information.
  %
  %   This file is part of 'project'

  out1 = [pr0fromK(k),pr1fromK(k),0,pv0fromK(k),pv1fromK(k),0];

end
function out1 = pr0fromK(k)
  %PR0FROMK  Autogenerated by SymPy
  %   Code generated with SymPy 1.12
  %
  %   See http://www.sympy.org/ for more information.
  %
  %   This file is part of 'project'

  out1 = (-log(k) + log(1.000694*k + 0.0012779579877598)).*log(k + 0.008227128);

end
function out1 = pr1fromK(k)
  %PR1FROMK  Autogenerated by SymPy
  %   Code generated with SymPy 1.12
  %
  %   See http://www.sympy.org/ for more information.
  %
  %   This file is part of 'project'

  out1 = 6.793992e-8*k.*exp(1.0*k) - 2.24543867849136e-6*k.*log(k) - 0.0056721005*k.*exp(-10.655896*k);

end
function out1 = pv0fromK(k)
  %PV0FROMK  Autogenerated by SymPy
  %   Code generated with SymPy 1.12
  %
  %   See http://www.sympy.org/ for more information.
  %
  %   This file is part of 'project'

  out1 = -1.9245615e-6*k + 2.19326050101555e-6 - 0.00079736474*exp(-8.374843*k);

end
function out1 = pv1fromK(k)
  %PV1FROMK  Autogenerated by SymPy
  %   Code generated with SymPy 1.12
  %
  %   See http://www.sympy.org/ for more information.
  %
  %   This file is part of 'project'

  out1 = (-log(k) + log(k + 0.0012918527) + 0.000685634993962657).*log(k + 0.008998305);

end