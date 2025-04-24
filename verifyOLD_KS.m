%проверяем тип решений в KS-переменных
load('Corolev2.mat')
%%
i = 10;
j = 113;
prx = PRevery_a(i,j,1);
pry = PRevery_a(i,j,2);
pvx = PVevery_a(i,j,1);
pvy = PVevery_a(i,j,2);

AN = ANevery_a(i,j)/(2*pi);
T = Tevery_a(i,j)/T_earth_days;
omega0 = Tevery_a(i,j);

r0=[1,0,0];
v0=[0,1,0];
pr_0=[prx,pry,0];
pv_0=[pvx,pvy,0];

y0 = cat(2,r0,v0,pr_0,pv_0,0,0)';
z0 = cat(2,pr_0,pv_0)';
tspan = linspace(t0,t0+T*2*pi, round(T*100));
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
[t,y_fromKS] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);

J_KS = y_fromKS(end,13);       %функционал
AN_KS = y_fromKS(end,14);%угловая дальность

%% решаем первым способом
minimize_delta_omega = @(delta_omega) find_close_solution(delta_omega, omega0, y0, z0, t0,T*2*pi);
angle_range = pi/8;
[delta_omega, Jt_end_min] = fminbnd(minimize_delta_omega, -angle_range, angle_range);
omega_min = delta_omega+omega0;


st.t = [t0, T_earth_days*T];
st.planet = 'Mars';
st.mode = 'Flat';
st.delta_omega = omega_min;

[mars_r_f, mars_v_f]=planetModel(st);
mars_r_f=mars_r_f'*1e+03/r_unit;
mars_v_f=mars_v_f'*1e+03/V_unit;

yf = [mars_r_f;mars_v_f]';

options_fsolve = optimoptions('fsolve','Display','off');
fsolve_traj_fun=@(z)fsolve_traj(z,y0,yf,t0,T*2*pi);
zf = fsolve(fsolve_traj_fun, z0, options_fsolve);
fsolve_traj_equinoctial_fun=@(z)fsolve_traj_equinoctial(z,y0,yf,t0,T*2*pi);
zf2 = fsolve(fsolve_traj_equinoctial_fun, z0, options_fsolve);


pr0=zf(1:3);
pv0=zf(4:6);
y0(7:9)=pr0;
y0(10:12)=pv0;

tspan = linspace(t0,t0+T*2*pi, round(T*2*pi*100));
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   

mug=1;
[~,eMag,iMag,O,~,nu,~,~,lonPer,p] = rv2orb(y0(1:3),y0(4:6),mug);
[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,iMag,O,lonPer,nu);
X = [p;ex;ey;ix;iy;L];
dxdX = equitoctial2decart_jacobian(X,mug);
P = dxdX'*[pr0;pv0];
y0_eq = [X;P];
[~,y_1_eq] = ode113(@(t,y) equitoctial_adjoint_for_integration(y,1), tspan,y0_eq,options);

p_L = y_1_eq(end,end)


[t,y_1] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);
J_1 = y_1(end,13);       %функционал
AN_1 = y_1(end,14);%угловая дальность
%% решаем вторым способом
t_end_0 = T*2*pi;
minimize_delta_t_end = @(delta_t_end) find_close_solution(delta_t_end, t_end_0, y0, z0, t0,AN);
angle_range = pi/8;
[delta_t_end, Jt_end_min] = fminbnd(minimize_delta_t_end, -angle_range, angle_range);
t_end = delta_t_end+t_end_0;


mars_r_f = 1.52*[cos(AN*2*pi);sin(AN*2*pi);0];
mars_v_f = (1/sqrt(1.52))*[cos(AN*2*pi+pi/2);sin(AN*2*pi+pi/2);0];

yf = [mars_r_f;mars_v_f]';

options_fsolve = optimoptions('fsolve','Display','off');
fsolve_traj_fun=@(z)fsolve_traj(z,y0,yf,t0,T*2*pi);
zf = fsolve(fsolve_traj_fun, z0, options_fsolve);

%fsolve_traj_equinoctial_fun=@(z)fsolve_traj_equinoctial(z,y0,yf,t0,T*2*pi);
%zf2 = fsolve(fsolve_traj_equinoctial_fun, z0, options_fsolve);

pr0=zf(1:3);
pv0=zf(4:6);
y0(7:9)=pr0;
y0(10:12)=pv0;

tspan = linspace(t0,t0+T*2*pi, round(T*2*pi*100));
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   

mug=1;
[~,eMag,iMag,O,~,nu,~,~,lonPer,p] = rv2orb(y0(1:3),y0(4:6),mug);
[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,iMag,O,lonPer,nu);
X = [p;ex;ey;ix;iy;L];
dxdX = equitoctial2decart_jacobian(X,mug);
P = dxdX'*[pr0;pv0];
y0_eq = [X;P];
[~,y_2_eq] = ode113(@(t,y) equitoctial_adjoint_for_integration(y,1), tspan,y0_eq,options);

p_L = y_2_eq(end,end)


[t,y_2] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);
J_2 = y_2(end,13);       %функционал
AN_2 = y_2(end,14);%угловая дальность
%% выводим графики
t_orbit = linspace(t0,t0+T_mars/(24*3600), 1000);

st.t = t_orbit';
st.planet = planet_end;
st.mode = orbits;
st.delta_omega = 0;

mars_traj = planetModel(st);
mars_traj=1e+3*mars_traj/ae;

plot(mars_traj(:,1), mars_traj(:,2), 'r', 'DisplayName','орбита Марса')
hold on;
plot(y_fromKS(:,1), y_fromKS(:,2), 'b', 'DisplayName','из KS', 'LineWidth',1.5)
plot(y_fromKS(end,1), y_fromKS(end,2), 'Ob', 'DisplayName','KS точка')
plot(y_1(:,1), y_1(:,2), '--k', 'DisplayName','способ 1', 'LineWidth',1.5)
plot(y_1(end,1), y_1(end,2), 'Ok', 'DisplayName','способ 1 точка')
plot(y_2(:,1), y_2(:,2), '--g', 'DisplayName','способ 2', 'LineWidth',1.5)
plot(y_2(end,1), y_2(end,2), 'Og', 'DisplayName','способ 2 точка')
hold off;
axis equal;
grid;
legend;
xlim([-2,2])
ylim([-2,2])

%% функции
function Jt_end = find_close_solution(delta_omega, omega0, y0, z0, t0,t_end)
    mug_0 = 132712.43994*(10^6)*(10^(3*3));
    ae = 149597870700;
    r_unit=ae;
    V_unit=sqrt(mug_0/ae);
    T_earth_days = 365.256363004;
    omega_min = delta_omega+omega0;
    st.t = [t0, T_earth_days*t_end/(2*pi)];
    st.planet = 'Mars';
    st.mode = 'Flat';
    st.delta_omega = omega_min;
    
    [mars_r_f, mars_v_f]=planetModel(st);
    mars_r_f=mars_r_f'*1e+03/r_unit;
    mars_v_f=mars_v_f'*1e+03/V_unit;
    
    yf = [mars_r_f;mars_v_f]';
    options_fsolve = optimoptions('fsolve','Display','off');
    fsolve_traj_fun=@(z)fsolve_traj(z,y0,yf,t0,t_end);
    zf = fsolve(fsolve_traj_fun, z0, options_fsolve);
    
    %fsolve_traj_equinoctial_fun=@(z)fsolve_traj_equinoctial(z,y0,yf,t0,t_end);
    %zf = fsolve(fsolve_traj_equinoctial_fun, z0, options_fsolve);

    pr0=zf(1:3);
    pv0=zf(4:6);
    y0(7:9)=pr0;
    y0(10:12)=pv0;
    
    tspan = linspace(t0,t0+t_end, round(t_end*100));
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);
    
    
    Jt_end = y(end,13);  
end
function Jt_end = find_close_solution_2(delta_t_end, t_end_0, y0, z0, t0,AN_i)
    
    mars_r_f = 1.52*[cos(AN_i*2*pi);sin(AN_i*2*pi);0];
    mars_v_f = (1/sqrt(1.52))*[cos(AN_i*2*pi+pi/2);sin(AN_i*2*pi+pi/2);0];
    t_end = t_end_0+delta_t_end;
    
    yf = [mars_r_f;mars_v_f]';
    options_fsolve = optimoptions('fsolve','Display','off');
    fsolve_traj_fun=@(z)fsolve_traj(z,y0,yf,t0,t_end);
    zf = fsolve(fsolve_traj_fun, z0, options_fsolve);
    
    %fsolve_traj_equinoctial_fun=@(z)fsolve_traj_equinoctial(z,y0,yf,t0,t_end);
    %zf = fsolve(fsolve_traj_equinoctial_fun, z0, options_fsolve);

    pr0=zf(1:3);
    pv0=zf(4:6);
    y0(7:9)=pr0;
    y0(10:12)=pv0;
    y0(13)=0;
    y0(14)=0;
    tspan = linspace(t0,t0+t_end, round(t_end*100));
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    %options = odeset(options, 'Events',@(t, y) eventIntegrationTrajStopTH(t, y, AN_i*2*pi));
    [t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);
    res = y(end,1:6)-yf(1:6);
    
    Jt_end = y(end,13);       %функционал
    AN_final = y(end,14);%угловая дальность
end
function res = fsolve_traj(z,y0,yf,t0,t_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tspan = linspace(t0,t0+t_end, 10);
pr0=z(1:3);
pv0=z(4:6);
y0(7:9)=pr0;
y0(10:12)=pv0;
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
[~,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);
res = y(end,1:6)-yf(1:6);
end
function res = fsolve_traj_equinoctial(z,y0,yf,t0,t_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tspan = linspace(t0,t0+t_end, 10);
pr0=z(1:3);
pv0=z(4:6);
y0(7:9)=pr0;
y0(10:12)=pv0;
mug=1;
[~,eMag,iMag,O,~,nu,~,~,lonPer,p] = rv2orb(y0(1:3),y0(4:6),mug);
[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,iMag,O,lonPer,nu);
X = [p;ex;ey;ix;iy;L];
dxdX = equitoctial2decart_jacobian(X,mug);
P = dxdX'*[pr0;pv0];
y0 = [X;P];
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
[~,y] = ode113(@(t,y) equitoctial_adjoint_for_integration(y,1), tspan,y0,options);

decart_rv = equitoctial2decart(y(end, 1:6)', mug)';
res = decart_rv-yf(1:6);
end

