% -*- coding: utf-8 -*-
% в этом скрипте изучаем пространственные перелёты с помощью афинных
% преобразований
%clear;
load('mat-files/Preprint2024_1_Mars.mat')
load('mat-files/Preprint2024_1_Venus.mat', 't_end_Venus')
load('mat-files/22-Nov-2024.mat');


ae = 149597870700;
mug_0 = 132712.43994*(10^6)*(10^(3*3));
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;
T_earth_days = 365.256363004;
T_mars_days = 365.256363004*1.8808476;

r_unit=ae;
V_unit=sqrt(mug_0/ae);
T_unit = T_earth/(2*pi);

PVevery_a_MARS_true = PVevery_a(:,113,:); %a = 1.52, Mars
PRevery_a_MARS_true = PRevery_a(:,113,:); %a = 1.52, Mars
Tevery_a_MARS_true = 3600*24*Tevery_a(:,113)/(T_unit*2*pi);
ANevery_a_MARS_true = ANevery_a(:,113)/(2*pi);

PVevery_a_Venus_true = PVevery_a(:,33,:); %a = 1.52, Mars
PRevery_a_Venus_true = PRevery_a(:,33,:); %a = 1.52, Mars
Tevery_a_Venus_true = 3600*24*Tevery_a(:,33)/(T_unit*2*pi);
ANevery_a_Venus_true = ANevery_a(:,33)/(2*pi);
%%    
B=0.2721831;
A=0.5433279;
T_a  = @(a, a_rel)(a_rel-a_rel*B*log(a_rel))*a+(a_rel*B^2*log(a_rel));
d=1;

k_low = 0.72;
k_up = 1.52;

an_low = 0.75;
an_up = 6;

N_points = 100;

k_array_1 = linspace(k_low,k_up,N_points);
k_array_2 = linspace(k_up,k_up,N_points);
k_array_3 = linspace(k_up,k_low,N_points);
k_array_4 = linspace(k_low,k_low,N_points);
k_coef_range = [k_array_1, k_array_2, k_array_3, k_array_4];

an_array_1 = linspace(an_low,an_low,N_points);
an_array_2 = linspace(an_low,an_up,N_points);
an_array_3 = linspace(an_up,an_up,N_points);
an_array_4 = linspace(an_up,an_low,N_points);
an_coef_range = [an_array_1, an_array_2, an_array_3, an_array_4];

N_count_1 = length(k_coef_range);




an_massive_before = zeros(N_count_1,1);
k_massive_before  = zeros(N_count_1,1);

e_massive_before  = zeros(N_count_1,1);
o_massive_before  = zeros(N_count_1,1);

i_massive_before  = zeros(N_count_1,1);
O_massive_before  = zeros(N_count_1,1);

an_massive = zeros(N_count_1,1);
k_massive = zeros(N_count_1,1);

e_massive = zeros(N_count_1,1);
o_massive = zeros(N_count_1,1);

i_massive = zeros(N_count_1,1);
O_massive = zeros(N_count_1,1);

t0 = 0;
start_pos = [1, 0, 0];
start_vel = [0, 1, 0];

for in = 1:N_count_1
    k = k_coef_range(in);
    an = an_coef_range(in);
    start_pos_2 = start_pos;
    start_vel_2 = start_vel;
    if k>1
        j = -1;
    else
        j=1;
    end
    O_target = 1*pi/6;
    %R-преобразование
    R_i = [
    1, 0, 0;
    0, 1/(1+j^2), -j/(1+j^2);
    0, j/(1+j^2), 1/(1+j^2);
    ];

    %восходящий узел
    R_rot = [cos(O_target),-sin(O_target),0;
            sin(O_target),cos(O_target),0;
            0,0,1];
    R_iO=R_rot*R_i*R_rot';
    

    if k==1.52 
        pr_0 = spline(ANevery_a_MARS_true,reshape(PRevery_a_MARS_true, [169, 3])',an);     
        pv_0 = spline(ANevery_a_MARS_true,reshape(PVevery_a_MARS_true, [169, 3])',an);
        t_end = 2*pi*spline(ANevery_a_MARS_true,Tevery_a_MARS_true,an);
    elseif k==0.72
        pr_0 = spline(ANevery_a_Venus_true,reshape(PRevery_a_Venus_true, [169, 3])',an);     
        pv_0 = spline(ANevery_a_Venus_true,reshape(PVevery_a_Venus_true, [169, 3])',an);
        t_end = 2*pi*spline(ANevery_a_Venus_true,Tevery_a_Venus_true,an);
    else
        [pr_0,pv_0] = get_initial_adjoint(an,k,d);
        t_end = 2*pi*(an*(0.03772064567906565*k + 1.405055577366462)*log(k + 1)+(0.08297934358955021*k + 0.009838514316752341).*log(k).^2);
    end
    %t_end = 2*pi*T_a(an, k);

    pr = R_iO^(-1)*pr_0;
    pv = R_iO^(-1)*pv_0;
    y0 = cat(2,start_pos_2,start_vel_2,pr',pv', 0, 0)';
    tspan = linspace(t0,t0+t_end, 10);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);
    
    
    rf = y(end,1:3);
    vf = y(end,4:6);
    [a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(rf',vf',1);
    [p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
    
    J = y(end,13);       %функционал
    AN_final = y(end,14)/(2*pi);%угловая дальность


    an_massive(in) = AN_final;
    k_massive(in) = a;
    
    e_massive(in) = eMag;
    o_massive(in) = lonPer;
    
    i_massive(in) = i;
    O_massive(in) = O;


    pr = pr_0;
    pv = pv_0;
    y0 = cat(2,start_pos_2,start_vel_2,pr',pv', 0, 0)';
    tspan = linspace(t0,t0+t_end, 10);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);
    
    
    rf = y(end,1:3);
    vf = y(end,4:6);
    [a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(rf',vf',1);
    [p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
    
    J = y(end,13);       %функционал
    AN_final = y(end,14)/(2*pi);%угловая дальность


    an_massive_before(in) = AN_final;
    k_massive_before(in) = a;
    
    e_massive_before(in) = eMag;
    o_massive_before(in) = lonPer;
    
    i_massive_before(in) = i;
    O_massive_before(in) = O;
end
%%
figure(1)
clf; % очистить содержимое figure, если нужно
set(gcf, 'Color', 'w'); % белый фон
plot(an_massive_before, k_massive_before, 'Color', 	"#4DBEEE", 'LineWidth', 2.0,'DisplayName', 'До R_i_\Omega преобразования')
hold on;
plot(an_massive, k_massive,'--', 'Color', "#D95319", 'LineWidth', 2.0,'DisplayName', 'После R_i_\Omega преобразования')
hold off;
legend('Location','best');
grid;
xlabel('Угловая дальность, витки')
ylabel('Большая полуось конечной орбиты, безразм.')
set(gca, 'FontSize', 11)

% figure(2)
% clf; % очистить содержимое figure, если нужно
% set(gcf, 'Color', 'w'); % белый фон
% scatter(o_massive_before, e_massive_before, 'o', 'Color', 	"#4DBEEE", 'LineWidth', 2.0,'DisplayName', 'до преобразования')
% hold on;
% scatter(o_massive, e_massive, '+', 'Color', "#D95319", 'LineWidth', 2.0,'DisplayName', 'после преобразования')
% hold off;
% legend('Location','best');
% grid;
% xlabel('Аргумент перицентра')
% ylabel('Эксцентриситет')
% set(gca, 'FontSize', 11)

figure(3)
clf; % очистить содержимое figure, если нужно
set(gcf, 'Color', 'w'); % белый фон
plot(O_massive_before, i_massive_before, 'Color', 	"#4DBEEE", 'LineWidth', 2.0,'DisplayName', 'до преобразования')
hold on;
plot(O_massive, i_massive, 'Color', "#D95319", 'LineWidth', 2.0,'DisplayName', 'после преобразования')
hold off;
legend('Location','best');
grid;
xlabel('Восходящий узел')
ylabel('Наклонение')
set(gca, 'FontSize', 11)
function res = fsolve_find_i(i_matrix,i_target, start_pos,start_vel,zf,t0,t_end,O_target)
%изменение наклонения
delta_i = i_matrix;
n_rot_i = [1; 0; 0];
q_rot_i = [cos(delta_i/2); n_rot_i*sin(delta_i/2)];
R_i0 = quat2rotm(q_rot_i');
R_i0(1:3,3)=R_i0(1:3,3)*R_i0(2,2);
R_i0(1:3,2)=R_i0(1:3,2)*R_i0(2,2);
%восходящий узел
%O_target = 0*pi/6;
delta_O = O_target-pi;
n_rot = [0; 0; 1];
q_rot = [cos(delta_O/2); n_rot*sin(delta_O/2)];
R_rot = quat2rotm(q_rot');
R_iO=R_rot*R_i0*R_rot';
[i,O,J,AN] = integrate_iO(R_iO, start_pos,start_vel,zf,t0,t_end);
res = i-i_target;
end
function [i,O,p,a,J,AN] = integrate_iO(R_iO, start_pos,start_vel,zf,t0,t_end)
pr = R_iO^(-1)*zf(1:3);
pv = R_iO^(-1)*zf(4:6);
y0 = cat(2,start_pos,start_vel,pr',pv',0,0)';
tspan = linspace(t0,t0+t_end, 1000);
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
[t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);

%J = y(end,13);       %функционал
%AN = y(end,14);%угловая дальность
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(y(end, 1:3)',y(end, 4:6)',1);
[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
end
function [ex,ey,lonPer,p,a] = integrate_iO_equinoctial(R_iO, start_pos,start_vel,zf,t0,t_end)
pr = R_iO^(-1)*zf(1:3);
pv = R_iO^(-1)*zf(4:6);

load('mat-files/Preprint2024_1_Mars.mat','dxdX');
P = dxdX'*[pr;pv];
y0 = [1;0;0;0;0;0;P];
tspan = linspace(t0,t0+t_end, 1000);
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
mug = 1;
[t,y] = ode113(@(t,y) equitoctial_adjoint_for_integration(y,mug), tspan,y0,options);

traj = arrayfun(@(p,ex,ey,ix,iy,L) equitoctial2decart([p;ex;ey;ix;iy;L], mug),...
    y(:, 1),y(:, 2),y(:, 3),y(:, 4),y(:, 5),y(:, 6),'UniformOutput',false);
traj = cell2mat(traj')';

[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(traj(end, 1:3)',traj(end, 4:6)',1);
[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
end
function t_end = get_t_end(k)
load('mat-files/Preprint2024_1_universal_2.mat','k_coef_for_spline','T_40_array')
t_end = spline(k_coef_for_spline,T_40_array',k);
end
function [pr_0,pv_0] = get_adjoint(k)
load('mat-files/Preprint2024_1_universal_2.mat','k_coef_for_spline','P_40_array')
z = spline(k_coef_for_spline,P_40_array',k);
pr_0 = z(1:3);
pv_0 = z(4:6);
end
function out1 = Ca0fromIK(i, k)
  out1 = k + i.*(i.*k.*(i.*(0.75708324 - 0.22989681*k).*exp(k) + 1.637357) + 0.031784486);
end
function out1 = Ca1fromIK(i, k)
  out1 = i.^2.*(k + 2.518653).*(0.050697707*k.*(-6.1506586*i + k) - 0.8663172);
end
function out1 = Ci0fromK(x0)

  out1 = 0.415541771668952*x0.*exp(-1.0*x0).*log(x0).^2 + 0.25060964*log(x0) + 0.0796329*exp(-5.9961753*x0).*log(x0).^2;

end
function out1 = Ci1fromK(x0)

  out1 = 0.16930881*x0.*(x0 - 1.1203467).*(x0 - 0.8851997).*exp(x0).*log(x0)./(x0.^2.*(x0 - 0.8851997).*exp(x0) + 1.3398796*exp(x0) + log(2*x0));

end
