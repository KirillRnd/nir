% в этом скрипте изучаем пространственные перелёты с помощью афинных
% преобразований
clear;
load('mat-files/Preprint2024_1_Mars.mat')
load('mat-files/Preprint2024_1_Venus.mat', 't_end_Venus')
%%
% k_coef_range=0.4:0.1:2.5;
% i_coef_range=0:pi/120:pi/6;
% O_coef_range=0:pi/15:2*pi;
% N_count_1 = length(k_coef_range);
% N_count_2 = length(i_coef_range);
% N_count_3 = length(O_coef_range);
% 
% T_right = zeros(N_count_1,1);
% a_massive = zeros(N_count_1,N_count_2,N_count_3);
% a_massive_error = zeros(N_count_1,N_count_2,N_count_3);
% i_massive = zeros(N_count_1,N_count_2,N_count_3);
% i_massive_error = zeros(N_count_1,N_count_2,N_count_3);
% for j1 = 1:N_count_1
%     for j2 = 1:1
%         for j3 = 1:1
%             disp([num2str(j1),'   ',num2str(j2),'   ',num2str(j3)]);
%             k_coef = k_coef_range(j1);
%             i_coef = i_coef_range(j2);
%             O_coef = O_coef_range(j3);
%             t_end_0 = get_t_end(k_coef);
%             [pr_0,pv_0] = get_adjoint(k_coef);
%             z = [pr_0;pv_0];
%             R_iO = eye(3);
%             AN_target = 40*2*pi;
% 
%             options_fsolve = optimoptions('fsolve','Display','off');
%             fsolve_fun=@(t_end)fsolve_find_t_end(AN_target, start_pos,start_vel,z,t0,t_end);
%             t_end = fsolve(fsolve_fun, t_end_0, options_fsolve);
%             T_right(j1) = t_end;
%             a_massive(j1,j2,j3) = a;
%         end
%     end
% end
%%
savefilename = 'mat-files/Preprint2024_5.mat';
load(savefilename)
for j1 = 1:N_count_1
    for j2 = 1:1%N_count_2
        for j3 = 1:1%N_count_3
        k_coef = k_coef_range(j1);
        i_coef = i_coef_range(j2);
        O_coef = O_coef_range(j3);
        
        if a_massive(j1,j2,j3) > 0
            continue
        end
        disp([num2str(j1),'   ',num2str(j2),'   ',num2str(j3)]);
        %изменение наклонения
        a_rel = k_coef;
        Ci0 = (0.38159616342652647-0.0243395353632856*a_rel-0.2767778010385592*exp(-a_rel^2))*log(a_rel);
        Ci1 = -a_rel*exp(-7.8328443*sinh(a_rel))+5.5964956*exp(-sinh((1.38935843796911*a_rel^2+1.09321695795059)/(0.7103882*a_rel^2+0.04272333)));
        delta_i_coef = Ci0+Ci1*cos(2*O_coef);
        if delta_i_coef==0
            delta_i_coef = 1;
        end
        delta_i_0 = i_coef/delta_i_coef;
        t_end = T_right(j1);
        x0 = [delta_i_0,a_rel];
        try
            fmincon_fun=@(x)fmincon_find_a_i(x, i_coef,k_coef, O_coef, start_pos,start_vel,t0,t_end);
            options = optimoptions('fmincon','Display','off');
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
            lb = [0, 0.2];
            ub = [30, 4];
            x = fmincon(@(x)1, x0, A, b, Aeq, beq, lb, ub, fmincon_fun, options);
            i_input = x(1);
            a_input = x(2);
            i_massive(j1,j2,j3) = i_input;
            a_massive(j1,j2,j3) = a_input;

            %наклонение
            delta_i = atan(i_input);
            n_rot_i = [1; 0; 0];
            q_rot_i = [cos(delta_i/2); n_rot_i*sin(delta_i/2)];
            R_i0 = quat2rotm(q_rot_i');
            R_i0(1:3,3)=R_i0(1:3,3)*R_i0(2,2);
            R_i0(1:3,2)=R_i0(1:3,2)*R_i0(2,2);
            %восходящий узел
            delta_O = O_coef-pi;
            n_rot = [0; 0; 1];
            q_rot = [cos(delta_O/2); n_rot*sin(delta_O/2)];
            R_rot = quat2rotm(q_rot');
            R_iO=R_rot*R_i0*R_rot';
            
            [pr_0,pv_0] = get_adjoint(a_input);
            z = [pr_0;pv_0];
            
            [i,O,p,a,J,AN] = integrate_iO(R_iO, start_pos,start_vel,z,t0,t_end);
            i_massive_error(j1,j2,j3) = i-i_coef;
            a_massive_error(j1,j2,j3) = a-k_coef;
        catch

            i_massive(j1,j2,j3) = -1000;
            a_massive(j1,j2,j3) = -1000;
        end
        %save(savefilename)
        end
    end
end
%% исправляем неправильный масштаб k
k_coef_range_fix=zeros(N_count_1,1);
for j1 = 1:N_count_1
    t_end = T_right(j1);
        
    k_coef = k_coef_range(j1);

    R_iO = eye(3);
    [pr_0,pv_0] = get_adjoint(k_coef);
    z = [pr_0;pv_0];
    
    [i,O,p,a,J,AN] = integrate_iO(R_iO, start_pos,start_vel,z,t0,t_end);
    k_coef_range_fix(j1) = a;
end
%%
%k = 2.5;
j1 = 21;
j2 = 1;
j3 = 1;
k_coef = k_coef_range(j1)
i_coef = i_coef_range(j2);
O_coef = O_coef_range(j3);

a_input = a_massive(j1,j2,j3)
a_input = k_coef;
i_input = i_massive(j1,j2,j3);

delta_i = atan(-i_input);
n_rot_i = [1; 0; 0];
q_rot_i = [cos(delta_i/2); n_rot_i*sin(delta_i/2)];
R_i0 = quat2rotm(q_rot_i');
R_i0(1:3,3)=R_i0(1:3,3)*R_i0(2,2);
R_i0(1:3,2)=R_i0(1:3,2)*R_i0(2,2);
%восходящий узел
delta_O = O_coef-pi;
n_rot = [0; 0; 1];
q_rot = [cos(delta_O/2); n_rot*sin(delta_O/2)];
R_rot = quat2rotm(q_rot');
R_iO=R_rot*R_i0*R_rot';
t_end = get_t_end(k_coef);
%t_end = T_right(j1);
[pr_0,pv_0] = get_adjoint(a_input);
pr = R_iO^(-1)*pr_0;
pv = R_iO^(-1)*pv_0;
t0 = 0;
start_pos = [1 0 0];
start_vel = [0 1 0];
y0 = cat(2,start_pos,start_vel,pr',pv',0,0)';
tspan = linspace(t0,t0+t_end, round(t_end*20));
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
[t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);

J = y(end,13);       %функционал
AN_final = y(end,14)/(2*pi)%угловая дальность
%рисуем
t_orbit = linspace(t0,t0+T_earth/(24*3600), 1000);

st.t = t_orbit';
st.planet = 'Earth';
st.mode = orbits;
st.delta_omega = 0;

earth_traj = planetModel(st);
earth_traj=1e+3*earth_traj/ae;

t_orbit = linspace(t0,t0+T_mars/(24*3600), 1000);

st.t = t_orbit';
st.planet = planet_end;
st.mode = orbits;
st.delta_omega = 0;

mars_traj = planetModel(st);
mars_traj=1e+3*mars_traj/ae;
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
%plot3(1.52*[cos(delta_O+pi),0*cos(delta_O)],1.52*[sin(delta_O+pi),0*sin(delta_O)],[0,0])
%plot3(orbit_final(:, 1), orbit_final(:, 2), orbit_final(:, 3), 'r');
%plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k')
plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'k')
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(y(end, 1:3)',y(end, 4:6)',1);
[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
disp(['p=', num2str(p)])
disp(['a=', num2str(a)])
disp(['O=', num2str(O)])
disp(['i=', num2str(i)])
%disp(['o=', num2str(lonPer)])
%disp(['ex=', num2str(ex)])
%disp(['ey=', num2str(ey)])
disp('-------------')
%disp(['ix=', num2str(ix)])
%disp(['iy=', num2str(iy)])
plot3(y(end, 1), y(end, 2), y(end, 3),'bO')

plot3(mars_r_f(1), mars_r_f(2), mars_r_f(3),'rO')

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
function [c, ceq] = fmincon_find_a_i(x, i_target,a_target, O_target, start_pos,start_vel,t0,t_end)
c =[];
i_input = x(1);
a_input = x(2);
%наклонение
delta_i = atan(i_input);
n_rot_i = [1; 0; 0];
q_rot_i = [cos(delta_i/2); n_rot_i*sin(delta_i/2)];
R_i0 = quat2rotm(q_rot_i');
R_i0(1:3,3)=R_i0(1:3,3)*R_i0(2,2);
R_i0(1:3,2)=R_i0(1:3,2)*R_i0(2,2);
%восходящий узел
delta_O = O_target-pi;
n_rot = [0; 0; 1];
q_rot = [cos(delta_O/2); n_rot*sin(delta_O/2)];
R_rot = quat2rotm(q_rot');
R_iO=R_rot*R_i0*R_rot';

[pr_0,pv_0] = get_adjoint(a_input);
z = [pr_0;pv_0];

[i,O,p,a,J,AN] = integrate_iO(R_iO, start_pos,start_vel,z,t0,t_end);
ceq = [i-i_target,a-a_target];
end
function res = fsolve_find_t_end(AN_target, start_pos,start_vel,zf,t0,t_end)
R_iO = eye(3);
[i,O,p,a,J,AN] = integrate_iO(R_iO, start_pos,start_vel,zf,t0,t_end);
res = AN-AN_target;
end
function res = fsolve_find_a_fix(k, a_target,R_iO, start_pos,start_vel,t0,t_end)
[pr_0,pv_0] = get_adjoint(k);
z = [pr_0;pv_0];
[i,O,p,a,J,AN] = integrate_iO(R_iO, start_pos,start_vel,z,t0,t_end);
res = a-a_target;
end
function [i,O,p,a,J,AN] = integrate_iO(R_iO, start_pos,start_vel,zf,t0,t_end)
pr = R_iO^(-1)*zf(1:3);
pv = R_iO^(-1)*zf(4:6);
y0 = cat(2,start_pos,start_vel,pr',pv',0,0)';
tspan = linspace(t0,t0+t_end, 1000);
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
[t,y] = ode113(@(t,y) internalIntegration3D_withJ_TH(t,y), tspan,y0,options);

J = y(end,13);       %функционал
AN = y(end,14);%угловая дальность
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(y(end, 1:3)',y(end, 4:6)',1);
[p,ex,ey,ix,iy,L] = orb2equintoctial(p,eMag,i,O,lonPer,nu);
end
function [i,O,p,a] = integrate_iO_equinoctial(R_iO, start_pos,start_vel,zf,t0,t_end)
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
load('mat-files/Preprint2024_1_Mars.mat','t_end_Mars')
load('mat-files/Preprint2024_1_Venus.mat','t_end_Venus')
k1=0.72;
k2=1.52;
CT1 = ((t_end_Venus-40*2*pi)/log(k1)-(t_end_Mars-40*2*pi)/log(k2))/(k1-k2);
CT0 = (t_end_Venus-40*2*pi)/log(k1)-CT1*k1;
t_end = 40*2*pi+(CT1*k+CT0)*log(k);
end
function [pr_0,pv_0] = get_adjoint(k)
%Adjoint variables approximation
% a - angular distance, revolutions
% k - radii ratio
% d - multiplier of D coefficients, 0 or 1
%MARS
load('mat-files/Preprint2024_1_Mars.mat','Z0_Mars','dxdX')
load('mat-files/Preprint2024_1_Venus.mat','Z0_Venus')
z_Mars = inv(dxdX)'*Z0_Mars;
z_Venus = inv(dxdX)'*Z0_Venus;
pr_x_M = z_Mars(1);
pr_x_V = z_Venus(1);
pr_y_M = z_Mars(2);
pr_y_V = z_Venus(2);

pv_x_M = z_Mars(4);
pv_x_V = z_Venus(4);
pv_y_M = z_Mars(5);
pv_y_V = z_Venus(5);

k1=0.72;
k2=1.52;

%prx
p1 = pr_x_V;
p2 = pr_x_M;
g = (p2*power(k2,1/4)/log(k2) - p1*k1*power(k2,-3/4)/log(k1))/(1-power(k1/k2,3/4));
f = p1*k1/log(k1)-g*power(k1,3/4);
pr_x = @(k)(log(k)*(f/k+g*power(k,-1/4)));
%pry
p1 = pr_y_V;
p2 = pr_y_M;
g = (p2*k2^(3/2)/log(k2) - p1*k1*sqrt(k2)/log(k1))/(1-sqrt(k2/k1));
f = p1*k1/log(k1)-g/sqrt(k1);
pr_y = @(k)(log(k)*(f/k+g/(k^(3/2))));
%pvx
p1 = pv_x_V;
p2 = pv_x_M;
g = (p2*k2^(3/2)/log(k2) - p1*k1*sqrt(k2)/log(k1))/(1-sqrt(k2/k1));
f = p1*k1/log(k1)-g/sqrt(k1);
pv_x = @(k)(log(k)*(f/k+g/(k^(3/2))));
%pvy
p1 = pv_y_V;
p2 = pv_y_M;
g = (p2*power(k2,1/4)/log(k2) - p1*k1*power(k2,-3/4)/log(k1))/(1-power(k1/k2,3/4));
f = p1*k1/log(k1)-g*power(k1,3/4);
pv_y = @(k)(log(k)*(f/k+g*power(k,-1/4)));

pr_0 = [pr_x(k);pr_y(k);0.];
pv_0 = [pv_x(k);pv_y(k);0.]; 
end