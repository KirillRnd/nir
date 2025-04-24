% в этом скрипте изучаем пространственные перелёты с помощью афинных
% преобразований
clear;
load('mat-files/Preprint2024_1_Mars.mat')
load('mat-files/Preprint2024_1_Venus.mat', 't_end_Venus')
%%
savefilename = 'mat-files/Preprint2024_ecc_3.mat';
k_coef_range=0.4:0.1:2.5;
dm_coef_range=-0.10:0.01:0.11;
o_coef_range=0:pi/15:2*pi;
N_count_1 = length(k_coef_range);
N_count_2 = length(dm_coef_range);
N_count_3 = length(o_coef_range);

ex_massive = zeros(N_count_1,N_count_2,N_count_3);
ey_massive = zeros(N_count_1,N_count_2,N_count_3);
o_massive = zeros(N_count_1,N_count_2,N_count_3);
p_massive = zeros(N_count_1,N_count_2,N_count_3);
a_massive = zeros(N_count_1,N_count_2,N_count_3);
%%
load(savefilename)
for j1 = 1:N_count_1
    for j2 = 1:N_count_2
        for j3 = 1:N_count_3
        k_coef = k_coef_range(j1);
        dm_coef = dm_coef_range(j2);
        o_target = o_coef_range(j3);
        
        if a_massive(j1,j2,j3) > 0
            continue
        end
        disp([num2str(j1),'   ',num2str(j2),'   ',num2str(j3)]);
        %a_input = Ca0fromIK(i_coef, k_coef)+Ca1fromIK(i_coef, k_coef);
        a_input=k_coef;

        theta = 0.642;
        m = dm_coef-sin(theta);%X_rotated
        n = cos(theta); %Y_rotated 
        

        ix = 0.0;
        iy = 0;
        L  = 0;
        ey = 0;
        ex = m*cos(theta)+n*sin(theta);
        a_param = -m*sin(theta)+n*cos(theta);
        p = a_param*(1-ex^2-ey^2);
        X_2 = equitoctial2decart([p;ex;ey;ix;iy;L], mug);
        start_pos_3 = X_2(1:3)';
        start_vel_3 = X_2(4:6)';
        
        R_e = calculateRotMatrix(start_pos_3,start_vel_3);


        R_rot = [cos(o_target),-sin(o_target),0;
                sin(o_target),cos(o_target),0;
                0,0,1];
        R_eo=R_rot*R_e*R_rot';

        t_end = get_t_end(a_input);
        [pr_0,pv_0] = get_adjoint(a_input);
        z = [pr_0;pv_0];
        [ex,ey,lonPer,p,a] = integrate_iO_equinoctial(R_eo, start_pos,start_vel,z,t0,t_end);
        ex_massive(j1,j2,j3) = ex;
        ey_massive(j1,j2,j3) = ey;
        o_massive(j1,j2,j3) = lonPer;
        p_massive(j1,j2,j3) = p;
        a_massive(j1,j2,j3) = a;
        end
    save(savefilename)
    end
    
end

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