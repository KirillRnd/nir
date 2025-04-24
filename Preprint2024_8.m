% в этом скрипте изучаем пространственные перелёты с помощью афинных
% преобразований
clear;
load('mat-files/Preprint2024_1_Mars.mat')
load('mat-files/Preprint2024_1_Venus.mat', 't_end_Venus')
%%
savefilename = 'mat-files/Preprint2024_8.mat';
k_coef_range=0.4:0.1:2.5;
i_coef_range=0:pi/120:pi/12;
O_coef_range=0:pi/15:2*pi;
N_count_1 = length(k_coef_range);
N_count_2 = length(i_coef_range);
N_count_3 = length(O_coef_range);

i_massive = zeros(N_count_1,N_count_2,N_count_3);
O_massive = zeros(N_count_1,N_count_2,N_count_3);
a_massive = zeros(N_count_1,N_count_2,N_count_3);
i_massive_error = zeros(N_count_1,N_count_2,N_count_3);
O_massive_error = zeros(N_count_1,N_count_2,N_count_3);
a_massive_error = zeros(N_count_1,N_count_2,N_count_3);
%%
load(savefilename)
for j1 = 1:N_count_1
    for j2 = 1:N_count_2
        for j3 = 1:N_count_3
        k_coef = k_coef_range(j1);
        i_coef = i_coef_range(j2);
        O_coef = O_coef_range(j3);
        
        if a_massive(j1,j2,j3) > 0
            continue
        end
        disp([num2str(j1),'   ',num2str(j2),'   ',num2str(j3)]);
        a_input = Ca0fromIK(i_coef, k_coef)+Ca1fromIK(i_coef, k_coef);
        delta_i_coef = Ci0fromK(a_input)+Ci1fromK(a_input)*cos(2*O_coef);
        delta_i = atan(i_coef/delta_i_coef);
        %изменение наклонения
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
        t_end = get_t_end(a_input);
        [pr_0,pv_0] = get_adjoint(a_input);
        z = [pr_0;pv_0];
        [i,O,p,a] = integrate_iO_equinoctial(R_iO, start_pos,start_vel,z,t0,t_end);
        i_massive(j1,j2,j3) = i;
        O_massive(j1,j2,j3) = O;
        a_massive(j1,j2,j3) = a;
        i_massive_error(j1,j2,j3) = i-i_coef;
        O_massive_error(j1,j2,j3) = O-O_coef;
        a_massive_error(j1,j2,j3) = a-k_coef;
        %save(savefilename)
        end
    end
end
%%
i_massive_error(isnan(i_massive_error))=0;
a_massive_error(isnan(a_massive_error))=0;
O_massive_error(isnan(O_massive_error))=0;

O_massive_error2=O_massive_error(:,2:end,:);
O_massive_error2(O_massive_error2<-pi)=2*pi+O_massive_error2(O_massive_error2<-pi);
O_massive_error2(O_massive_error2>pi)=-2*pi+O_massive_error2(O_massive_error2>pi);

i_std = sqrt(sum(i_massive_error.^2, 'all')/(numel(i_massive_error)-1));
a_std = sqrt(sum(a_massive_error.^2, 'all')/(numel(a_massive_error)-1));
O_std = sqrt(sum(O_massive_error2.^2, 'all')/(numel(O_massive_error2)-1));
disp(['i_std=', num2str(i_std)]);
disp(['a_std=', num2str(a_std)]);
disp(['O_std=', num2str(O_std)]);
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