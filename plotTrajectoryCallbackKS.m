function [] = plotTrajectoryCallbackKS(i_one, j_one, F_one, mega_st)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%Скрипт для вывода ровно одного графика для фикисрованной точки
%i_one = 184;
%j_one = 30;
%global i_one j_one;
% F_one = 2;
fig_number = 2;
if F_one>=0
    px=reshape(mega_st.PXevery(i_one,j_one,:,F_one),[1,8]);
    %неважен phi
    phi_new=mega_st.PHIevery(i_one,j_one, F_one);
else
    px=reshape(mega_st.PXevery(i_one,j_one,:),[1,8]);
    %неважен phi
    phi_new=mega_st.PHIevery(i_one,j_one);
end
%ds(j)
disp([i_one,j_one])
delta_s=mega_st.ds(j_one)*2*pi;
omega=mega_st.omega_space(i_one);
modifier_p=10^(-4-sqrt(delta_s));
x0_sec = [px delta_s/(2*pi) phi_new/(2*pi)];
integration_acc=1e-12;
calculate_condition=0;
t_start = 0;
[~, ~, ~, ~, ~, phi, t_end, ~, uu, rr, ~, ~, ~, ~, ~] = calculateTimeKS(t_start,delta_s, ...
    mega_st.rad,mega_st.UorR,mega_st.decreaseNonPsysical,modifier_p, ...
    mega_st.modifier_f,x0_sec,mega_st.eta, mega_st.case_traj, ...
    mega_st.planet_end, mega_st.display,mega_st.terminal_state,integration_acc,calculate_condition, mega_st.orbits, omega);
% 
%Выводим траекторию в параметрических переменных"
figure(fig_number);
ae = 149597870700;
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;

eul = [pi/12 pi/4 pi/12];
rotmZYX = eul2rotm(eul);
plot3(0, 0, 0, 'k--o')
set(gca,'FontSize',14)
hold on;
t0 = t_start;
t_orbit = linspace(t0,t0+T_earth/(24*3600), 1000);

st.t = t_orbit';
st.planet = 'Earth';
st.mode = mega_st.orbits;
st.delta_omega = omega;

earth_traj = planetModel(st);
earth_traj=earth_traj*1e+03/ae;
%Для KS
earth_traj_New = arrayfun(@(x,y,z)rotmZYX*[x, y, z]', earth_traj(:, 1),earth_traj(:, 2),earth_traj(:, 3),'UniformOutput',false);
earth_traj_New = cell2mat(earth_traj_New')';

t_orbit = linspace(t0,t0+T_mars/(24*3600), 1000);

st.t = t_orbit';
st.planet = mega_st.planet_end;
st.mode = mega_st.orbits;
st.delta_omega = omega;

mars_traj = planetModel(st);
mars_traj=mars_traj*1e+03/ae;

%Для KS
mars_traj_New = arrayfun(@(x,y,z)rotmZYX*[x, y, z]', mars_traj(:, 1),mars_traj(:, 2),mars_traj(:, 3),'UniformOutput',false);
mars_traj_New = cell2mat(mars_traj_New')';
phi_MARS=atan2(uu(end,4),uu(end,1));
mars_traj_ks = arrayfun(@(r1, r2, r3) rToU([r1,r2,r3], phi_MARS), mars_traj_New(:, 1),mars_traj_New(:, 2),mars_traj_New(:, 3),'UniformOutput',false);
mars_traj_ks = cell2mat(mars_traj_ks')';

earth_traj_ks = arrayfun(@(r1, r2, r3) rToU([r1,r2,r3], phi_new), earth_traj_New(:, 1),earth_traj_New(:, 2),earth_traj_New(:, 3),'UniformOutput',false);
earth_traj_ks = cell2mat(earth_traj_ks')';
plot3(earth_traj_ks(:, 1), earth_traj_ks(:, 2), earth_traj_ks(:, 3), 'k')
plot3(mars_traj_ks(:, 1), mars_traj_ks(:, 2), mars_traj_ks(:, 3), 'r')
plot3(uu(:, 1), uu(:, 2), uu(:, 3), 'b', 'LineWidth', 2.5);
axis equal

title('Траектория КА KS')
xlabel('u1')
ylabel('u2')
zlabel('u3')
grid on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
hold off;
end