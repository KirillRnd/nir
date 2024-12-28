%Скрипт для вывода ровно одного графика для фикисрованной точки
%i_one = 184;
%j_one = 30;
%global i_one j_one;
i_one = 181;
j_one = 73;
%j_one = 72;
F_one = 1;
fig_number = 2;
px_new=reshape(PXevery(i_one,j_one,:,F_one),[1,8]);
%неважен phi
phi_new=PHIevery(i_one,j_one, F_one);
%ds(j)
disp([i_one,j_one, F_one])
delta_s=ds(j_one)*2*pi;
omega=omega_space(i_one);
modifier_p=10^(-4-sqrt(delta_s));
x0_sec = [px_new delta_s/(2*pi) phi_new/(2*pi)];
integration_acc=1e-12;
calculate_condition=0;
t_start = 0;
[~, ~, ~, ~, ~, ~, t_end, ~, ~, rr, ~, ~, ~, ~, ~] = calculateTimeKS(t_start,delta_s, ...
    rad,UorR,decreaseNonPsysical,modifier_p, ...
    modifier_f,x0_sec,eta, case_traj, ...
    planet_end, display,terminal_state,integration_acc,calculate_condition, orbits, omega);
% 
% 
ae = 149597870700;
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;

eul = [pi/12 pi/4 pi/12];
rotmZYX = eul2rotm(eul);

%Проверка "на глаз"
figure(fig_number);
plot3(0, 0, 0, 'k--o', 'HandleVisibility','off')
set(gca,'FontSize',14)
hold on;

%вывод траектории Земли
t0 = t_start;
t_orbit = linspace(t0,t0+T_earth/(24*3600), 1000);

st.t = t_orbit';
st.planet = 'Earth';
st.mode = orbits;
st.delta_omega = omega;

earth_traj = planetModel(st);
earth_traj=earth_traj*1e+03/ae;
%Для KS
% earth_traj_New = arrayfun(@(x,y,z)rotmZYX*[x, y, z]', earth_traj(:, 1),earth_traj(:, 2),earth_traj(:, 3),'UniformOutput',false);
% earth_traj_New = cell2mat(earth_traj_New')';

plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k', 'HandleVisibility','off')

%вывод траектории Марса
t_orbit = linspace(t0,t0+T_mars/(24*3600), 1000);

st.t = t_orbit';
st.planet = planet_end;
st.mode = orbits;
st.delta_omega = omega;

mars_traj = planetModel(st);
mars_traj=mars_traj*1e+03/ae;

% mars_traj_New = arrayfun(@(x,y,z)rotmZYX*[x, y, z]', mars_traj(:, 1),mars_traj(:, 2),mars_traj(:, 3),'UniformOutput',false);
% mars_traj_New = cell2mat(mars_traj_New')';

plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'k', 'HandleVisibility','off')

%положение Марса в конечной точке
st.t = [t_start, t_end];
st.planet = planet_end;
st.mode = orbits;
st.delta_omega = omega;

[mars_r_f, ~]=planetModel(st);
mars_r_f=mars_r_f'*1e+03;
%mars_v_f=mars_v_f'*1e+03;

rr_old = arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', rr(:, 1),rr(:, 2),rr(:, 3),'UniformOutput',false);
rr_old = cell2mat(rr_old')';

% VV_old = arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', VV(:, 1),VV(:, 2),VV(:, 3),'UniformOutput',false);
% VV_old = cell2mat(VV_old')';


FaceAlpha=0.6;
red_color = [1 0 0 FaceAlpha];
cyan_color = [0 1 1 FaceAlpha];
plot3(rr_old(:, 1), rr_old(:, 2), rr_old(:, 3), 'Color',  cyan_color, 'LineWidth', 2, 'DisplayName','Семейство 2');
plot3(rr_old_for_draw(:, 1), rr_old_for_draw(:, 2), rr_old_for_draw(:, 3), 'Color',  red_color, 'LineWidth', 2, 'DisplayName','Семейство 3');


% a_ks_old= arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', a_ks(:, 1),a_ks(:, 2),a_ks(:, 3),'UniformOutput',false);
% a_ks_old = cell2mat(a_ks_old')';
% 
% a_scale=3e-01/mean(vecnorm(a_ks, 2, 2));
% %a_scale=0;
% d = 24*3600;
% idxes=1;
% for i=1:ceil(t(end)/d)
%     ix = find(t>d*i*10, 1);
%     idxes=[idxes, ix];
% end   
%тяга
% for i = idxes
%     plot3([rr_old(i, 1), rr_old(i, 1)+a_scale*a_ks_old(i, 1)], [rr_old(i, 2), rr_old(i, 2)+...
%         a_scale*a_ks_old(i, 2)],[rr_old(i, 3), rr_old(i, 3)+a_scale*a_ks_old(i, 3)],'k')
% end

%plot3(rr_old(end, 1), rr_old(end, 2), rr_old(end, 3),'bO')
%plot3(mars_r_f(1)/ae, mars_r_f(2)/ae,mars_r_f(3)/ae,'rO')

%plot3(rr_cont(:, 1)/ae, rr_cont(:, 2)/ae, rr_cont(:, 3)/ae, 'g', 'LineWidth', 2.5);
%эти две точки должны находиться рядом
%plot3(rr_cont(500, 1)/ae, rr_cont(500, 2)/ae, rr_cont(500, 3)/ae, 'gO', 'LineWidth', 2.5);
%plot3(rr_old(500, 1), rr_old(500, 2), rr_old(500, 3), 'bO', 'LineWidth', 2.5);
axis equal

title('Траектория КА')
xlabel('x, a.e.')
ylabel('y, a.e.')
zlabel('z, a.e.')
view(0,90)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
grid on;
box off;
hold off;
legend;