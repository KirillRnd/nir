%plot checkMethod_ver2 result
figure(1);
plot3(0, 0, 0, 'k--o')
set(gca,'FontSize',14)
hold on;

t0 = t_start;
t_orbit = linspace(t0,t0+T_earth/(24*3600), 1000);

st.t = t_orbit';
st.planet = 'Earth';
st.mode = orbits;
st.delta_omega = omega;

earth_traj = planetModel(st);
earth_traj=earth_traj*1e+03/ae;
%Для KS
earth_traj_New = arrayfun(@(x,y,z)rotmZYX*[x, y, z]', earth_traj(:, 1),earth_traj(:, 2),earth_traj(:, 3),'UniformOutput',false);
earth_traj_New = cell2mat(earth_traj_New')';

t_orbit = linspace(t0,t0+T_mars/(24*3600), 1000);

st.t = t_orbit';
st.planet = planet_end;
st.mode = orbits;
st.delta_omega = omega;

mars_traj = planetModel(st);
mars_traj=mars_traj*1e+03/ae;

%Для KS
mars_traj_New = arrayfun(@(x,y,z)rotmZYX*[x, y, z]', mars_traj(:, 1),mars_traj(:, 2),mars_traj(:, 3),'UniformOutput',false);
mars_traj_New = cell2mat(mars_traj_New')';

plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k')
plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'k')

st.t = [t_start, T_end];
st.planet = 'Mars';
st.mode = orbits;
st.delta_omega = omega;

[mars_r_f, mars_v_f]=planetModel(st);
mars_r_f=mars_r_f'*1e+03;
mars_v_f=mars_v_f'*1e+03;
rr = checkMethod_results.rr;
rr_old = arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', rr(:, 1),rr(:, 2),rr(:, 3),'UniformOutput',false);
rr_old = cell2mat(rr_old')';

VV_old = arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', VV(:, 1),VV(:, 2),VV(:, 3),'UniformOutput',false);
VV_old = cell2mat(VV_old')';


plot3(rr_old(:, 1), rr_old(:, 2), rr_old(:, 3), 'cyan', 'LineWidth', 1);


plot3(rr_old(end, 1), rr_old(end, 2), rr_old(end, 3),'bO')
plot3(mars_r_f(1)/ae, mars_r_f(2)/ae,mars_r_f(3)/ae,'rO')

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