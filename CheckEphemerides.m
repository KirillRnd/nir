%отладочный скрипт
t0 = juliandate(2001,0,0);
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;
mug = 132712.43994*(10^6)*(10^(3*3));

plot3(0, 0, 0, 'y--o')
set(gca,'FontSize',14)
hold on;

a = 0.5;

th = linspace(0 , a*2*pi,10000);
t_orbit_E = linspace(t0,t0+T_earth/(24*3600), 10000);
[earth_traj, earth_traj_V] = planetEphemeris(t_orbit_E','SolarSystem','Earth','430');
earth_traj=earth_traj*1e+03;
earth_traj_V=earth_traj_V*1e+03;

eul = [0 pi/4 0];
rotmZYX = eul2rotm(eul);

earth_traj_New = arrayfun(@(x,y,z)rotmZYX*[x, y, z]', earth_traj(:, 1),earth_traj(:, 2),earth_traj(:, 3),'UniformOutput',false);
earth_traj_New = cell2mat(earth_traj_New')';

t_orbit_m = linspace(t0,t0+T_mars/(24*3600), 10000);
mars_traj = planetEphemeris(t_orbit_m','SolarSystem','Mars','430');
mars_traj=mars_traj*1e+03;

mars_traj_New = arrayfun(@(x,y,z)rotmZYX*[x, y, z]', mars_traj(:, 1),mars_traj(:, 2),mars_traj(:, 3),'UniformOutput',false);
mars_traj_New = cell2mat(mars_traj_New')';


mars_traj_ks = arrayfun(@(x, y, z) rToU([x,y,z]), mars_traj_New(:, 1),mars_traj_New(:, 2),mars_traj_New(:, 3),'UniformOutput',false);
mars_traj_ks = cell2mat(mars_traj_ks')';
% 

% 
earth_traj_ks = arrayfun(@(x, y, z) rToU([x,y,z]), earth_traj_New(:, 1),earth_traj_New(:, 2),earth_traj_New(:, 3),'UniformOutput',false);
earth_traj_ks = cell2mat(earth_traj_ks')';
plot3(earth_traj_ks(:, 1), earth_traj_ks(:, 2), earth_traj_ks(:, 3), 'k')
plot3(mars_traj_ks(:, 1), mars_traj_ks(:, 2), mars_traj_ks(:, 3), 'r')


axis equal

title('Проверка орбит в KS-переменных')
xlabel('u1')
ylabel('u2')
zlabel('u3')

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
hold off;