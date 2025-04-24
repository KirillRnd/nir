% Параметры GIF
%load('Corolev2.mat');
filename = 'NGPM2.gif'; % Имя файла GIF
delay = 0.05; % Задержка между кадрами (в секундах)

% Данные для графиков
x = linspace(0, 2*pi, 100);
%%
B=0.2721831;
A=0.5433279;
T_a  = @(a, a_rel)(a_rel-a_rel*B*log(a_rel))*a+(a_rel*B^2*log(a_rel));
a_OM = @(OM, a_rel,k)(OM+k+A*exp(-a_rel)*log(a_rel)*log(a_rel))/(A*(1+exp(-a_rel))*log(a_rel));
OM_a = @(a, a_rel)A*(a*(exp(a_rel)+1)-log(a_rel))*exp(-a_rel)*log(a_rel);

%%
i_max = 2;
J_array_1 = linspace(0,i_max,100);
J_array_2 = linspace(i_max,i_max,100);
J_array_3 = linspace(i_max,-i_max,100);
J_array_4 = linspace(-i_max,-i_max,100);
J_array_5 = linspace(-i_max,-0,100);
J_array = [J_array_1, J_array_2, J_array_3, J_array_4, J_array_5];
O_array_1 = linspace(0,0,100);
O_array_2 = linspace(0,2*pi,100);
O_array_3 = linspace(2*pi,2*pi,100);
O_array_4 = linspace(0,2*pi,100);
O_array_5 = linspace(2*pi,2*pi,100);
O_array = [O_array_1, O_array_2, O_array_3, O_array_4, O_array_5];
% Создание GIF
for k = 1:500
    AN_i = 10;
    %i_target = J_array(k);
    j = J_array(k);
    O_target = O_array(k);
    t0 = 0;
    d_coef = 1;
    a_rel = 1.52;
    
    a_rel_fix = Ca0fromIK(i_target, a_rel)+ Ca1fromIK(i_target, a_rel)*cos(2*O_target);
    % %a_rel_fix = a_rel;
    % a_input = a_rel_fix;
    Ci0 = (0.38159616342652647-0.0243395353632856*a_input-0.2767778010385592*exp(-a_input^2))*log(a_input);
    Ci1 = -a_input*exp(-7.8328443*sinh(a_input))+5.5964956*exp(-sinh((1.38935843796911*a_input^2+1.09321695795059)/(0.7103882*a_input^2+0.04272333)));
    delta_i_coef = Ci0+Ci1*cos(2*O_target);
    % if abs(delta_i_coef) > 0
    %     j = -i_target/delta_i_coef;
    % else
    %     j = 0;
    % end
    %T_i = 2*pi*(AN_i*(0.03772064567906565*a_rel_fix + 1.405055577366462)*log(a_rel_fix + 1)+(0.08297934358955021*a_rel_fix + 0.009838514316752341).*log(a_rel_fix).^2);
    T_i = 2*pi*T_a(AN_i, a_rel);
    OM_i = 2*pi*OM_a(AN_i, a_rel);
    
    [pr0,pv0] = get_initial_adjoint(AN_i,a_rel,d_coef);
    
    R_i = [
        1, 0, 0;
        0, 1/(1+j^2), -j/(1+j^2);
        0, j/(1+j^2), 1/(1+j^2);
    ];
    R_rot = [
        cos(O_target),-sin(O_target),0;
        sin(O_target),cos(O_target),0;
        0,0,1
        ];
    R_iO=R_rot*R_i*R_rot';
    pr_0 = (R_iO^(-1)*R^(-1)*pr0)';
    pv_0 = (R_iO^(-1)*R^(-1)*pv0)';
    y0 = cat(2,start_pos,start_vel,pr_0,pv_0)';
    tspan = linspace(t0,t0+T_i, AN_i*400);
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    st.t = [t_start, T_i];
    st.planet = planet_end;
    st.mode = orbits;
    st.delta_omega = OM_i;
    
    [mars_r_f, mars_v_f]=planetModel(st);
    mars_r_f=mars_r_f'*1e+03/r_unit;
    mars_v_f=mars_v_f'*1e+03/V_unit;
    res_r_f = norm(y(end, 1:3)-mars_r_f');
    res_v_f = norm(y(end, 4:6)-mars_v_f');
    
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
    
    st.t = t_orbit';
    st.planet = planet_end;
    st.mode = orbits;
    st.delta_omega = 0;
    
    mars_traj = planetModel(st);
    mars_traj=1e+3*mars_traj/ae;
    y0 = cat(2,y(1, 1:3),y(1, 4:6),[0,0,0],[0,0,0])';
    [t,orbit_initial] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    y0 = cat(2,y(end, 1:3),y(end, 4:6),[0,0,0],[0,0,0])';
    [t,orbit_final] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    THspan = linspace(0,2*pi,400);
    a_rel_traj = zeros(400,3);
    a_rel_traj(:,1) = a_rel*cos(THspan);
    a_rel_traj(:,2) = a_rel*sin(THspan);
    a_rel_traj(:,3) = 0*THspan;
    orbit_initial = zeros(400,3);
    orbit_initial(:,1) = 1*cos(THspan);
    orbit_initial(:,2) = 1*sin(THspan);
    orbit_initial(:,3) = 0*THspan;
    %выводим график
    figure(5);
    plot(0, 0, 'k--o');
    set(gca,'FontSize',14);
    hold on;
    set(gcf, 'color', 'white')
    
    %plot3(a_rel_traj(:, 1), a_rel_traj(:, 2), a_rel_traj(:, 3), 'k')
    plot3(y(1, 1), y(1, 2), y(1, 3), 'bO', 'LineWidth', 1);
    plot3(y(:, 1), y(:, 2),y(:, 3), 'red', 'LineWidth', 1.5);
    plot3(orbit_initial(:, 1), orbit_initial(:, 2), orbit_initial(:, 3), 'k--', 'LineWidth', 1.5);
    plot3(orbit_final(:, 1), orbit_final(:, 2),orbit_final(:, 3), 'k', 'LineWidth', 1.5);
    %plot(mars_traj(:, 1), mars_traj(:, 2), 'k')
    plot3(y(end, 1), y(end, 2), y(end, 3),'bO')
    % Настройка вида под углом
    azimuth = 45;      % Азимутальный угол (в градусах)
    elevation = 30;    % Угол возвышения (в градусах)
    view(azimuth, elevation);
    %plot3(mars_r_f(1), mars_r_f(2), mars_r_f(3),'rO')
    
    hold off;
    axis equal
    
    %title('Траектория КА')
    xlabel('x, AU')
    ylabel('y, AU')
    zlabel('z, AU')
    % if i_target <0
    %     i_target = -i_target;
    %     O_target = O_target-pi;
    % end
    title(['j=' sprintf('%.2f', j) ', \Omega=' sprintf('%.2f', O_target)]);
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    grid on;
    box off;
    xlim([-2,2])
    ylim([-2,2])
    zlim([-1,1])

    % Захват текущего кадра
    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);

    % Запись кадра в GIF
    if k == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', delay);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay);
    end
end

disp(['GIF сохранён как ' filename]);

function R = calculateRotMatrix(r0,v0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
e1 = r0;
e2 = v0;
e3 = cross(e1,e2);
R = [e1;e2;e3];
end
function out1 = Ca0fromIK(i, k)
  out1 = k + i.*(i.*k.*(i.*(0.75708324 - 0.22989681*k).*exp(k) + 1.637357) + 0.031784486);
end
function out1 = Ca1fromIK(i, k)
  out1 = i.^2.*(k + 2.518653).*(0.050697707*k.*(-6.1506586*i + k) - 0.8663172);
end
