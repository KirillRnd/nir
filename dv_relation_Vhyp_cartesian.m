%Ищем оптимум по гипизбытку
load('D:\MATLAB\ipm\mat-files\16-Jun-2024.mat')
%%
phi0 = 0;
mug = 1;

st.t = t_start;
st.planet = 'Earth';
st.mode = 'Flat';
st.delta_omega = omega_space(1,1);
st.a_rel = 1;

ae = 149597870700; %км
%mug = 132712.43994*(10^6)*(10^9); %км3с−2

[r0, v0]=planetModel(st);
eul = [pi/12 pi/4 pi/12];
rotmZYX = eul2rotm(eul);

R0 = [rotmZYX*r0'/ae; 0]*1e+03;

PVevery_Vhyp=zeros(L2,L4, 3);
PRevery_Vhyp=zeros(L2,L4, 3);
for i = 1:L2
    for j = 1:L4
        hi = HIevery_Vhyp(i,j);
        dV_value=dV_range(j);
        dV_add = dV_value*rotmZYX*[cos(hi), sin(hi), 0]'/V_unit;
        V0 = [rotmZYX*v0'/V_unit; 0]*1e+03+[dV_add; 0];
        phi0 = PHIevery_Vhyp(i,j);
        u0 = rToU(R0, phi0);
        w0 = vFromV(V0,R0,mug,phi0);
        ortdgduv = get_ortdgduv(u0,w0);
        dFdX = get_dgduv(u0,w0);
        vecPX1 = reshape(PXevery_Vhyp(i,j,:), [1,8]);
        vecPV_cartesian = a_reactive(u0,w0,vecPX1(1:4)',vecPX1(5:8)');
        b = vecPX1'-dFdX' * [0; 0; 0; vecPV_cartesian(1:3)]; %TODO переписать в символьных
        vecPR_cartesian = dFdX(1:3,:)'\b;

        vecPV_cartesian = rotmZYX^(-1)*vecPV_cartesian(1:3);
        vecPR_cartesian = rotmZYX^(-1)*vecPR_cartesian;
        PVevery_Vhyp(i,j,:) = vecPV_cartesian';
        PRevery_Vhyp(i,j,:) = vecPR_cartesian';
    end
end
%% визуализируем pr и pv
figure(1)
[X, Y] = meshgrid(ds(1:105), dV_range/V_unit);
s = surf(X,Y,PVevery_Vhyp(1:105,:,1)');
s.EdgeColor = 'none';
xlabel("Угловая дальность, витки")
ylabel("Добавка скорости, безразм.")
zlabel("p_v_x")
figure(2)
s = surf(X,Y,PVevery_Vhyp(1:105,:,2)');
s.EdgeColor = 'none';
xlabel("Угловая дальность, витки")
ylabel("Добавка скорости, безразм.")
zlabel("p_v_y")
figure(3)
s = surf(X,Y,PRevery_Vhyp(1:105,:,1)');
s.EdgeColor = 'none';
xlabel("Угловая дальность, витки")
ylabel("Добавка скорости, безразм.")
zlabel("p_r_x")
figure(4)
s = surf(X,Y,PRevery_Vhyp(1:105,:,2)');
s.EdgeColor = 'none';
xlabel("Угловая дальность, витки")
ylabel("Добавка скорости, безразм.")
zlabel("p_r_y")

%%
figure(26);
dV_mesh = repmat(dV_range/V_unit, 105, 1);%.*sin(HIevery_Vhyp(1:105,:));
s = surf(PRevery_Vhyp(1:105,:,1), PVevery_Vhyp(1:105,:,2), dV_mesh, 'FaceColor','cyan', 'FaceAlpha',FaceAlpha, 'HandleVisibility','off');
s.EdgeColor = 'none';
hold on;
for i = [1,41,105]
     pos_j = 20;
     plot3(PRevery_Vhyp(i,1:pos_j,1), PVevery_Vhyp(i,1:pos_j,2), dV_mesh(i,1:pos_j), 'black', 'HandleVisibility','off')
     plot3(PRevery_Vhyp(i,pos_j+3:end,1), PVevery_Vhyp(i,pos_j+3:end,2), dV_mesh(i,pos_j+3:end), 'black', 'HandleVisibility','off')
     text_label = num2str(ds(i));
     h = text((PRevery_Vhyp(i,pos_j,1)+PRevery_Vhyp(i,pos_j+3,1))/2, ...
         (PVevery_Vhyp(i,pos_j,2)+PVevery_Vhyp(i,pos_j+3,2))/2, ...
         (dV_mesh(i,pos_j)+dV_mesh(i,pos_j+3))/2, text_label,'HorizontalAlignment','center',...
    'VerticalAlignment','Bottom','FontName','consolas','FontSize',11); 

end
plot3(PR_best_minimum_line(:,1),PV_best_minimum_line(:,2),vhyp_best_minimum_line(6,1:105)/V_unit,'r','LineWidth', 1.5)
hold off;
xlabel('Pr x')
ylabel('Pv y')
zlabel('Отлётная скорость')
%%
figure(27);
dV_mesh = repmat(dV_range/V_unit, 105, 1);%.*sin(HIevery_Vhyp(1:105,:));
s = surf(PRevery_Vhyp(1:105,:,2), PVevery_Vhyp(1:105,:,1), dV_mesh, 'FaceColor','cyan', 'FaceAlpha',FaceAlpha, 'HandleVisibility','off');
s.EdgeColor = 'none';
hold on;
for i = [1,41,105]
     pos_j = 20;
     plot3(PRevery_Vhyp(i,1:pos_j,2), PVevery_Vhyp(i,1:pos_j,1), dV_mesh(i,1:pos_j), 'black', 'HandleVisibility','off')
     plot3(PRevery_Vhyp(i,pos_j+3:end,2), PVevery_Vhyp(i,pos_j+3:end,1), dV_mesh(i,pos_j+3:end), 'black', 'HandleVisibility','off')
     text_label = num2str(ds(i));
     h = text((PRevery_Vhyp(i,pos_j,2)+PRevery_Vhyp(i,pos_j+3,2))/2, ...
         (PVevery_Vhyp(i,pos_j,1)+PVevery_Vhyp(i,pos_j+3,1))/2, ...
         (dV_mesh(i,pos_j)+dV_mesh(i,pos_j+3))/2, text_label,'HorizontalAlignment','center',...
    'VerticalAlignment','Bottom','FontName','consolas','FontSize',11); 

end
plot3(PR_best_minimum_line(:,2),PV_best_minimum_line(:,1),vhyp_best_minimum_line(6,1:105)/V_unit,'r','LineWidth', 1.5)
hold off;
xlabel('Pr y')
ylabel('Pv x')
zlabel('Отлётная скорость')
%%
figure(5)
[X, Y] = meshgrid(ds(1:105), dV_range/V_unit);
s = surf(X,Y,PRevery_Vhyp(1:105,:,2)');
s.EdgeColor = 'none';
xlabel("Угловая дальность, витки")
ylabel("Добавка скорости, безразм.")
zlabel("p_r_y")
%% ищем гладкое оптимальное решение
T_earth_days = 365.256363004;
N_count = 105;
PRevery_Vhyp_best=zeros(N_count, 3); %базис-вектор
Jevery_Vhyp_best=zeros(N_count, 1);     %функционал
Tevery_Vhyp_best = zeros(N_count, 1);   %время перелёта
ANevery_Vhyp_best = zeros(N_count, 1);  %угловая дальность
OMevery_Vhyp_best = zeros(N_count, 1);  %разность фаз
Vadd_Vhyp_best = zeros(N_count, 3);  %добавка скорости
%% ищем ещё более гладкое оптимальное решение
PRevery_Vhyp_best_2=PRevery_Vhyp_best; %базис-вектор
Jevery_Vhyp_best_2=Jevery_Vhyp_best;     %функционал
Tevery_Vhyp_best_2 = Tevery_Vhyp_best;   %время перелёта
ANevery_Vhyp_best_2 = ANevery_Vhyp_best;  %угловая дальность
OMevery_Vhyp_best_2 = OMevery_Vhyp_best;  %разность фаз
Vadd_Vhyp_best_2 = Vadd_Vhyp_best;  %добавка скорости
%%
for i=1:105
    disp(i)
    % j_local = Jevery_Vhyp(i,2:end);
    % j_min = min(j_local,[],'all');
    % j_min = find(j_local==j_min);
    % j_min = j_min+1;
    % 
    % j=j_min;
    % hi = HIevery_Vhyp(i,j);
    dV_value=dV_range(j-1);
    %dV_add = dV_value*[cos(hi), sin(hi), 0]'/V_unit;
    dV_add = Vadd_Vhyp_best(i,:)';
    R0 = [1;0;0];
    V0 = [0;1;0]+dV_add;
    %PR0 = PRevery_Vhyp(i,j,:);
    PR0 = PRevery_Vhyp_best(i,:)';
    %PR0 = reshape(PR0,1,[])';
    PR0(3)=0;
    PV0 = [0;0;0]; %из понтрягина для задачи пролёта
    %omega0 = OMevery_Vhyp(i,j);

    omega0 = OMevery_Vhyp_best(i);
    t0=0;
    %t_end = 2*pi*Tevery_Vhyp(i,j)/T_earth_days;
    t_end = Tevery_Vhyp_best(i);
    y0 = [R0;V0;PR0;PV0];
    z0 = [V0;PR0];
    
    minimize_delta_omega = @(delta_omega) find_close_solution(delta_omega, omega0, y0, z0, t0,t_end);
    angle_range = pi/8;
    [delta_omega, Jt_end_min] = fminbnd(minimize_delta_omega, -angle_range, angle_range);
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
    
    v0=zf(1:3);
    pr0=zf(4:6);
    y0(4:6)=v0;
    y0(7:9)=pr0;
    
    tspan = linspace(t0,t0+t_end, round(t_end*100));
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    % plot(y(:,1),y(:,2))
    % hold on;
    % plot(0,0, 'kO')
    % plot(y(end,1),y(end,2), 'bO')
    % plot(yf(end,1),yf(end,2), 'rO')
    % hold off;
    % axis equal
    % grid;
    
    a_vec=vecnorm(y(:,10:12), 2, 2).^2;
    Jt = cumtrapz(t, a_vec)/(2);
    J = Jt(end);
    AN = calculate_angular_distance(y(:,1:3));
    dV_add = v0-[0;1;0];

    J_old = Jevery_Vhyp_best(i);
    %if J<J_old
    %записываем в файл
    PRevery_Vhyp_best(i,:)=pr0; %базис-вектор
    Jevery_Vhyp_best(i)=J;     %функционал
    Tevery_Vhyp_best(i) = t_end;   %время перелёта
    ANevery_Vhyp_best(i) = AN;  %угловая дальность
    OMevery_Vhyp_best(i) = omega_min;  %разность фаз
    Vadd_Vhyp_best(i,:) = dV_add;  %добавка скорости 
    %end
    savefilename = join(['mat-files/special',date]);
    %save(savefilename);
end
%%
load('Corolev.mat')
% Порог для определения точки разрыва
threshold = 1;

% Вычисление разностей между соседними элементами
diffs = abs(diff(ANevery_Vhyp_best));

% Индексы, где разность превышает порог
break_points = find(diffs > threshold);

figure(1);

% Построение функции с учётом разрывов
for i = 1:length(break_points) + 1
    if i == 1
        % Первая часть до первого разрыва
        plot(Tevery_Vhyp_best(1:break_points(i))/(2*pi), ANevery_Vhyp_best(1:break_points(i))/(2*pi), 'b', 'LineWidth', 1.5);
    elseif i == length(break_points) + 1
        % Последняя часть после последнего разрыва
        plot(Tevery_Vhyp_best(break_points(i-1)+1:N_count)/(2*pi), ANevery_Vhyp_best(break_points(i-1)+1:end)/(2*pi), 'b', 'LineWidth', 1.5);
    else
        hold on;
        % Промежуточные части
        plot(Tevery_Vhyp_best(break_points(i-1)+1:break_points(i))/(2*pi), ANevery_Vhyp_best(break_points(i-1)+1:break_points(i))/(2*pi), 'b', 'LineWidth', 1.5);
    end
end
hold off;
grid;
xlabel("Время перелёта, годы")
ylabel("Угловая дальность, витки")
%%
figure(2)

% Построение функции с учётом разрывов
for i = 1:length(break_points) + 1
    if i == 1
        % Первая часть до первого разрыва
        plot(Tevery_Vhyp_best(1:break_points(i))/(2*pi), Vadd_Vhyp_best(1:break_points(i),1), 'b', 'LineWidth', 1.5, 'DisplayName', '\deltav_x');
    elseif i == length(break_points) + 1
        % Последняя часть после последнего разрыва
        plot(Tevery_Vhyp_best(break_points(i-1)+1:N_count)/(2*pi), Vadd_Vhyp_best(break_points(i-1)+1:end,1), 'b', 'LineWidth', 1.5, 'HandleVisibility','off');
    else
        hold on;
        % Промежуточные части
        plot(Tevery_Vhyp_best(break_points(i-1)+1:break_points(i))/(2*pi), Vadd_Vhyp_best(break_points(i-1)+1:break_points(i),1), 'b', 'LineWidth', 1.5, 'HandleVisibility','off');
    end
end

%plot(Tevery_Vhyp_best/(2*pi),Vadd_Vhyp_best(:,1), 'DisplayName', 'v_x')
hold on;
Vadd_norm = sqrt(Vadd_Vhyp_best(:,1).^2+Vadd_Vhyp_best(:,2).^2);
% Построение функции с учётом разрывов
for i = 1:length(break_points) + 1
    if i == 1
        % Первая часть до первого разрыва
        plot(Tevery_Vhyp_best(1:break_points(i))/(2*pi), Vadd_Vhyp_best(1:break_points(i),2), 'green', 'LineWidth', 1.5, 'DisplayName', '\deltav_y');
    elseif i == length(break_points) + 1
        % Последняя часть после последнего разрыва
        plot(Tevery_Vhyp_best(break_points(i-1)+1:N_count)/(2*pi), Vadd_Vhyp_best(break_points(i-1)+1:end,2), 'green', 'LineWidth', 1.5, 'HandleVisibility','off');
    else
        % Промежуточные части
        plot(Tevery_Vhyp_best(break_points(i-1)+1:break_points(i))/(2*pi), Vadd_Vhyp_best(break_points(i-1)+1:break_points(i),2), 'green', 'LineWidth', 1.5, 'HandleVisibility','off');
    end
end
%plot(Tevery_Vhyp_best/(2*pi),Vadd_Vhyp_best(:,2), 'DisplayName', 'v_y')
for i = 1:length(break_points) + 1
    if i == 1
        % Первая часть до первого разрыва
        plot(Tevery_Vhyp_best(1:break_points(i))/(2*pi), Vadd_norm(1:break_points(i)), 'k--', 'LineWidth', 1.5, 'DisplayName', '||\deltav|| пролёт');
    elseif i == length(break_points) + 1
        % Последняя часть после последнего разрыва
        plot(Tevery_Vhyp_best(break_points(i-1)+1:N_count)/(2*pi), Vadd_norm(break_points(i-1)+1:end), 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');
    else
        hold on;
        % Промежуточные части
        plot(Tevery_Vhyp_best(break_points(i-1)+1:break_points(i))/(2*pi), Vadd_norm(break_points(i-1)+1:break_points(i)), 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');
    end
end
%plot(Tevery_Vhyp_best/(2*pi),Vadd_norm, 'DisplayName', '|v|')
plot(vhyp_best_minimum_line(1,1:105)/T_earth_days, vhyp_best_minimum_line(6,1:105)/V_unit, 'r', 'LineWidth', 1.5, 'DisplayName', '||\deltav|| опт. гипызбыток')
hold off;
grid;
legend('Location','best');
xlabel("Время перелёта, годы")
ylabel("Отлётная скорость, безразм.")
%%
figure(4)
Vadd_angle = atan2(Vadd_Vhyp_best(:,2),Vadd_Vhyp_best(:,1));
plot(Tevery_Vhyp_best/(2*pi),Vadd_angle)
grid;
xlabel("Время перелёта, годы")
ylabel("Угол скорости, рад.")
%%
J_unit = 176.62545106129272;
figure(5)
plot(Tevery_a(:,113)/T_earth_days,Jevery_a(:,113)/J_unit,'b--', 'LineWidth', 1.5, 'DisplayName', 'J задача 1*')
hold on;
plot(vhyp_best_minimum_line(1,1:105)/T_earth_days, vhyp_best_minimum_line(5,1:105)/J_unit,'cyan', 'LineWidth', 1.5, 'DisplayName', 'J задача 2')
plot(Tevery_Vhyp_best/(2*pi),Jevery_Vhyp_best,'k--', 'LineWidth', 1.5, 'DisplayName', 'J задача 3')
hold off;
xlim([1,8])
legend('Location','best');
grid;
xlabel("Время перелёта, годы")
ylabel("Целевой функционал, безразм")
%%
figure(6)
plot(vhyp_best_minimum_line(1,1:105)/T_earth_days,vhyp_best_minimum_line(7,1:105), 'DisplayName', 'OM old')
hold on;
plot(Tevery_Vhyp_best/(2*pi),OMevery_Vhyp_best, 'DisplayName', 'OM new')
hold off;
legend('Location','best');
grid;
xlabel("Время перелёта, годы")
%%
figure(7)
%plot(Tevery_Vhyp_best/(2*pi),Vadd_Vhyp_best(:,2), 'DisplayName', 'v_y')
for i = 1:length(break_points) + 1
    if i == 1
        % Первая часть до первого разрыва
        plot(Tevery_Vhyp_best(1:break_points(i))/(2*pi), PRevery_Vhyp_best(1:break_points(i),1), 'k', 'LineWidth', 1.5, 'DisplayName', 'p_r_x');
    elseif i == length(break_points) + 1
        % Последняя часть после последнего разрыва
        plot(Tevery_Vhyp_best(break_points(i-1)+1:N_count)/(2*pi), PRevery_Vhyp_best(break_points(i-1)+1:end,1), 'k', 'LineWidth', 1.5, 'HandleVisibility','off');
    else
        hold on;
        % Промежуточные части
        plot(Tevery_Vhyp_best(break_points(i-1)+1:break_points(i))/(2*pi), PRevery_Vhyp_best(break_points(i-1)+1:break_points(i),1), 'k', 'LineWidth', 1.5, 'HandleVisibility','off');
    end
end
%plot(Tevery_Vhyp_best/(2*pi),PRevery_Vhyp_best(:,1), 'DisplayName', 'p_r_x')
hold on;
plot(Tevery_Vhyp_best/(2*pi),PRevery_Vhyp_best(:,2),'cyan', 'LineWidth', 1.5, 'DisplayName', 'p_r_y')
hold off;
legend('Location','best');
grid;
ylabel("Сопряжённые переменные в t_0")
xlabel("Время перелёта, годы")
%%

figure(8)
plot(ANevery_a(:,113)/(2*pi), PRevery_a(:,113,2))
hold on;
for i = 1:length(break_points) + 1
    if i == 1
        % Первая часть до первого разрыва
        plot(ANevery_Vhyp_best(1:break_points(i))/(2*pi), PRevery_Vhyp_best(1:break_points(i),1), 'k', 'LineWidth', 1.5, 'DisplayName', 'p_r_x');
    elseif i == length(break_points) + 1
        % Последняя часть после последнего разрыва
        plot(ANevery_Vhyp_best(break_points(i-1)+1:N_count)/(2*pi), PRevery_Vhyp_best(break_points(i-1)+1:end,1), 'k', 'LineWidth', 1.5, 'HandleVisibility','off');
    else
        % Промежуточные части
        plot(ANevery_Vhyp_best(break_points(i-1)+1:break_points(i))/(2*pi), PRevery_Vhyp_best(break_points(i-1)+1:break_points(i),1), 'k', 'LineWidth', 1.5, 'HandleVisibility','off');
    end
end


for i = 1:length(break_points)
    plot([ANevery_Vhyp_best(break_points(i))/(2*pi), ANevery_Vhyp_best(break_points(i))/(2*pi)], [-1e-3,4e-3], 'r--')
end
hold off;
grid;
%%
figure(8)
plot(vhyp_best_minimum_line(1,1:105)/T_earth_days, PR_best_minimum_line(:,1), 'LineWidth', 1.5, 'DisplayName', 'p_r_x')
hold on;
plot(vhyp_best_minimum_line(1,1:105)/T_earth_days, PR_best_minimum_line(:,2), 'LineWidth', 1.5, 'DisplayName', 'p_r_y')
plot(vhyp_best_minimum_line(1,1:105)/T_earth_days, PV_best_minimum_line(:,1), 'LineWidth', 1.5, 'DisplayName', 'p_v_x')
plot(vhyp_best_minimum_line(1,1:105)/T_earth_days, PV_best_minimum_line(:,2), 'LineWidth', 1.5, 'DisplayName', 'p_v_y')
hold off;
grid;
legend('Location','best');
ylabel("Сопряжённые переменные в t_0")
xlabel("Время перелёта, годы")

%%
i = 10;
j = 113;
prx = PRevery_a(i,j,1);
pry = PRevery_a(i,j,2);
pvx = PVevery_a(i,j,1);
pvy = PVevery_a(i,j,2);

vx0 = 0;
vy0 = 1;
vx10 = 0;
vy10 = 1;
rx10 = 1;
ry10 = 0;

prx_new = PR_best_minimum_line(i,1);
p = [prx, pry, pvx, pvy];
v_param_0 = [vx0;vy0;vx10;vy10;rx10;ry10];
fsolve_transform_fun=@(v_param)fsolve_transform(v_param, p);
[v_param_f, res]= fsolve(fsolve_transform_fun, v_param_0)
%%
function res = fsolve_transform(v_param, p)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
vx = v_param(1);
vy = v_param(2);
vx1 = v_param(3);
vy1 = v_param(4);
rx1 = v_param(5);
ry1 = v_param(6);
pr = p(1:2)';
pv = p(3:4)';
R1 = [1, vx; 0, vy];
R2 = [1, vx1; 0, vy1];

pr_new = inv(R2)'*inv(R1)'*pr;
pv_new = inv(R2)'*inv(R1)'*pv;
res = [pr_new(1); pv_new(1); pv_new(2);];
end
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
    
    v0=zf(1:3);
    pr0=zf(4:6);
    y0(4:6)=v0;
    y0(7:9)=pr0;
    
    tspan = linspace(t0,t0+t_end, round(t_end*100));
    options = odeset('AbsTol',1e-12);
    options = odeset(options,'RelTol',1e-12);   
    [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    
    
    a_vec=vecnorm(y(:,10:12), 2, 2).^2;
    Jt = cumtrapz(t, a_vec)/(2);
    Jt_end = Jt(end);
end
function res = fsolve_traj(z,y0,yf,t0,t_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tspan = linspace(t0,t0+t_end, 10);
v0=z(1:3);
pr0=z(4:6);
y0(4:6)=v0;
y0(7:9)=pr0;
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
[~,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
res = y(end,1:6)-yf(1:6);
%res = [norm(y(end,1:3)-yf(1:3));norm(y(end,4:6)-yf(4:6))];
%res = [norm(y(end,1:2)-yf(1:2));y(end,3)-yf(3); norm(y(end,4:5)-yf(4:5));y(end,6)-yf(6);];
%res = sign(res).*log(sign(res).*res+1);
end