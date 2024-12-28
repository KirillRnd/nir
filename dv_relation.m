%в этом скрипте перебираем различные отлётные скорости

t_start = 0;
orbits='Flat';
planet_end = 'Mars';
N=1350;
m0=367;
eta=0.45;

step = 1/32;
ds = 3/4:step:12/2; %угловая дальность

dV_range = [0:10:50,100:100:4000]; %добавка к скорости в м/с
rad_dV_range = 0:pi/32:2*pi; %угол добавки к скорости в плоскости в радианах
L2 = length(ds);
L4 = length(dV_range);
L5 = length(rad_dV_range);

%время старта, угловое расстояние, сопряжённая переменная, семейство
PXevery_V=zeros(L2,L4,L5, 8); %базис-вектор
Mevery_V=zeros(L2,L4,L5);     %расход массы
Jevery_V=zeros(L2,L4,L5);     %функционал
PHIevery_V = zeros(L2,L4,L5); %фи в начальный момент
CONDevery_V = ones(L2,L4,L5); %число обусловленности
Tevery_V = zeros(L2,L4,L5);   %время перелёта
ANevery_V = zeros(L2,L4,L5);  %угловая дальность
OMevery_V = zeros(L2,L4,L5);  %разность фаз



%3 означает "никогда не проверялось"
%2 означает "метод не сошёлся"
%0 означает "метод сошёлся"
skipGrid_V=3*ones(L2,L4,L5);

%% Заполняем данные по Марсу (из plotUnited)
PXevery_V(:,1,1,:) = PXevery_a_MARS; %нулевая добавка скорости
PHIevery_V(:,1,1) = PHIevery_a_MARS;
OMevery_V(:,1,1) = OMevery_a_MARS;
skipGrid_V(:,1,1) = 0;
%%
for i=165:165 %угловая дальность
    for j=1:L4 %величина скорости
        for k=1:L5 %угол поворота скорости
            [skip, i_nearest, j_nearest, k_nearest] = checkNear_ijk(skipGrid_V, i, j, k);
            if skip>0
                px_new=reshape(PXevery_V(i_nearest,j_nearest,k_nearest,:),[1,8]);
                if skip > 1 && norm(px_new) == 0
                    continue
                end
                phi_new=PHIevery_V(i_nearest,j_nearest,k_nearest);
                omega=OMevery_V(i_nearest,j_nearest,k_nearest);
                disp([i,j,k])
                px=px_new;
                delta_s=ds(i)*2*pi;
                dV_value=dV_range(j);
                dV_rad=rad_dV_range(k);
                dV_add = dV_value*[cos(dV_rad), sin(dV_rad), 0];
                modifier_p=10^(-4-sqrt(delta_s));
                x0_sec = [px delta_s/(2*pi) phi/(2*pi)];
                integration_acc=1e-12;
                calculate_condition=1;
                
                checkMethod_params.t_start = t_start;
                checkMethod_params.delta_s = delta_s;
                checkMethod_params.rad = rad;
                checkMethod_params.UorR = UorR;
                checkMethod_params.decreaseNonPsysical = decreaseNonPsysical;
                checkMethod_params.modifier_p = modifier_p;
                checkMethod_params.modifier_f = modifier_f;
                checkMethod_params.x0_sec = x0_sec;
                checkMethod_params.eta = eta;
                checkMethod_params.case_traj = case_traj;
                checkMethod_params.planet_end = planet_end;
                checkMethod_params.display = display;
                checkMethod_params.terminal_state = terminal_state;
                checkMethod_params.integration_acc = integration_acc;
                checkMethod_params.calculate_condition = calculate_condition;
                checkMethod_params.orbits = orbits;
                checkMethod_params.omega = omega;
                checkMethod_params.a_rel = 1.52;
                checkMethod_params.dV_add = dV_add;
                
                minimize_delta_omega = @(delta_omega) find_close_solution(delta_omega, checkMethod_params);
    
                try 
                    delta_omega = fminbnd(minimize_delta_omega, -pi/16, pi/16);
                    omega_min = delta_omega+omega;
                    checkMethod_params.omega = omega_min;
                    [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] = checkMethod(checkMethod_params);
                catch
                    skipGrid_V(i,j,k)=2;
                    continue
                end
            
                if dr < 10000 && length(t)==1000
                    skipGrid_V(i,j,k)=0;
                else
                    skipGrid_V(i,j,k)=2;
                    continue
                end
                m=massLP(Jt, m0, N);
                Mevery_V(i,j,k)=m(1)-m(end);
                CONDevery_V(i,j,k)=C;
                PXevery_V(i,j,k,:)=px;
                PHIevery_V(i,j,k)=phi;
                T_end = t(end)/(24*60*60); %в днях
                Tevery_V(i,j, k)=T_end;
                Jevery_V(i,j,k)=Jt(end)*eta;
                ANevery_V(i,j,k) = calculate_angular_distance(rr(:,1:3));
    
                OMevery_V(i,j,k)=omega_min;
            end
        end
    end
end
%%
phi0 = 0;
mug = 1;

st.t = t_start;
st.planet = 'Earth';
st.mode = 'Flat';
st.delta_omega = omega_space(1,1);
st.a_rel = 1;

ae = 149597870700; %км
mug = 132712.43994*(10^6)*(10^9); %км3с−2

[r0, v0]=planetModel(st);
eul = [pi/12 pi/4 pi/12];
rotmZYX = eul2rotm(eul);

R0 = [rotmZYX*r0'/ae; 0]*1e+03;
V0 = [rotmZYX*v0'/V_unit; 0]*1e+03;
u0 = rToU(R0, phi0);
w0 = vFromV(V0,R0,mug,phi0);
ortdgduv = get_ortdgduv(u0,w0);
dFdX = get_dgduv(u0,w0);
PVevery_V=zeros(L2,L4,L5, 3);
PRevery_V=zeros(L2,L4,L5, 3);
for i = 1:L2
    for j = 1:L4
        for k = 1:L5
            vecPX1 = reshape(PXevery_V(i,j,k,:), [1,8]);
            vecPV_cartesian = a_reactive(u0,w0,vecPX1(1:4)',vecPX1(5:8)');
            b = vecPX1'-dFdX' * [0; 0; 0; vecPV_cartesian(1:3)]; %TODO переписать в символьных
            vecPR_cartesian = dFdX(1:3,:)'\b;
    
            vecPV_cartesian = rotmZYX^(-1)*vecPV_cartesian(1:3)/(ae/sqrt(mug_0)).^2;
            vecPR_cartesian = rotmZYX^(-1)*vecPR_cartesian/(ae/sqrt(mug_0)).^2;
            PVevery_V(i,j,k,:) = vecPV_cartesian';
            PRevery_V(i,j,k,:) = vecPR_cartesian';
        end
    end
end
%% выводим результат
figure(33);
Jevery_V_fix = Jevery_V;
Jevery_V_fix(Jevery_V_fix==0)=nan;

rad_dV_range_cos_sin = [cos(rad_dV_range); sin(rad_dV_range)];

dV_circle = repmat(reshape(dV_range, 1, 1, L4), 2, L5, 1).*repmat(rad_dV_range_cos_sin, 1, 1, L4);
X = dV_range'*cos(rad_dV_range)/1000;
Y = dV_range'*sin(rad_dV_range)/1000;
Z_slice = 1;
Z = reshape(Jevery_V_fix(Z_slice,:,:), L4, 65);
minimum = min(min(Z));
[i_idx,j_idx]=find(Z==minimum);
s = surf(X, Y, Z, 'FaceAlpha',FaceAlpha, 'HandleVisibility','off');
s.EdgeColor = 'none';
hold on;
contour3(X, Y, Z,linspace(0, 0.3,11), 'ShowText','on', 'HandleVisibility','off', 'LabelSpacing', 150);
contour3(X, Y, Z,linspace(0.5, 3,11), 'ShowText','on', 'HandleVisibility','off', 'LabelSpacing', 150);
% plot3(X(i_idx, j_idx),Y(i_idx, j_idx),Z(i_idx, j_idx), 'r*')
% %contour3(repmat(da, L2,1), ANevery_V/(2*pi), Jevery_V_fix,'ShowText','on', 'HandleVisibility','off');

%s = surf(scaler*reshape(PVevery_V(Z_slice,:,:,1), L4, 65),scaler*reshape(PVevery_V(Z_slice,:,:,2), L4, 65),reshape(Jevery_V_fix(Z_slice,:,:), L4, 65), 'FaceColor','green', 'FaceAlpha',FaceAlpha);
%s.EdgeColor = 'none';
%из dv_relation_best_angle
Jevery_Van_fix = Jevery_Van;
Jevery_Van_fix(Jevery_Van_fix==0)=nan;
dV_circle_best_angle = [dV_range.*cos(HIevery_Van(Z_slice,:)); dV_range.*sin(HIevery_Van(Z_slice,:))];
%plot3(dV_circle_best_angle(1,:),dV_circle_best_angle(2,:), Jevery_Van_fix(Z_slice,:), 'r','LineWidth',3)
scaler = 1500*1e+2;
% for i = 1:L4
%     xy_pv = scaler*PVevery_Van(Z_slice,i,1:2);
%     %plot3(X(i_idx, j_idx),Y(i_idx, j_idx),Z(i_idx, j_idx), 'r+')
%     plot3([dV_circle_best_angle(1,i), dV_circle_best_angle(1,i)+xy_pv(1)],...
%         [dV_circle_best_angle(2,i),dV_circle_best_angle(2,i)+xy_pv(2)],[Jevery_Van_fix(Z_slice,i), Jevery_Van_fix(Z_slice,i)], 'g');

%end

Jevery_Vhyp_fix = Jevery_Vhyp;
Jevery_Vhyp_fix(Jevery_Vhyp_fix==0)=nan;
dV_circle_best_angle_Vhyp = [dV_range.*cos(HIevery_Vhyp(Z_slice,:)); dV_range.*sin(HIevery_Vhyp(Z_slice,:))]/1000;
plot3(dV_circle_best_angle_Vhyp(1,:),dV_circle_best_angle_Vhyp(2,:), Jevery_Vhyp_fix(Z_slice,:), 'g','LineWidth',3,'DisplayName','Локальный минимум')
plot3(dV_circle_best_angle_Vhyp(1,32),dV_circle_best_angle_Vhyp(2,32), Jevery_Vhyp_fix(Z_slice,32)+0.1, 'r*','LineWidth',3,'DisplayName','Глобальный минимум')

dcm = datacursormode;
dcm.Enable = 'on';
set(dcm,'UpdateFcn',@datatipWithSubscript);
% scaler = 1500*1e+2;
% for i = 1:L4
%     xy_pv = scaler*PVevery_Vhyp(Z_slice,i,1:2);
%     %plot3(X(i_idx, j_idx),Y(i_idx, j_idx),Z(i_idx, j_idx), 'r+')
%     plot3([dV_circle_best_angle_Vhyp(1,i), dV_circle_best_angle_Vhyp(1,i)+xy_pv(1)],...
%         [dV_circle_best_angle_Vhyp(2,i),dV_circle_best_angle_Vhyp(2,i)+xy_pv(2)],[Jevery_Vhyp_fix(1,i), Jevery_Vhyp_fix(1,i)], 'g');
% 
% end
hold off;
axis equal;
legend('Location', 'best')
zlabel('Значение функционала, безразм.')
xlabel('X компонента отлётной скорости, км/с')
ylabel('Y компонента отлётной скорости, км/с')
%% выводим результат
figure(34);
Jevery_V_fix = Jevery_V;
Jevery_V_fix(Jevery_V_fix==0)=nan;


[X, Y] = meshgrid(rad_dV_range, dV_range);
s = surf(X, Y, reshape(Jevery_V_fix(1,:,:), 46, 65), 'FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';
hold on;
%contour3(repmat(da, L2,1), ANevery_V/(2*pi), Jevery_V_fix,'ShowText','on', 'HandleVisibility','off');
hold off;
zlabel('Значение функционала, безразм.')
ylabel('Величина скорости')
xlabel('Угол скорости')
%%
function [skip, i_nearest, j_nearest, k_nearest] = checkNear_ijk(arr, i, j, k)
    % checkNear
    s = size(arr);
    if arr(i,j,k)==0
        %уже есть решение
        skip = 0;
        i_nearest = i;
        j_nearest = j;
        k_nearest = k;
    else
        i_j_range = [i,i-1,i,i+1,i-1,i-1,j+1,j+1;
                    j-1,j,j+1,j,j-1,j+1,j+1,j-1];
        i_j_range_ext = [i,i,i-1,i,i+1,i-1,i-1,j+1,j+1;
                    j,j-1,j,j+1,j,j-1,j+1,j+1,j-1];
        i_j_k_range_1 = [i_j_range;repmat(k, 1, 8)];
        i_j_k_range_2 = [i_j_range_ext;repmat(k, 1, 9)+1];
        i_j_k_range_3 = [i_j_range_ext;repmat(k, 1, 9)-1];
        i_j_k_range = [i_j_k_range_1, i_j_k_range_2, i_j_k_range_3];
        idx_rearrange = [1, 2, 3, 4, 9, 18, 5, 6, 7, 8, ...
            10, 11, 12, 13, 19, 20, 21, 22,...
             14, 15, 16, 17, 23, 24, 25, 26]; %сначала перебираем ближайшие элементы
        i_j_k_range = i_j_k_range(:,idx_rearrange);
        for q=1:26
            i_nearest=i_j_k_range(1,q);
            j_nearest=i_j_k_range(2,q);
            k_nearest=i_j_k_range(3,q);
            if i_nearest>=1 && i_nearest <= s(1) && j_nearest>=1 && j_nearest <= s(2) && k_nearest>=1 && k_nearest <= s(3)
                if arr(i_nearest, j_nearest, k_nearest) == 0
                    skip=2;
                    return
                end
            end
        end
        %нет соседних решений
        skip = -1;
        i_nearest = 0;
        j_nearest = 0;
        k_nearest = 0;
    end
end

function Jt_end = find_close_solution(delta_omega, checkMethod_params)
    checkMethod_params.omega = checkMethod_params.omega + delta_omega;
    [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] = checkMethod(checkMethod_params);
    Jt_end = Jt(end)*checkMethod_params.eta;
end
function output_txt = datatipWithSubscript(obj,event_obj)
% Display the position of the data cursor with index/subscript information
%
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
%
% A simple datatip text update function that additionally displays the
% linear index or the subscripts of an object's selected point.
pos = get(event_obj,'Position');
% Get the linear index
dataIndex = get(event_obj,'DataIndex');
%output_txt = dataIndex;
output_txt = {['X: ',num2str(pos(1),4)],...
              ['Y: ',num2str(pos(2),4)],...
              [num2str(dataIndex)]};
end