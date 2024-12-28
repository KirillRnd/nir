% %в этом скрипте перебираем различные отлётные скорости
% 
% t_start = 0;
% orbits='Flat';
% planet_end = 'Mars';
% N=1350;
% m0=367;
% eta=0.45;
% 
% step = 1/32;
% ds = 3/4:step:12/2; %угловая дальность
% 
% dV_range = [0:10:50,100:100:4000]; %добавка к скорости в м/с
% L2 = length(ds);
% L4 = length(dV_range);
% 
% %время старта, угловое расстояние, сопряжённая переменная, семейство
% PXevery_Vhyp_end=zeros(L2,L4, 8); %базис-вектор
% Mevery_Vhyp_end=zeros(L2,L4);     %расход массы
% Jevery_Vhyp_end=zeros(L2,L4);     %функционал
% PHIevery_Vhyp_end = zeros(L2,L4); %фи в начальный момент
% CONDevery_Vhyp_end = ones(L2,L4); %число обусловленности
% Tevery_Vhyp_end = zeros(L2,L4);   %время перелёта
% ANevery_Vhyp_end = zeros(L2,L4);  %угловая дальность
% OMevery_Vhyp_end = zeros(L2,L4);  %разность фаз
% HIevery_Vhyp_end = zeros(L2,L4);  %оптимальное направление скорости
% 
% display = 1;
% 
% %3 означает "никогда не проверялось"
% %2 означает "метод не сошёлся"
% %0 означает "метод сошёлся"
% skipGrid_Vhyp_end=3*ones(L2,L4);
% UorR = 'u_hat_v_hyp_end';
% %UorR = 'u_hat';
% %% Заполняем данные по Марсу (из plotUnited)
% PXevery_Vhyp_end(:,1,:) = PXevery_a_MARS; %нулевая добавка скорости
% PHIevery_Vhyp_end(:,1) = PHIevery_a_MARS;
% OMevery_Vhyp_end(:,1) = OMevery_a_MARS;
% skipGrid_Vhyp_end(:,1) = 0;
% HIevery_Vhyp_end(:,1) = 1.4322; %надо добыть из pv
%%
%load('D:\MATLAB\ipm\mat-files\04-Nov-2024.mat')
display = 1; 
for i=8:8%105%L2 %угловая дальность
    for j= 30:30%:L4%L4 %величина скорости
        [skip, i_nearest, j_nearest] = checkNear(skipGrid_Vhyp_end, i, j);
        if skip>0
            px_new=reshape(PXevery_Vhyp_end(i_nearest,j_nearest,:),[1,8]);
            if skip > 1 && norm(px_new) == 0
                continue
            end
            phi_new=PHIevery_Vhyp_end(i_nearest,j_nearest);
            omega=OMevery_Vhyp_end(i_nearest,j_nearest); %оптимальная разность фаз
            hi=HIevery_Vhyp_end(i_nearest,j_nearest);    %оптимальное направление
            disp([i,j])
            i_last = i;
            j_last = j;
            px=px_new;
            phi=phi_new;
            delta_s=ds(i)*2*pi;
            dV_value_end=dV_range(j);
            Jt_old = Jevery_Vhyp_end(i,j);
            %dV_add = dV_value*[cos(dV_rad), sin(dV_rad), 0];
            modifier_p=10^(-4-sqrt(delta_s));
            x0_sec = [px delta_s/(2*pi) phi/(2*pi)];
            integration_acc=1e-12;
            calculate_condition=1;

            checkMethod_params = struct();
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
            checkMethod_params.hi = hi;
            checkMethod_params.a_rel = 1.52;
            %checkMethod_params.dV_add = dV_add;
            %checkMethod_params.dV_value = dV_value;
            checkMethod_params.dV_value_end = dV_value_end;
            reg_p = 0.0;
            minimize_delta_omega = @(delta_omega) find_close_solution(delta_omega, checkMethod_params, omega, reg_p);
    

            try
                angle_range = 2*pi/16;
                % if (j>=20) && (j<=27) 
                %     reg_p = 0.001;
                %     angle_range = 3*pi/16;
                % end

                [delta_omega, Jt_end_min] = fminbnd(minimize_delta_omega, -angle_range, angle_range);
                omega_min = delta_omega+omega;
                checkMethod_params.omega = omega_min;
                checkMethod_results = checkMethod_ver2(checkMethod_params);
                %[dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] 
            catch
                if skip>0
                    skipGrid_Vhyp_end(i,j)=2;
                end
                
                continue
            end
            Jt_new = checkMethod_results.Jt(end)*eta;
            % если новое решение хуже старого
            if skip==0 && Jt_new > Jt_old
                continue
            end

            if checkMethod_results.dr < 10000 && length(checkMethod_results.t)==1000
                m=massLP(checkMethod_results.Jt, m0, N);
                Mevery_Vhyp_end(i,j)=m(1)-m(end);
                CONDevery_Vhyp_end(i,j)=checkMethod_results.C;
                PXevery_Vhyp_end(i,j,:)=checkMethod_results.px;
                PHIevery_Vhyp_end(i,j)=checkMethod_results.phi;
                T_end = checkMethod_results.t(end)/(24*60*60); %в днях
                Tevery_Vhyp_end(i,j)=T_end;
                Jevery_Vhyp_end(i,j)=checkMethod_results.Jt(end)*eta;
                ANevery_Vhyp_end(i,j) = calculate_angular_distance(checkMethod_results.rr(:,1:3));
    
                OMevery_Vhyp_end(i,j)=omega_min;
                HIevery_Vhyp_end(i,j)=checkMethod_results.hi_opt;
                skipGrid_Vhyp_end(i,j)=0; %set status only in the end
                savefilename = join(['mat-files/',date]);
                save(savefilename);
                disp('OK')
            else
                skipGrid_Vhyp_end(i,j)=2;
                continue
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
%mug = 132712.43994*(10^6)*(10^9); %км3с−2

[r0, v0]=planetModel(st);
eul = [pi/12 pi/4 pi/12];
rotmZYX = eul2rotm(eul);

R0 = [rotmZYX*r0'/ae; 0]*1e+03;

PVevery_Vhyp_end=zeros(L2,L4, 3);
PRevery_Vhyp_end=zeros(L2,L4, 3);
for i = 1:L2
    for j = 1:L4
        hi = HIevery_Vhyp_end(i,j);
        dV_value_end=dV_range(j);
        V0 = [rotmZYX*v0'/V_unit; 0]*1e+03;
        phi0 = PHIevery_Vhyp_end(i,j);
        u0 = rToU(R0, phi0);
        w0 = vFromV(V0,R0,mug,phi0);
        ortdgduv = get_ortdgduv(u0,w0);
        dFdX = get_dgduv(u0,w0);
        vecPX1 = reshape(PXevery_Vhyp_end(i,j,:), [1,8]);
        vecPV_cartesian = a_reactive(u0,w0,vecPX1(1:4)',vecPX1(5:8)');
        b = vecPX1'-dFdX' * [0; 0; 0; vecPV_cartesian(1:3)]; %TODO переписать в символьных
        vecPR_cartesian = dFdX(1:3,:)'\b;

        vecPV_cartesian = rotmZYX^(-1)*vecPV_cartesian(1:3);
        vecPR_cartesian = rotmZYX^(-1)*vecPR_cartesian;
        PVevery_Vhyp_end(i,j,:) = vecPV_cartesian';
        PRevery_Vhyp_end(i,j,:) = vecPR_cartesian';
    end
end
%% выводим результат

figure(35);
Jevery_Vhyp_end_fix = Jevery_Vhyp_end;%PXevery_Vhyp_end(:,:,1);
%Jevery_Vhyp_end_fix = OMevery_Vhyp_end;
%Jevery_Vhyp_end_fix = PXevery_Vhyp_end(:,:,1);
Jevery_Vhyp_end_fix(Jevery_Vhyp_end_fix==0)=nan;
Jevery_Vhyp_end_fix(skipGrid_Vhyp_end>0)=nan;
%Jevery_Vhyp_end_fix(11:end,:)=nan;

FaceAlpha = 0.4;

[X, Y] = meshgrid(ds, dV_range);
s = surf(X, Y, Jevery_Vhyp_end_fix', 'FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';
hold on;
%contour3(repmat(da, L2,1), ANevery_V/(2*pi), Jevery_V_fix,'ShowText','on', 'HandleVisibility','off');
xlabel('Угловая дальность')
ylabel('Величина скорости')
zlabel('Значение функционала, безразм.')

contour3(X, Y, Jevery_Vhyp_end_fix',...
    'ShowText','on', 'HandleVisibility','off', 'Color', 'k');
hold off;

dcm = datacursormode;
dcm.Enable = 'on';
dcm.UpdateFcn = @(obj,event_obj)datatipWithSubscript(obj,event_obj);
xlim([0.75,4])
%%
figure(35);
[X, Y] = meshgrid(ds, dV_range);
PXevery_Vhyp_end_fix = PXevery_Vhyp_end;
%PXevery_Vhyp_end_fix(PXevery_Vhyp_end_fix==0)=nan;
X_fit = reshape(X,[L2*L4, 1]);
Y_fit = reshape(Y,[L2*L4, 1]);
Z_fit = reshape(PXevery_Vhyp_end_fix(:,:,1)',[L2*L4, 1]);
X_fit(Z_fit==0) = [];
Y_fit(Z_fit==0) = [];
Z_fit(Z_fit==0) = [];
PX_fpoly = fit([X_fit Y_fit],Z_fit,'poly55');
s = surf(X, Y,PXevery_Vhyp_end_fix(:,:,1)', 'FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';
hold on;
%surf(X, Y,PX_fpoly(X, Y), 'FaceAlpha',0);
hold off;
%%
function [skip, i_nearest, j_nearest] = checkNear(arr, i, j)
    % checkNear
    s = size(arr);
    if arr(i,j)==0
        %уже есть решение
        skip = 0;
        i_nearest = i;
        j_nearest = j;
    else
        i_j_range = [i,i-1,i,i+1,i-1,i-1,j+1,j+1;
                    j-1,j,j+1,j,j-1,j+1,j+1,j-1];
        idx_rearrange = [2, 4, 1, 3, 5, 6, 7, 8];
        %idx_rearrange = [4, 2, 3, 1, 5, 6, 7, 8];
        %idx_rearrange = randperm(8);
        i_j_range = i_j_range(:,idx_rearrange);
        for k=1:8
            i_nearest=i_j_range(1,k);
            j_nearest=i_j_range(2,k);
            if i_nearest>=1 && i_nearest <= s(1) && j_nearest>=1 && j_nearest <= s(2)
                if arr(i_nearest, j_nearest) == 0
                    skip=2;
                    return
                end
            end
        end
        %нет соседних решений
        skip = -1;
        i_nearest = 0;
        j_nearest = 0;
    end
end

function J_mix = find_close_solution(delta_omega, checkMethod_params, omega, reg_p)
    omega_min = delta_omega+omega;
    checkMethod_params_COPY =checkMethod_params;
    checkMethod_params_COPY.omega = omega_min;
    checkMethod_results = checkMethod_ver2(checkMethod_params_COPY);
    Jt_end = checkMethod_results.Jt(end)*checkMethod_params_COPY.eta;
    J_mix = reg_p*(log10(checkMethod_results.dr)/11.1749+log10(checkMethod_results.dv)/4.4740)+Jt_end;
    J_mix = J_mix*100;
    %J_mix=Jt_end;
    % if checkMethod_results.dr > 10000
    %     Jt_end=Jt_end*10;
    % end
end
