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
L2 = length(ds);
L4 = length(dV_range);

%время старта, угловое расстояние, сопряжённая переменная, семейство
PXevery_Van=zeros(L2,L4, 8); %базис-вектор
Mevery_Van=zeros(L2,L4);     %расход массы
Jevery_Van=zeros(L2,L4);     %функционал
PHIevery_Van = zeros(L2,L4); %фи в начальный момент
CONDevery_Van = ones(L2,L4); %число обусловленности
Tevery_Van = zeros(L2,L4);   %время перелёта
ANevery_Van = zeros(L2,L4);  %угловая дальность
OMevery_Van = zeros(L2,L4);  %разность фаз
HIevery_Van = zeros(L2,L4);  %оптимальное направление скорости

display = 0;

%3 означает "никогда не проверялось"
%2 означает "метод не сошёлся"
%0 означает "метод сошёлся"
skipGrid_V=3*ones(L2,L4);
UorR = 'u_hat';
%% Заполняем данные по Марсу (из plotUnited)
PXevery_Van(:,1,:) = PXevery_a_MARS; %нулевая добавка скорости
PHIevery_Van(:,1) = PHIevery_a_MARS;
OMevery_Van(:,1) = OMevery_a_MARS;
skipGrid_V(:,1) = 0;
HIevery_Van(:,1) = 1.4322; %надо добыть из pv
%%
for i=1:1 %угловая дальность
    for j=1:L4 %величина скорости
        [skip, i_nearest, j_nearest] = checkNear(skipGrid_V, i, j);
        if skip==0
            px_new=reshape(PXevery_Van(i_nearest,j_nearest,:),[1,8]);
            if skip > 1 && norm(px_new) == 0
                continue
            end
            phi_new=PHIevery_Van(i_nearest,j_nearest);
            omega=OMevery_Van(i_nearest,j_nearest); %оптимальная разность фаз
            hi=HIevery_Van(i_nearest,j_nearest);    %оптимальное направление
            disp([i,j])
            px=px_new;
            delta_s=ds(i)*2*pi;
            dV_value=dV_range(j);
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
            checkMethod_params.dV_value = dV_value;
            
            minimize_delta_omega_hi = @(delta) find_close_solution_omega_hi(delta(1), delta(2), checkMethod_params, omega, hi);
            x0 = [0;0];
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            nonlcon=[];
            lb = -pi*[1/16, 1/16];
            ub = pi*[1/16, 1/16];

            try
                omega_hi = fmincon(minimize_delta_omega_hi, x0, A, b, Aeq, beq, lb, ub, nonlcon);
                delta_omega = omega_hi(1);
                delta_hi    = omega_hi(2);
                omega_min = delta_omega+omega;
                hi_min = delta_hi+hi;
                checkMethod_params.omega = omega_min;
                checkMethod_params.hi    = hi_min;
                [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] = checkMethod(checkMethod_params);
            catch
                skipGrid_V(i,j)=2;
                continue
            end
        
            if dr < 10000 && length(t)==1000
                skipGrid_V(i,j)=0;
            else
                skipGrid_V(i,j)=2;
                continue
            end
            m=massLP(Jt, m0, N);
            Mevery_Van(i,j)=m(1)-m(end);
            CONDevery_Van(i,j)=C;
            PXevery_Van(i,j,:)=px;
            PHIevery_Van(i,j)=phi;
            T_end = t(end)/(24*60*60); %в днях
            Tevery_Van(i,j)=T_end;
            Jevery_Van(i,j)=Jt(end)*eta;
            ANevery_Van(i,j) = calculate_angular_distance(rr(:,1:3));

            OMevery_Van(i,j)=omega_min;
            HIevery_Van(i,j)=hi_min;
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

PVevery_Van=zeros(L2,L4, 3);
PRevery_Van=zeros(L2,L4, 3);
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
        vecPX1 = reshape(PXevery_Van(i,j,:), [1,8]);
        vecPV_cartesian = a_reactive(u0,w0,vecPX1(1:4)',vecPX1(5:8)');
        b = vecPX1'-dFdX' * [0; 0; 0; vecPV_cartesian(1:3)]; %TODO переписать в символьных
        vecPR_cartesian = dFdX(1:3,:)'\b;

        vecPV_cartesian = rotmZYX^(-1)*vecPV_cartesian(1:3);
        vecPR_cartesian = rotmZYX^(-1)*vecPR_cartesian;%/(ae/sqrt(mug_0)).^2;
        PVevery_Van(i,j,:) = vecPV_cartesian';
        PRevery_Van(i,j,:) = vecPR_cartesian';
    end
end
%% выводим результат

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
        idx_rearrange = [1, 2, 3, 4, 5, 6, 7, 8];
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

function Jt_end = find_close_solution_omega_hi(delta_omega, delta_hi, checkMethod_params, omega, hi)
    omega_min = delta_omega+omega;
    hi_min = delta_hi+hi;
    checkMethod_params.omega = omega_min;
    checkMethod_params.hi    = hi_min;
    [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] = checkMethod(checkMethod_params);
    Jt_end = Jt(end)*checkMethod_params.eta;
%     if dr > 10000
%         Jt_end=Jt_end*1000;
%     end
end