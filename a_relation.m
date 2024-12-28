%в этом скрипте перебираем различные соотношения полуосей

t_start = 0;
orbits='Flat_a';
planet_end = 'other';
N=1350;
m0=367;
eta=0.45;

step = 1/32;
ds = 3/4:step:12/2; %угловая дальность

step2 = 0.01;
da = 0.4:step2:2.5; %соотношение полуосей
L2 = length(ds);
L3 = length(da);

%время старта, угловое расстояние, сопряжённая переменная, семейство
PXevery_a=zeros(L2,L3, 8); %базис-вектор
Mevery_a=zeros(L2,L3);     %расход массы
Jevery_a=zeros(L2,L3);     %функционал
PHIevery_a = zeros(L2,L3); %фи в начальный момент
CONDevery_a = ones(L2,L3); %число обусловленности
Tevery_a = zeros(L2,L3);   %время перелёта
ANevery_a = zeros(L2,L3);  %угловая дальность
OMevery_a = zeros(L2,L3);  %разность фаз



%3 означает "никогда не проверялось"
%2 означает "метод не сошёлся"
%0 означает "метод сошёлся"
skipGrid_a=3*ones(L2,L3);
%%
% PXevery_a = [zeros(L2,70, 8), PXevery_a];
% Mevery_a = [zeros(L2,70), Mevery_a];
% Jevery_a = [zeros(L2,70), Jevery_a];
% PHIevery_a = [zeros(L2,70), PHIevery_a];
% CONDevery_a = [zeros(L2,70), CONDevery_a];
% Tevery_a = [zeros(L2,70), Tevery_a];
% ANevery_a = [zeros(L2,70), ANevery_a];
% OMevery_a = [zeros(L2,70), OMevery_a];
% skipGrid_a=[3*ones(L2,70), skipGrid_a];
%% Заполняем данные по Марсу (из plotUnited)
PXevery_a(:,43,:) = PXevery_a_MARS; %a = 1.52, Mars
PHIevery_a(:,43) = PHIevery_a_MARS;
OMevery_a(:,43) = OMevery_a_MARS;
skipGrid_a(:,43) = 0;
%%
display = 0;
for i=1:L2
    for j=113:113%L3
        [skip, i_nearest, j_nearest] = checkNear(skipGrid_a, i, j);
        if skip==0
            px_new=reshape(PXevery_a(i_nearest,j_nearest,:),[1,8]);
            if skip > 1 && norm(px_new) == 0
                continue
            end
            phi_new=PHIevery_a(i_nearest,j_nearest);
            omega=OMevery_a(i_nearest,j_nearest);
            disp([i,j])
            px=px_new;
            delta_s=ds(i)*2*pi;
            a_rel=da(j);
            modifier_p=10^(-4-sqrt(delta_s));
            x0_sec = [px delta_s/(2*pi) phi/(2*pi)];
            integration_acc=1e-12;
            calculate_condition=1;
            
            checkMethod_params = struct();
            checkMethod_params.t_start = t_start;
            checkMethod_params.delta_s = delta_s;
            checkMethod_params.rad = rad;
            checkMethod_params.UorR = 'u_hat';
            checkMethod_params.decreaseNonPsysical = decreaseNonPsysical;
            checkMethod_params.modifier_p = modifier_p;
            checkMethod_params.modifier_f = modifier_f;
            checkMethod_params.x0_sec = x0_sec;
            checkMethod_params.eta = eta;
            checkMethod_params.case_traj = case_traj;
            checkMethod_params.planet_start = planet_start;
            checkMethod_params.planet_end = planet_end;
            checkMethod_params.display = display;
            checkMethod_params.terminal_state = terminal_state;
            checkMethod_params.integration_acc = 1e-14;
            checkMethod_params.calculate_condition = calculate_condition;
            checkMethod_params.orbits = orbits;
            checkMethod_params.omega = omega;
            checkMethod_params.a_rel = a_rel;
            
            minimize_delta_omega = @(delta_omega) find_close_solution(delta_omega, checkMethod_params);

            try 
                delta_omega = fminbnd(minimize_delta_omega, -pi/32, pi/32);
                omega_min = delta_omega+omega;
                checkMethod_params.omega = omega_min;
                [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] = checkMethod(checkMethod_params);
            catch
                skipGrid_a(i,j)=2;
                continue
            end
            
            Jt_new = Jt(end)*eta;
            Jt_old = Jevery_a(i,j);
            % если новое решение хуже старого
            if skip==0 && Jt_new > Jt_old
                continue
            end

            if dr < 10000 && length(t)==1000
                skipGrid_a(i,j)=0;
            else
                skipGrid_a(i,j)=2;
                continue
            end
            m=massLP(Jt, m0, N);
            Mevery_a(i,j)=m(1)-m(end);
            CONDevery_a(i,j)=C;
            PXevery_a(i,j,:)=px;
            PHIevery_a(i,j)=phi;
            T_end = t(end)/(24*60*60); %в днях
            Tevery_a(i,j)=T_end;
            Jevery_a(i,j)=Jt(end)*eta;
            ANevery_a(i,j) = calculate_angular_distance(rr(:,1:3));

            OMevery_a(i,j)=omega_min;
        end
    end
end
%%
phi0 = 1.7981e-07;
mug = 1;

st.t = t_start;
st.planet = 'Earth';
st.mode = 'Flat';
st.delta_omega = omega_space(1,1);
st.a_rel = 1;

ae = 149597870700; %км
mug_0 = 132712.43994*(10^6)*(10^9); %км3с−2

[r0, v0]=planetModel(st);
eul = [pi/12 pi/4 pi/12];
rotmZYX = eul2rotm(eul);

R0 = [rotmZYX*r0'/ae; 0]*1e+03;
V0 = [rotmZYX*v0'/V_unit; 0]*1e+03;

PVevery_a=zeros(L2,L3, 3);
PRevery_a=zeros(L2,L3, 3);
for i = 1:L2
    for j = 1:L3
        phi0 = PHIevery_a(i,j);
        u0 = rToU(R0, phi0);
        w0 = vFromV(V0,R0,mug,phi0);
        ortdgduv = get_ortdgduv(u0,w0);
        dFdX = get_dgduv(u0,w0);
        vecPX1 = reshape(PXevery_a(i,j,:), [1,8]);
        vecPV_cartesian = a_reactive(u0,w0,vecPX1(1:4)',vecPX1(5:8)');
        b = vecPX1'-dFdX' * [0; 0; 0; vecPV_cartesian(1:3)]; %TODO переписать в символьных
        vecPR_cartesian = dFdX(1:3,:)'\b;

        vecPV_cartesian = rotmZYX^(-1)*vecPV_cartesian(1:3);%/(ae/sqrt(mug_0)).^2;
        vecPR_cartesian = rotmZYX^(-1)*vecPR_cartesian;%/(ae/sqrt(mug_0)).^2;
        PVevery_a(i,j,:) = vecPV_cartesian';
        PRevery_a(i,j,:) = vecPR_cartesian';
    end
end
%%

J_unit = 176.62545106129272;

figure(24);
Jevery_a_fix = Jevery_a/J_unit;
Jevery_a_fix(Jevery_a_fix==0)=nan;

s = surf(da, ANevery_a/(2*pi), Jevery_a_fix, 'FaceAlpha',FaceAlpha, 'HandleVisibility','off');
s.EdgeColor = 'none';
hold on;
Jevery_a_contour = [0.0001, 0.01];
contour3(repmat(da, L2,1), ANevery_a/(2*pi), Jevery_a_fix,Jevery_a_contour,...
    'ShowText','on', 'HandleVisibility','off', 'Color', 'k', 'LabelSpacing', 150);
Jevery_a_contour = [0.001, 0.03];
contour3(repmat(da, L2,1), ANevery_a/(2*pi), Jevery_a_fix,Jevery_a_contour,...
    'ShowText','on', 'HandleVisibility','off', 'Color', 'k', 'LabelSpacing', 300);
plot3([1,1],[0.75,3],[0.001,0.001], 'HandleVisibility','off', 'Color', 'k')
plot3([1,1],[3.75,6],[0.001,0.001], 'HandleVisibility','off', 'Color', 'k')
h = text(1, (3.5+3) /2,0.001, '0.0','HorizontalAlignment','center',...
    'VerticalAlignment','Bottom','FontSize',9); 
%plot3([1,1.001],[2,2.001],[0,0.001], 'DisplayName', 'Functional isolines', 'Color', 'k')
hold off;
%legend
zlabel('Значение функционала, безразм.')
xlabel('Radii ratio')
ylabel('Angular distance, revolutions')
%xlim([0.4,1])
xlim([1,2.5])
zlim([0,0.02])
%%
figure(25);

ae = 149597870.700; %км
mug = 132712.43994*(10^6); %км3с−2
Tunit_earth = (2*pi/sqrt(mug))*((ae).^(3/2))/(3600*24);
T_planets = repmat(Tunit_earth, L2,L3);
s = surf(repmat(da, L2,1), ANevery_a/(2*pi), Tevery_a./T_planets, 'FaceAlpha',FaceAlpha, 'HandleVisibility','off');
s.EdgeColor = 'none';
hold on;
%contour3(repmat(da, L2,1), ANevery_a/(2*pi), Tevery_a./T_planets, 'ShowText','on', 'DisplayName', 'Flight time isolines, years', 'Color', 'k');
contour3(repmat(da, L2,1), ANevery_a/(2*pi), Tevery_a./T_planets, 'ShowText','on', 'HandleVisibility','off', 'Color', 'k');
plot3([1,1.001],[2,2.001],[0,0.001], 'DisplayName', 'Flight time isolines, years', 'Color', 'k')
hold off;
legend
zlabel('Время перелёта, безразм.')
xlabel('Соотношение полуосей')
ylabel('Угловая дальность, безразм.')
%%
figure(26);

s = surf(PVevery_a(:,:,1), PVevery_a(:,:,2), ANevery_a/(2*pi), 'FaceAlpha',FaceAlpha, 'HandleVisibility','off');
s.EdgeColor = 'none';
hold on;
for j = 1:20:61
     pos_i = 5;
     plot3(PVevery_a(1:pos_i,j,1), PVevery_a(1:pos_i,j,2), ANevery_a(1:pos_i,j)/(2*pi), 'black', 'HandleVisibility','off')
     plot3(PVevery_a(pos_i+3:end,j,1), PVevery_a(pos_i+3:end,j,2), ANevery_a(pos_i+3:end,j)/(2*pi), 'black', 'HandleVisibility','off')
     text_label = num2str(da(j));
     h = text((PVevery_a(pos_i,j,1)+PVevery_a(pos_i+3,j,1))/2, (PVevery_a(pos_i,j,2)+PVevery_a(pos_i+3,j,2))/2, (ANevery_a(pos_i,j)+ANevery_a(pos_i+3,j))/(4*pi), text_label,'HorizontalAlignment','center',...
    'VerticalAlignment','Bottom','FontName','consolas','FontSize',11); 
     
end
Poses_h = [5,30,5,60,5,30];
Poses_h_i = 1;
for j = 61:40:L3
     pos_i = Poses_h(Poses_h_i);
     plot3(PVevery_a(1:pos_i,j,1), PVevery_a(1:pos_i,j,2), ANevery_a(1:pos_i,j)/(2*pi), 'black', 'HandleVisibility','off')
     plot3(PVevery_a(pos_i+3:end,j,1), PVevery_a(pos_i+3:end,j,2), ANevery_a(pos_i+3:end,j)/(2*pi), 'black', 'HandleVisibility','off')
     text_label = num2str(da(j));
     h = text((PVevery_a(pos_i,j,1)+PVevery_a(pos_i+3,j,1))/2, (PVevery_a(pos_i,j,2)+PVevery_a(pos_i+3,j,2))/2, (ANevery_a(pos_i,j)+ANevery_a(pos_i+3,j))/(4*pi), text_label,'HorizontalAlignment','center',...
    'VerticalAlignment','Bottom','FontName','consolas','FontSize',11); 
     Poses_h_i = Poses_h_i+1;
end
plot3(PVevery_a(pos_i+3:end,j,1), PVevery_a(pos_i+3:end,j,2), ANevery_a(pos_i+3:end,j)/(2*pi), 'black', 'DisplayName', 'Semimajor axis ratio isolines')
legend
hold off;
zlabel('Угловая дальность, безразм.')
xlabel('Pv x')
ylabel('Pv y')
%%
figure(27);

s = surf(PRevery_a(:,:,1), PRevery_a(:,:,2), ANevery_a/(2*pi), 'FaceAlpha',FaceAlpha, 'HandleVisibility','off');
s.EdgeColor = 'none';
hold on;
for j = 1:20:61
     pos_i = 5;
     plot3(PRevery_a(1:pos_i,j,1), PRevery_a(1:pos_i,j,2), ANevery_a(1:pos_i,j)/(2*pi), 'black', 'HandleVisibility','off')
     plot3(PRevery_a(pos_i+3:end,j,1), PRevery_a(pos_i+3:end,j,2), ANevery_a(pos_i+3:end,j)/(2*pi), 'black', 'HandleVisibility','off')
     text_label = num2str(da(j));
     h = text((PRevery_a(pos_i,j,1)+PRevery_a(pos_i+3,j,1))/2, (PRevery_a(pos_i,j,2)+PRevery_a(pos_i+3,j,2))/2, (ANevery_a(pos_i,j)+ANevery_a(pos_i+3,j))/(4*pi), text_label,'HorizontalAlignment','center',...
    'VerticalAlignment','Bottom','FontName','consolas','FontSize',11); 
end
Poses_h = [5,30,5,60,5,30];
Poses_h_i = 1;
for j = 61:40:L3
     pos_i = Poses_h(Poses_h_i);
     plot3(PRevery_a(1:pos_i,j,1), PRevery_a(1:pos_i,j,2), ANevery_a(1:pos_i,j)/(2*pi), 'black', 'HandleVisibility','off')
     plot3(PRevery_a(pos_i+3:end,j,1), PRevery_a(pos_i+3:end,j,2), ANevery_a(pos_i+3:end,j)/(2*pi), 'black', 'HandleVisibility','off')
     text_label = num2str(da(j));
     h = text((PRevery_a(pos_i,j,1)+PRevery_a(pos_i+3,j,1))/2, (PRevery_a(pos_i,j,2)+PRevery_a(pos_i+3,j,2))/2, (ANevery_a(pos_i,j)+ANevery_a(pos_i+3,j))/(4*pi), text_label,'HorizontalAlignment','center',...
    'VerticalAlignment','Bottom','FontName','consolas','FontSize',11); 
     Poses_h_i = Poses_h_i+1;
end
plot3(PRevery_a(pos_i+3:end,j,1), PRevery_a(pos_i+3:end,j,2), ANevery_a(pos_i+3:end,j)/(2*pi), 'black', 'DisplayName', 'Semimajor axis ratio isolines')
legend
hold off;
zlabel('Угловая дальность, безразм.')
xlabel('Pr x')
ylabel('Pr y')
%%

figure(24);
OMevery_a_fix = OMevery_a/(2*pi);
OMevery_a_fix(OMevery_a_fix==0)=nan;

for i = 100:L3
    for j = 2:L2
        if OMevery_a_fix(j,i) < OMevery_a_fix(j-1,i)
            OMevery_a_fix(j:end, i) = OMevery_a_fix(j:end, i) +1;
        end
    end
end

s = surf(da, ANevery_a/(2*pi), OMevery_a_fix, 'FaceAlpha',FaceAlpha, 'HandleVisibility','off');
s.EdgeColor = 'none';
hold on;
contour3(repmat(da, L2,1), ANevery_a/(2*pi), OMevery_a_fix,'ShowText','on', 'HandleVisibility','off','Color','k');
plot3([1,1.001],[2,2.001],[0,0.001], 'DisplayName', 'Phase difference isolines', 'Color', 'k')
%legend
hold off;
zlabel('Фазовая долгота, рад.')
xlabel('Соотношение радиусов')
ylabel('Угловая дальность, витки.')
view(0,90)
%ylim([0.75, 6])
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

function Jt_end = find_close_solution(delta_omega, checkMethod_params)
    checkMethod_params.omega = checkMethod_params.omega + delta_omega;
    [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] = checkMethod(checkMethod_params);
    Jt_end = Jt(end)*checkMethod_params.eta;
end