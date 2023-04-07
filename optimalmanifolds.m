t_start = juliandate(2022,1,1);
t_end = juliandate(2022,12,31);
dt = linspace(t_start, t_end, 365);

N=1350;
m0=367;
eta=0.45;

step = 1/32;
ds = 3/4:step:12/2;
L1 = length(dt);
L2 = length(ds);

%время старта, угловое расстояние, сопряжённая переменная, семейство
PXevery=zeros(L1,L2, 8, 4);
Mevery=zeros(L1,L2, 4);
PHIevery = zeros(L1,L2, 4);
CONDevery = ones(L1,L2, 4);
%3 означает "никогда не проверялось"
%0 означает "метод сошёлся"
skipGrid=3*ones(L1,L2, 4);
[X,Y] = meshgrid(dt, ds);


PX_1 = PX(:,1:36);
PXevery(1,1:36,:,1)=PX_1';
PHIevery(1,1:36, 1)=PHI(1:36);
skipGrid(1,1:36, 1)=skipPoint(1:36);
skipGrid(skipGrid==-1)=0;
Mevery(1,1:36, 1)=M(1:36);
%F===family

PX_2 = PX(:,63:100);

PXevery(1,63:100,:,2)=PX_2';
PHIevery(1,63:100, 2)=PHI(63:100);
skipGrid(1,63:100, 2)=skipPoint(63:100);
skipGrid(skipGrid==-1)=0;
Mevery(1,63:100, 2)=M(63:100);
F=2;
%%
for i = 1:L1
    for j = 63:-1:40
        [skip, i_nearest, j_nearest] = checkNear(skipGrid, i, j, F);
        %[skip, px_new] = approximateNear(skipGrid, PXevery,dt, ds, i, j, F,step);
        if skip>1
            px_new=reshape(PXevery(i_nearest,j_nearest,:,F),[1,8]);
            if skip > 1 && norm(px_new) == 0
                continue
            end
            %неважен phi
            phi_new=PHIevery(i_nearest,j_nearest, F);
            %ds(j)
            disp([i,j])

            delta_s=ds(j)*2*pi;
            modifier_p=10^(-4-sqrt(delta_s));
            px=px_new;
            homotopy_steps = linspace(0,1,5);
            for k = 2:4
                delta_s=(ds(j_nearest)+(ds(j)-ds(j_nearest))*homotopy_steps(k))*2*pi;
                %modifier_p=10^(-4-sqrt(delta_s));
                integration_acc=1e-12;
                %Одиночный запуск метода и получение всех необходимых для графиков
                %переменных
                display = 1;
                terminal_state = 's';
                UorR = 'u_hat';
                
                x0_sec = [px delta_s/(2*pi) phi_new/(2*pi)];
                rad=0;
                calculate_condition=1;
    
                delta_t=dt(i_nearest)+(dt(i)-dt(i_nearest))*homotopy_steps(k);
                try 
                    [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] = checkMethod(delta_t,delta_s,rad,UorR,decreaseNonPsysical,modifier_p,modifier_f,x0_sec,eta, case_traj,planet_end, display,terminal_state,integration_acc,calculate_condition);
                    disp(C);
                catch
                    skipGrid(i,j, F)=2;
                    continue
                end
            end

            delta_s=ds(j)*2*pi;
            delta_t=dt(i);
            modifier_p=10^(-4-sqrt(delta_s));
            x0_sec = [px delta_s/(2*pi) phi/(2*pi)];
            integration_acc=1e-12;
            calculate_condition=1;

            try 
                [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] = checkMethod(delta_t,delta_s,rad,UorR,decreaseNonPsysical,modifier_p,modifier_f,x0_sec,eta, case_traj,planet_end, display,terminal_state,integration_acc,calculate_condition);
            catch
                skipGrid(i,j, F)=2;
                continue
            end
                %evaluation_time=evaluation_time+evaluation_time_2;
        
            if dr < 10000 && length(t)==1000
                skipGrid(i,j, F)=0;
            else
                skipGrid(i,j, F)=2;
                continue
            end
            m=massLP(Jt, m0, N);
            Mevery(i,j, F)=m(1)-m(end);
            CONDevery(i,j, F)=C;
            PXevery(i,j,:,F)=px;
            PHIevery(i,j, F)=phi;
        end
    end
end
%%

Mevery(Mevery==0)=nan;
ax = figure(11);
M1 = Mevery(1:L1, 1:100,1)';
cmap1 = colormap(ax, summer);
cmap2 = colormap(ax, cool);
C1 = rescale(M1,0,255);
mixed_cmap = cat(1,cmap1, cmap2);
colormap(ax,mixed_cmap);

Mscaled1 = 100*M1/m0;
s = surf(dt(1:L1)-dt(1),ds(1:100),Mscaled1,C1,'DisplayName','Семейство 1');
s.EdgeColor = 'none';
hold on;
C2 = rescale(M1,256,511);
M2 = Mevery(1:L1, 1:100,2)';
Mscaled2 = 100*M2/m0;
s = surf(dt(1:L1)-dt(1),ds(1:100),Mscaled2,C2,'DisplayName','Семейство 2');
s.EdgeColor = 'none';
%colormap(s,spring)
xlabel('Дата старта, дни (от 01.01.2022) ')
ylabel('Угловая дальность, витков')
zlabel('Затраты топлива, %')
hold off;
legend;
%%
CONDevery(CONDevery==0)=nan;
CONDevery(CONDevery>7000)=7000;
ax = figure(12);

cmap1 = colormap(ax,winter);
cmap2 = colormap(ax,summer);
mixed_cmap = cat(1,cmap1, cmap2);
colormap(ax,mixed_cmap);
COND1 = CONDevery(1:L1, 1:100,1)';
C1 = rescale(COND1,0,255);
s = surf(dt(1:L1)-dt(1),ds(1:100),COND1,C1);
s.EdgeColor = 'none';

hold on;
COND2=CONDevery(1:L1, 1:100,2)';
C2 = rescale(COND2,256,511);
s = surf(dt(1:L1)-dt(1),ds(1:100),COND2,C2);
s.EdgeColor = 'none';
xlabel('Time, days')
ylabel('Angle distance')
zlabel('Condition number')
hold off;

%% вычисление чисел обусловленности на прямоугольнике
t_point = 222;
s_point = 100;
vecPX1 = reshape(PXevery(t_point,s_point,:,1), [1,8]);
vecPX2 = reshape(PXevery(t_point,s_point,:,2), [1,8]);
N_points_cond = 10;
CONDparallelogram = zeros(N_points_cond,N_points_cond);
CONDparallelogram_points = linspace(0,1,N_points_cond);
for i = 1:N_points_cond
    for j = 1:N_points_cond
        px_mixed = vecPX1*CONDparallelogram_points(i)+vecPX2*CONDparallelogram_points(j);
        C = calculateCondOnly(dt(t_point), px_mixed, ds(s_point), PHIevery(t_point,s_point));
        CONDparallelogram(i,j)=C;
    end
end
figure(13);
surf(CONDparallelogram_points,CONDparallelogram_points,CONDparallelogram);
xlabel('1 family')
ylabel('2 family')

%%
cosineDistEvery=zeros(L1,L2);
cosineSimEvery=zeros(L1,L2);
for i = 1:L1
    for j = 1:L2
        vecPX1 = reshape(PXevery(i,j,:,1), [1,8]);
        %vecPX2 = reshape(PXevery(1,41,:,1), [1,8]);
        vecPX2 = reshape(PXevery(i,j,:,2), [1,8]);
        cosSim = vecPX1*vecPX2'/norm(vecPX1)/norm(vecPX2);
        %cosDist= 1-abs(cosSim);
        cosDist=norm(vecPX1)/norm(vecPX2);
        cosineDistEvery(i,j)=cosDist;
        cosineSimEvery(i,j)=cosSim;
    end
end
cosineDistEvery(cosineDistEvery==0)=NaN;
cosineSimEvery(cosineSimEvery==0)=NaN;
figure(14);
s = surf(dt(1:L1)-dt(1),ds(41:100),cosineSimEvery(:, 41:end)');
s.EdgeColor = 'none';
xlabel('Дата старта, дни (от 01.01.2022)')
ylabel('Угловая дальность, витков')
zlabel('Косинусное сходство')
%%
function [skip, i_nearest, j_nearest] = checkNear(arr, i, j, F)
    % checkNear
    s = size(arr);
    if arr(i,j,F)==0
        %уже есть решение
        skip = 1;
        i_nearest = i;
        j_nearest = j;
    else
%         i_range = i-1:i+1;
%         j_range = j-1:j+1;
        
        i_j_range = [i,i-1,i,i+1,i-1,i-1,j+1,j+1;
                    j-1,j,j+1,j,j-1,j+1,j+1,j-1];
%         i_j_range = [i-1,i,i+1,i,i-1,i-1,j+1,j+1;
%                      j,j-1,j,j+1,j-1,j+1,j+1,j-1];
        for k=1:8
            i_nearest=i_j_range(1,k);
            j_nearest=i_j_range(2,k);
            if i_nearest>=1 && i_nearest <= s(1) && j_nearest>=1 && j_nearest <= s(2)
                if arr(i_nearest, j_nearest, F) == 0
                    skip=2;
                    return
                end
            end
        end
%         for i_nearest=i_range(i_range >= 1 & i_range <= s(1))
%             for j_nearest=j_range(j_range >= 1 & j_range <= s(2))
%                 if arr(i_nearest, j_nearest) == 0
%                     skip=2;
%                     return
%                 end
%             end
%         end
        %нет соседних решений
        skip = 0;
        i_nearest = 0;
        j_nearest = 0;
    end
end

function [skip, px_approx] = approximateNear(skipGrid, PXevery,dt, ds, i, j, F,step)

    px_approx=zeros([1,8]);
    if skipGrid(i,j)==0
        %уже есть решение
        skip=1;
        px_approx=PXevery(i,i,:,F);
        return
    end
    R = 3;
    dt_local=dt(i)-R<=dt & dt<=dt(i)+R;
    ds_local=ds(j)-R*step<=ds & ds<=ds(j)+R*step;
    localskip=skipGrid(dt_local,ds_local,F);
    if sum(sum(localskip==0))<6
        %нет возможности решить
        skip=0;
        return
    end
    localPX = PXevery(dt_local,ds_local,:,F);
    DTlocal=dt(dt_local)-dt(1);
    DSlocal=ds(ds_local);
    [DTmesh, DSmesh] = meshgrid(DTlocal,DSlocal);
    DT_points = DTmesh(localskip==0);
    DS_points = DSmesh(localskip==0);
    for k=1:8
        PX_layer = localPX(:,:,k);
        PX_layer_points = PX_layer(localskip==0);
        approx_fun = fit([DT_points, DS_points], PX_layer_points, 'poly22');
        px_approx(k)=approx_fun(dt(i)-dt(1),ds(j));
    end
    skip=2;
end