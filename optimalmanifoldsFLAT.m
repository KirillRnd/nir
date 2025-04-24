%t_start = juliandate(2022,1,1);
%t_end = juliandate(2022,12,31);
%dt = linspace(t_start, t_end, 365);

t_start = 0;
omega_space = linspace(-pi, pi, 361); %шаг 1 градус

orbits='Flat';
N=1350;
m0=367;
eta=0.45;

step = 1/32;
ds = 3/4:step:12/2;
L1 = length(omega_space);
L2 = length(ds);

%время старта, угловое расстояние, сопряжённая переменная, семейство
PXevery=zeros(L1,L2, 8, 4);
Mevery=zeros(L1,L2, 4);
Jevery=zeros(L1,L2, 4);
PHIevery = zeros(L1,L2, 4);
CONDevery = ones(L1,L2, 4);
Tevery = zeros(L1,L2, 4);
%3 означает "никогда не проверялось"
%0 означает "метод сошёлся"
skipGrid=3*ones(L1,L2, 4);
[X,Y] = meshgrid(omega_space, ds);


PX_1 = PX(:,1:45);
PXevery(1,1:45,:,1)=PX_1';
PHIevery(1,1:45, 1)=PHI(1:45);
skipGrid(1,1:45, 1)=skipPoint(1:45);
skipGrid(skipGrid==-1)=0;
Mevery(1,1:45, 1)=M(1:45);
%F===family

PX_2 = PX(:,50:90);

PXevery(1,50:90,:,2)=PX_2';
PHIevery(1,50:90, 2)=PHI(50:90);
skipGrid(1,50:90, 2)=skipPoint(50:90);
skipGrid(skipGrid==-1)=0;
Mevery(1,50:90, 2)=M(50:90);

PX_3 = PX(:,143:162);

PXevery(1,143:162,:,3)=PX_3';
PHIevery(1,143:162, 3)=PHI(143:162);
skipGrid(1,143:162, 3)=skipPoint(143:162);
skipGrid(skipGrid==-1)=0;
Mevery(1,143:162, 3)=M(143:162);

PX_4 = PX(:,163:169);

PXevery(1,163:169,:,4)=PX_4';
PHIevery(1,163:169, 4)=PHI(163:169);
skipGrid(1,163:169, 4)=skipPoint(163:169);
skipGrid(skipGrid==-1)=0;
Mevery(1,163:169, 4)=M(163:169);
%%
Mevery(skipGrid>0)=0;
Jevery(skipGrid>0)=0;
%%
F=5;
apply_homotopy = true;
for i = L1:-1:180%L1:-1:350
    for j = 1:L2%50:L2
        [skip, i_nearest, j_nearest] = checkNear(skipGrid, i, j, F);
        if skip>0
            px_new=reshape(PXevery(i_nearest,j_nearest,:,F),[1,8]);
            if skip > 1 && norm(px_new) == 0
                continue
            end
            %неважен phi
            phi_new=PHIevery(i_nearest,j_nearest, F);
            %ds(j)
            disp([i,j])

            px=px_new;
            homotopy_fail = false;
            if apply_homotopy
                delta_s=ds(j)*2*pi;
                modifier_p=10^(-4-sqrt(delta_s));
                homotopy_steps = linspace(0,1,5);
                for k = 2:4%2:4
                    if homotopy_fail
                        continue;
                    end
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
                    omega=omega_space(i_nearest)+(omega_space(i)-omega_space(i_nearest))*homotopy_steps(k);
                    %delta_t=dt(i_nearest)+(dt(i)-dt(i_nearest))*homotopy_steps(k);
                    try 
                        [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] =...
                            checkMethod(t_start,delta_s,rad,UorR,decreaseNonPsysical,modifier_p,modifier_f,x0_sec,eta, case_traj,planet_end, display,terminal_state,integration_acc,calculate_condition, orbits, omega);
                        disp(C);
                        if k > 2
                            disp(['Величина шага продолжения: ', num2str(norm(rr_old - rr(end,:)),'%10.2e\n'), ' а.е.'])
                        end
                        rr_old = rr(end,:);
                    catch
                        skipGrid(i,j, F)=2; 
                        continue
                    end
                    if dr >= 10000 || length(t)~=1000
                        homotopy_fail = true;
                    end
                end
            end
            if homotopy_fail
                skipGrid(i,j, F)=2;
                continue
            end
            delta_s=ds(j)*2*pi;
            omega=omega_space(i);
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

            try 
                [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] = checkMethod(checkMethod_params);
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
            T_end = t(end)/(24*60*60); %в днях
            Tevery(i,j, F)=T_end;
            Jevery(i,j, F)=Jt(end)*eta;
        end
    end
end
%% заполнить время И угловое расстояние
ANevery = zeros(361, 169, 5);
for F=1:5
    for i = 1:L1
        for j = 1:L2
            [skip, i_nearest, j_nearest] = checkNear(skipGrid, i, j, F);
            if skip <= 0 && i_nearest>0 && j_nearest>0
                px_new=reshape(PXevery(i_nearest,j_nearest,:,F),[1,8]);
                %неважен phi
                phi_new=PHIevery(i_nearest,j_nearest, F);
                %ds(j)
                disp([i,j])
                px=px_new;
    
                delta_s=ds(j)*2*pi;
                omega=omega_space(i);
                modifier_p=10^(-4-sqrt(delta_s));
                x0_sec = [px delta_s/(2*pi) phi/(2*pi)];
                integration_acc=1e-12;
                calculate_condition=0;
    
                [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] = calculateTimeKS(t_start, ...
                    delta_s,rad,UorR,decreaseNonPsysical,modifier_p,modifier_f,x0_sec,eta, case_traj,planet_end, ...
                    display,terminal_state,integration_acc,calculate_condition, orbits, omega);
                
                ANevery(i,j, F) = calculate_angular_distance(rr(:,1:3));
                T_end = t(end)/(24*60*60); %в днях
                Tevery(i,j, F)=T_end;
                Jevery(i,j, F)=Jt(end)*eta;
            end
        end
    end
end
%% одна поверхность
J_unit =  176.6255; %есть в plotUnited.m
Jevery_fix = Jevery/J_unit;
Jevery_fix(Jevery_fix==0)=nan;

figure(11);
J1 = Jevery_fix(1:L1, 1:L2,1)';
drop1Fam_j = 1:10;
drop1Fam_i = 181:L1;
J1(drop1Fam_j,drop1Fam_i) = nan;
for i=181:L1
    for j=10:20
        if j <= 10+(i-181)*7/110
            J1(j,i) = nan;
        end
    end
end

X_axis = 180*omega_space(1:L1)/pi;
Y_axis = ds(1:L2);
FaceAlpha=0.4;
s = surf(X_axis,Y_axis,J1-0.01,'DisplayName','Семейство 1', 'FaceColor','cyan','FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';
hold on;
J1_isolines = round(logspace(-4,-0.7328,101)/J_unit,5);
s_contour = contour3(X_axis,Y_axis,J1,J1_isolines(2:2:end),'ShowText','on', 'Color', 'k', 'LabelSpacing', 450);
J1_isolines = round(logspace(-0.5850,3,11)/J_unit,5);
s_contour = contour3(X_axis,Y_axis,J1,J1_isolines(1:2:end),'ShowText','on', 'Color', 'k', 'LabelSpacing', 300);
plot3(X_axis(1),Y_axis(1),J1(1,1),'r*')
hold off;
ylabel('Angular distance, revolutions')
xlabel('Initial synodic phase')
%set(gca, 'ZScale', 'log')
%set(gca,'ColorScale','log')
grid;
view(0,90)
xticks([X_axis(1) (X_axis(1)+X_axis(end))/2 X_axis(end)])
xticklabels({'-\pi','0','\pi'})
%colorbar;
%% оптимальные поверхности
Mevery_fix = Mevery;
Mevery_fix(Mevery_fix==0)=nan;
figure(11);
M1 = Mevery_fix(1:L1, 1:L2,1)';
M2 = Mevery_fix(1:L1, 1:L2,2)';
M3 = Mevery_fix(1:L1, 1:L2,3)';
M4 = Mevery_fix(1:L1, 1:L2,4)';
M5 = Mevery_fix(1:L1, 1:L2,5)';


drop1Fam_j = 1:10;
drop1Fam_i = 181:L1;
M1(drop1Fam_j,drop1Fam_i) = nan;
for i=181:L1
    for j=10:20
        if j <= 10+(i-181)*7/110
            M1(j,i) = nan;
        end
    end
end

drop2Fam_j = 10:39;
drop2Fam_i = 161:L1;
M2(drop2Fam_j,drop2Fam_i) = nan;

drop3Fam_j = 10:71;
drop3Fam_i = 160:L1;
M3(drop3Fam_j,drop3Fam_i) = nan;

X_axis = 180*omega_space(1:L1)/pi;
Y_axis = ds(1:L2);

Mscaled5 = 100*M5/m0;
%Mscaled1 = M1;

FaceAlpha=0.4;
s = surf(180*omega_space(1:L1)/pi,ds(1:L2),Mscaled5,'DisplayName','Family 0', 'FaceColor','blue','FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';

hold on;

Tevery_fix = Tevery;
Tevery_fix(Tevery_fix==0)=nan;
T1 = Tevery_fix(1:L1, 1:L2,1)';
T2 = Tevery_fix(1:L1, 1:L2,2)';
T3 = Tevery_fix(1:L1, 1:L2,3)';
T4 = Tevery_fix(1:L1, 1:L2,4)';
T5 = Tevery_fix(1:L1, 1:L2,5)';

%Mscaled5 = M5;


pos_shift = 15;

figure(19);
s_contour = contour3(180*omega_space(1:L1)/pi,ds(1:L2),T5,(1:10)*365);
hold off;
figure(11)
K_contours = 1;
while K_contours<length(s_contour)
    H_contours = s_contour(1,K_contours);
    N_contours = s_contour(2,K_contours);
    X_contours = s_contour(1,K_contours+1:K_contours+N_contours);
    Y_contours = s_contour(2,K_contours+1:K_contours+N_contours);

    Z_contours = interp2(X_axis,Y_axis,Mscaled5, X_contours, Y_contours);
    K_contours = K_contours+N_contours+1;
    %plot3(X_contours, Y_contours,Z_contours,'black','HandleVisibility','off')
    pos_i = 15;
    plot3(X_contours(1:pos_i), Y_contours(1:pos_i),Z_contours(1:pos_i),'black','HandleVisibility','off')
    plot3(X_contours(pos_i+pos_shift:end), Y_contours(pos_i+pos_shift:end),Z_contours(pos_i+pos_shift:end),'black','HandleVisibility','off')
    %  plot3(PVevery_a(1:pos_i,j,1), PVevery_a(1:pos_i,j,2), ANevery_a(1:pos_i,j)/(2*pi), 'black', 'HandleVisibility','off')
    %  plot3(PVevery_a(pos_i+3:end,j,1), PVevery_a(pos_i+3:end,j,2), ANevery_a(pos_i+3:end,j)/(2*pi), 'black', 'HandleVisibility','off')
    text_label = num2str(H_contours/365);
    h = text((X_contours(pos_i+pos_shift)+X_contours(pos_i))/2, (Y_contours(pos_i+pos_shift)+Y_contours(pos_i))/2, (Z_contours(pos_i+pos_shift)+Z_contours(pos_i))/2, text_label,'HorizontalAlignment','center',...
     'VerticalAlignment','Bottom','FontName','consolas','FontSize',9); 
    %disp(K_contours)
end

Mscaled1 = 100*M1/m0;
s = surf(180*omega_space(1:L1)/pi,ds(1:L2),Mscaled1,'DisplayName','Family 1', 'FaceColor','cyan','FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';

figure(19);
s_contour = contour3(180*omega_space(1:L1)/pi,ds(1:L2),T1,(1:10)*365);
figure(11)
hold on;
K_contours = 1;
while K_contours<length(s_contour)
    H_contours = s_contour(1,K_contours);
    N_contours = s_contour(2,K_contours);
    X_contours = s_contour(1,K_contours+1:K_contours+N_contours);
    Y_contours = s_contour(2,K_contours+1:K_contours+N_contours);
    Z_contours = interp2(X_axis,Y_axis,Mscaled1, X_contours, Y_contours);
    K_contours = K_contours+N_contours+1;
    %plot3(X_contours, Y_contours,Z_contours,'black', 'HandleVisibility','off')
    pos_i = 15;
    plot3(X_contours(1:pos_i), Y_contours(1:pos_i),Z_contours(1:pos_i),'black','HandleVisibility','off')
    plot3(X_contours(pos_i+pos_shift:end), Y_contours(pos_i+pos_shift:end),Z_contours(pos_i+pos_shift:end),'black','HandleVisibility','off')
    %  plot3(PVevery_a(1:pos_i,j,1), PVevery_a(1:pos_i,j,2), ANevery_a(1:pos_i,j)/(2*pi), 'black', 'HandleVisibility','off')
    %  plot3(PVevery_a(pos_i+3:end,j,1), PVevery_a(pos_i+3:end,j,2), ANevery_a(pos_i+3:end,j)/(2*pi), 'black', 'HandleVisibility','off')
    text_label = num2str(H_contours/365);
    h = text((X_contours(pos_i)+X_contours(pos_i))/2, (Y_contours(pos_i)+Y_contours(pos_i))/2, (Z_contours(pos_i)+Z_contours(pos_i))/2, text_label,'HorizontalAlignment','center',...
     'VerticalAlignment','Bottom','FontName','consolas','FontSize',9); 
    %disp(K_contours)
end
%C2 = rescale(M1,256,511);
Mscaled2 = 100*M2/m0;

%Mscaled2 = M2;
s = surf(180*omega_space(1:L1)/pi,ds(1:L2),Mscaled2,'DisplayName','Family 2', 'FaceColor','red','FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';
figure(19);
hold on;
s_contour = contour3(180*omega_space(1:L1)/pi,ds(1:L2),T2,(1:10)*365);
figure(11)
K_contours = 1;
while K_contours<length(s_contour)
    H_contours = s_contour(1,K_contours);
    N_contours = s_contour(2,K_contours);
    X_contours = s_contour(1,K_contours+1:K_contours+N_contours);
    Y_contours = s_contour(2,K_contours+1:K_contours+N_contours);
    Z_contours = interp2(X_axis,Y_axis,Mscaled2, X_contours, Y_contours);
    K_contours = K_contours+N_contours+1;
    %plot3(X_contours, Y_contours,Z_contours,'black','HandleVisibility','off')
    pos_i = 15;
    if length(X_contours) >= pos_i+pos_shift
        plot3(X_contours(1:pos_i), Y_contours(1:pos_i),Z_contours(1:pos_i),'black','HandleVisibility','off')
        plot3(X_contours(pos_i+pos_shift:end), Y_contours(pos_i+pos_shift:end),Z_contours(pos_i+pos_shift:end),'black','HandleVisibility','off')
        text_label = num2str(H_contours/365);
       h = text((X_contours(pos_i)+X_contours(pos_i))/2, (Y_contours(pos_i)+Y_contours(pos_i))/2, (Z_contours(pos_i)+Z_contours(pos_i))/2, text_label,'HorizontalAlignment','center',...
     'VerticalAlignment','Bottom','FontName','consolas','FontSize',9); 
    else
        plot3(X_contours, Y_contours,Z_contours,'black','HandleVisibility','off')
    end%  plot3(PVevery_a(1:pos_i,j,1), PVevery_a(1:pos_i,j,2), ANevery_a(1:pos_i,j)/(2*pi), 'black', 'HandleVisibility','off')
    %  plot3(PVevery_a(pos_i+3:end,j,1), PVevery_a(pos_i+3:end,j,2), ANevery_a(pos_i+3:end,j)/(2*pi), 'black', 'HandleVisibility','off')
    
    
    %disp(K_contours)
end

Mscaled3 = 100*M3/m0;

%Mscaled3 = M3;
s = surf(180*omega_space(1:L1)/pi,ds(1:L2),Mscaled3,'DisplayName','Family 3', 'FaceColor','magenta','FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';

figure(19);
s_contour = contour3(180*omega_space(1:L1)/pi,ds(1:L2),T3,(1:10)*365);
figure(11)
K_contours = 1;
while K_contours<length(s_contour)
    H_contours = s_contour(1,K_contours);
    N_contours = s_contour(2,K_contours);
    X_contours = s_contour(1,K_contours+1:K_contours+N_contours);
    Y_contours = s_contour(2,K_contours+1:K_contours+N_contours);
    Z_contours = interp2(X_axis,Y_axis,Mscaled3, X_contours, Y_contours);
    K_contours = K_contours+N_contours+1;
    %plot3(X_contours, Y_contours,Z_contours,'black','HandleVisibility','off')
    pos_i = 15;
    if length(X_contours) >= pos_i+pos_shift
        plot3(X_contours(1:pos_i), Y_contours(1:pos_i),Z_contours(1:pos_i),'black','HandleVisibility','off')
        plot3(X_contours(pos_i+pos_shift:end), Y_contours(pos_i+pos_shift:end),Z_contours(pos_i+pos_shift:end),'black','HandleVisibility','off')
        text_label = num2str(H_contours/365);
        h = text((X_contours(pos_i)+X_contours(pos_i))/2, (Y_contours(pos_i)+Y_contours(pos_i))/2, (Z_contours(pos_i)+Z_contours(pos_i))/2, text_label,'HorizontalAlignment','center',...
     'VerticalAlignment','Bottom','FontName','consolas','FontSize',9); 
    else
        plot3(X_contours, Y_contours,Z_contours,'black','HandleVisibility','off')
    end
    %disp(K_contours)
end

Mscaled4 = 100*M4/m0;

%Mscaled4 = M4;
s = surf(180*omega_space(1:L1)/pi,ds(1:L2),Mscaled4,'DisplayName','Family 4', 'FaceColor','green','FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';

figure(19);
s_contour = contour3(180*omega_space(1:L1)/pi,ds(1:L2),T4,(1:10)*365);
hold off;
figure(11)
K_contours = 1;
while K_contours<length(s_contour)
    H_contours = s_contour(1,K_contours);
    N_contours = s_contour(2,K_contours);
    X_contours = s_contour(1,K_contours+1:K_contours+N_contours);
    Y_contours = s_contour(2,K_contours+1:K_contours+N_contours);
    Z_contours = interp2(X_axis,Y_axis,Mscaled4, X_contours, Y_contours);
    K_contours = K_contours+N_contours+1;
    %plot3(X_contours, Y_contours,Z_contours,'black','HandleVisibility','off')
    pos_i = 15;
    if length(X_contours) >= pos_i+pos_shift
        plot3(X_contours(1:pos_i), Y_contours(1:pos_i),Z_contours(1:pos_i),'black','HandleVisibility','off')
        plot3(X_contours(pos_i+pos_shift:end), Y_contours(pos_i+pos_shift:end),Z_contours(pos_i+pos_shift:end),'black','HandleVisibility','off')
        text_label = num2str(H_contours/365);
        h = text((X_contours(pos_i)+X_contours(pos_i))/2, (Y_contours(pos_i)+Y_contours(pos_i))/2, (Z_contours(pos_i)+Z_contours(pos_i))/2, text_label,'HorizontalAlignment','center',...
     'VerticalAlignment','Bottom','FontName','consolas','FontSize',9);  
    else
        plot3(X_contours, Y_contours,Z_contours,'black','HandleVisibility','off')
    end
    %disp(K_contours)
end


%plot3(X_contours(1:2), Y_contours(1:2),Z_contours(1:2),'black','DisplayName','Изолинии времени перелёта, годы')

%colormap(s,spring)
%xlabel('Разность фаз, градусов')
xticks([X_axis(1) (X_axis(1)+X_axis(end))/2 X_axis(end)])
xticklabels({'-\pi','0','\pi'})
xlabel('Initial synodic phase')
ylabel('Angular distance, revolutions')
zlabel('Propellant consumption, %')
%zlabel('J')

%set(gca,'zscale','log')
hold off;
xlim([-180, 180])
legend;
% dcm = datacursormode;
% dcm.Enable = 'on';
% dcm.UpdateFcn = @(obj,event_obj)datatipWithSubscript(obj,event_obj);

figure(19);
close;
%% числа обусловленности
CONDevery_fix = CONDevery;
CONDevery_fix(CONDevery_fix==1)=nan;
%CONDevery_fix(CONDevery_fix>1e5)=1e5;
ax = figure(12);

%cmap1 = colormap(ax,winter);
%cmap2 = colormap(ax,summer);
%mixed_cmap = cat(1,cmap1, cmap2);
%colormap(ax,mixed_cmap);
COND1 = CONDevery_fix(1:L1, 1:L2,1)';
%C1 = rescale(COND1,0,255);
s = surf(omega_space(1:L1),ds(1:L2),COND1,'DisplayName','Семейство 1', 'FaceColor','cyan');
s.EdgeColor = 'none';

hold on;
COND2=CONDevery_fix(1:L1, 1:L2,2)';
%C2 = rescale(COND2,256,511);
s = surf(omega_space(1:L1),ds(1:L2),COND2,'DisplayName','Семейство 2', 'FaceColor','red');
s.EdgeColor = 'none';

COND3=CONDevery_fix(1:L1, 1:L2,3)';
%C2 = rescale(COND2,256,511);
s = surf(omega_space(1:L1),ds(1:L2),COND3,'DisplayName','Семейство 3', 'FaceColor','magenta');
s.EdgeColor = 'none';

COND4=CONDevery_fix(1:L1, 1:L2,4)';
%C2 = rescale(COND2,256,511);
s = surf(omega_space(1:L1),ds(1:L2),COND4,'DisplayName','Семейство 4', 'FaceColor','green');
s.EdgeColor = 'none';
set(gca,'zscale','log')
xlabel('Долгота')
ylabel('Угловая дальность, витков')
zlabel('Число обусловленности')
hold off;
%% цилиндр
Mevery_fix = Jevery;
Mevery_fix(Mevery_fix==0)=nan;
ax = figure(20);
M1 = Mevery_fix(1:L1, 1:L2,1)';
%cmap1 = colormap(ax, summer);
%cmap2 = colormap(ax, cool);
%C1 = rescale(M1,0,255);
%mixed_cmap = cat(1,cmap1, cmap2);
%colormap(ax,mixed_cmap);
%Mscaled1 = 100*M1/m0;
M_transform = @(x)x;
%M_transform = @(x)x;
Mscaled1 = M_transform(M1);
z = repmat(ds(1:L2), L1,1)';
x = Mscaled1.*cos(omega_space(1:L1));
y = Mscaled1.*sin(omega_space(1:L1));
s = surf(x,y,z,'DisplayName','Семейство 1', 'FaceColor','cyan');
s.EdgeColor = 'none';
X_axis=omega_space(1:L1);
figure(19);
s_contour = contour3(X_axis,ds(1:L2),T1,(1:10)*365);
hold on;
figure(20)
hold on;
K_contours = 1;
while K_contours<length(s_contour)
    H_contours = s_contour(1,K_contours);
    N_contours = s_contour(2,K_contours);
    X_contours = s_contour(1,K_contours+1:K_contours+N_contours);
    Y_contours = s_contour(2,K_contours+1:K_contours+N_contours);
    Z_contours = interp2(X_axis,Y_axis,Mscaled1, X_contours, Y_contours);
    K_contours = K_contours+N_contours+1;
    plot3(Z_contours.*cos(X_contours), Z_contours.*sin(X_contours), Y_contours,'black', 'HandleVisibility','off')
    %disp(K_contours)
end
%C2 = rescale(M1,256,511);
M2 = Mevery_fix(1:L1, 1:L2,2)';
%Mscaled2 = 100*M2/m0;
Mscaled2 = M_transform(M2);
z = repmat(ds(1:L2), L1,1)';
x = Mscaled2.*cos(omega_space(1:L1));
y = Mscaled2.*sin(omega_space(1:L1));
s = surf(x,y,z,'DisplayName','Семейство 2', 'FaceColor','red');
s.EdgeColor = 'none';
figure(19);
s_contour = contour3(X_axis,ds(1:L2),T2,(1:10)*365);
figure(20)
K_contours = 1;
while K_contours<length(s_contour)
    H_contours = s_contour(1,K_contours);
    N_contours = s_contour(2,K_contours);
    X_contours = s_contour(1,K_contours+1:K_contours+N_contours);
    Y_contours = s_contour(2,K_contours+1:K_contours+N_contours);
    Z_contours = interp2(X_axis,Y_axis,Mscaled2, X_contours, Y_contours);
    K_contours = K_contours+N_contours+1;
    plot3(Z_contours.*cos(X_contours), Z_contours.*sin(X_contours), Y_contours,'black', 'HandleVisibility','off')
    %disp(K_contours)
end
% 
M3 = Mevery_fix(1:L1, 1:L2,3)';
Mscaled3 = M_transform(M3);
%Mscaled3 = 100*M3/m0;
z = repmat(ds(1:L2), L1,1)';
x = Mscaled3.*cos(omega_space(1:L1));
y = Mscaled3.*sin(omega_space(1:L1));
s = surf(x,y,z,'DisplayName','Семейство 3', 'FaceColor','magenta');
s.EdgeColor = 'none';
figure(19);
s_contour = contour3(X_axis,ds(1:L2),T3,(1:10)*365);
figure(20)
K_contours = 1;
while K_contours<length(s_contour)
    H_contours = s_contour(1,K_contours);
    N_contours = s_contour(2,K_contours);
    X_contours = s_contour(1,K_contours+1:K_contours+N_contours);
    Y_contours = s_contour(2,K_contours+1:K_contours+N_contours);
    Z_contours = interp2(X_axis,Y_axis,Mscaled3, X_contours, Y_contours);
    K_contours = K_contours+N_contours+1;
    plot3(Z_contours.*cos(X_contours), Z_contours.*sin(X_contours), Y_contours,'black', 'HandleVisibility','off')
    %disp(K_contours)
end
% 
M4 = Mevery_fix(1:L1, 1:L2,4)';
%Mscaled4 = 100*M4/m0;
Mscaled4 = M_transform(M4);
z = repmat(ds(1:L2), L1,1)';
x = Mscaled4.*cos(omega_space(1:L1));
y = Mscaled4.*sin(omega_space(1:L1));
s = surf(x,y,z,'DisplayName','Семейство 4', 'FaceColor','green');
s.EdgeColor = 'none';
figure(19);
s_contour = contour3(X_axis,ds(1:L2),T4,(1:10)*365);
figure(20)
K_contours = 1;
while K_contours<length(s_contour)
    H_contours = s_contour(1,K_contours);
    N_contours = s_contour(2,K_contours);
    X_contours = s_contour(1,K_contours+1:K_contours+N_contours);
    Y_contours = s_contour(2,K_contours+1:K_contours+N_contours);
    Z_contours = interp2(X_axis,Y_axis,Mscaled4, X_contours, Y_contours);
    K_contours = K_contours+N_contours+1;
    plot3(Z_contours.*cos(X_contours), Z_contours.*sin(X_contours), Y_contours,'black', 'HandleVisibility','off')
    %disp(K_contours)
end
%colormap(s,spring)
%xlabel('Разность фаз, градусов')
%ylabel('Угловая дальность, витков')
zlabel('Угловая дальность, витков')
hold off;
%xlim([-180, 180])
legend;
%% вычисление чисел обусловленности на прямоугольнике
% t_point = 222;
% s_point = 100;
% vecPX1 = reshape(PXevery(t_point,s_point,:,1), [1,8]);
% vecPX2 = reshape(PXevery(t_point,s_point,:,2), [1,8]);
% N_points_cond = 10;
% CONDparallelogram = zeros(N_points_cond,N_points_cond);
% CONDparallelogram_points = linspace(0,1,N_points_cond);
% for i = 1:N_points_cond
%     for j = 1:N_points_cond
%        px_mixed = vecPX1*CONDparallelogram_points(i)+vecPX2*CONDparallelogram_points(j);
%        C = calculateCondOnly(dt(t_point), px_mixed, ds(s_point), PHIevery(t_point,s_point));
%         CONDparallelogram(i,j)=C;
%     end
% end
% figure(13);
% surf(CONDparallelogram_points,CONDparallelogram_points,CONDparallelogram);
% xlabel('1 family')
% ylabel('2 family')

%%
phi0 = 0;
mug = 1;

st.t = t_start;
st.planet = 'Earth';
st.mode = orbits;
st.delta_omega = omega_space(1,1);

[r0, v0]=planetModel(st);
eul = [pi/12 pi/4 pi/12];
rotmZYX = eul2rotm(eul);

r0 = [rotmZYX*r0'/ae; 0]*1e+03;
V0 = [rotmZYX*V0'/V_unit; 0]*1e+03;
u0 = rToU(r0, phi0);
w0 = vFromV(V0,r0,mug,phi0);
ortdgduv = get_ortdgduv(u0,w0);
dFdX = get_dgduv(u0,w0);
%%
cosineDistEvery=zeros(L1,L2);
cosineSimEvery=zeros(L1,L2);
PVevery=zeros(L1,L2, 3, 5);
PRevery=zeros(L1,L2, 3, 5);
for i = 1:L1
    for j = 1:L2
        for F = 1:5
            vecPX1 = reshape(PXevery(i,j,:,F), [1,8]);
            vecPV_cartesian = a_reactive(u0,w0,vecPX1(1:4)',vecPX1(5:8)');
            b = vecPX1'-dFdX' * [0; 0; 0; vecPV_cartesian(1:3)]; %TODO переписать в символьных
            vecPR_cartesian = dFdX(1:3,:)'\b;

            vecPV_cartesian = rotmZYX^(-1)*vecPV_cartesian(1:3);%/(ae/sqrt(mug_0)).^2;
            vecPR_cartesian = rotmZYX^(-1)*vecPR_cartesian;%/(ae/sqrt(mug_0)).^2;
            PVevery(i,j,:,F) = vecPV_cartesian';
            PRevery(i,j,:,F) = vecPR_cartesian';
        end
    end
end
%%
for i = 1:L1
    for j = 1:L2
        vecPX1 = reshape(PVevery(i,j,:,1), [1,3]);
        vecPX2 = reshape(PVevery(i,j,:,2), [1,3]);
        cosSim = vecPX1*vecPX2'/norm(vecPX1)/norm(vecPX2);
        %cosDist= 1-abs(cosSim);
        cosDist=norm(vecPX1)/norm(vecPX2);
        cosineDistEvery(i,j)=cosDist;
        cosineSimEvery(i,j)=cosSim;
    end
end
cosineDistEvery(cosineDistEvery==0)=NaN;
cosineSimEvery(cosineSimEvery==0)=NaN;
figure(13);
s = surf(omega_space(1:L1),ds(1:L2),cosineSimEvery(:, 1:L2)');
s.EdgeColor = 'none';
xlabel('Разность фаз')
ylabel('Угловая дальность, витков')
zlabel('Косинусное сходство')

PV_1 = squeeze(PVevery(1,:,1:2,1));
PV_1_norm = abs(vecnorm(PV_1,2,2));
PV_1_angle = atan2(PV_1(:,2), PV_1(:,1)); %rad

PV_2 = squeeze(PVevery(1,:,1:2,2));
PV_2_norm = abs(vecnorm(PV_2,2,2));
PV_2_angle = atan2(PV_2(:,2), PV_2(:,1)); %rad

PV_3 = squeeze(PVevery(1,:,1:2,3));
PV_3_norm = abs(vecnorm(PV_3,2,2));
PV_3_angle = atan2(PV_3(:,2), PV_3(:,1)); %rad

PV_4 = squeeze(PVevery(1,:,1:2,4));
PV_4_norm = abs(vecnorm(PV_4,2,2));
PV_4_angle = atan2(PV_4(:,2), PV_4(:,1)); %rad

figure(14)
plot(THETA(1:L2),PV_1_norm,'cyan', 'DisplayName','Семейство 1')
hold on
plot(THETA(50:L2),PV_2_norm(50:L2),'red', 'DisplayName','Семейство 2')
plot(THETA(100:L2),PV_3_norm(100:L2),'magenta','DisplayName','Семейство 3')
plot(THETA(140:L2),PV_4_norm(140:L2),'green','DisplayName','Семейство 4')
hold off
title('Pv norm')
xlabel('Угловая дальность, витков')
ylabel('Норма базис-вектора')
legend;
xlim([0 6.5])
grid on;

figure(15)
plot(THETA(1:L2),PV_1_angle,'cyan', 'DisplayName','Семейство 1')
hold on
plot(THETA(50:L2),PV_2_angle(50:L2),'red', 'DisplayName','Семейство 2')
plot(THETA(100:L2),PV_3_angle(100:L2),'magenta','DisplayName','Семейство 3');
plot(THETA(140:L2),PV_4_angle(140:L2),'green','DisplayName','Семейство 4');
hold off
title('Pv angle')
xlabel('Угловая дальность, витков')
ylabel('Угол базис-вектора\n относительно Ox, радиан')
legend;
xlim([0 6.5])
grid on;

PR_1 = squeeze(PRevery(1,:,1:2,1));
PR_1_norm = abs(vecnorm(PR_1,2,2));
PR_1_angle = atan2(PR_1(:,2), PR_1(:,1)); %rad
PR_1_angle(PR_1_angle<-1)=PR_1_angle(PR_1_angle<-1)+2*pi;

PR_2 = squeeze(PRevery(1,:,1:2,2));
PR_2_norm = abs(vecnorm(PR_2,2,2));
PR_2_angle = atan2(PR_2(:,2), PR_2(:,1)); %rad
PR_2_angle(PR_2_angle<-1)=PR_2_angle(PR_2_angle<-1)+2*pi;

PR_3 = squeeze(PRevery(1,:,1:2,3));
PR_3_norm = abs(vecnorm(PR_3,2,2));
PR_3_angle = atan2(PR_3(:,2), PR_3(:,1)); %rad
PR_3_angle(PR_3_angle>1)=PR_3_angle(PR_3_angle>1)-2*pi;

PR_4 = squeeze(PRevery(1,:,1:2,4));
PR_4_norm = abs(vecnorm(PR_4,2,2));
PR_4_angle = atan2(PR_4(:,2), PR_4(:,1)); %rad

figure(16)
plot(THETA(1:L2),PR_1_norm,'cyan', 'DisplayName','Семейство 1')
hold on
plot(THETA(50:L2),PR_2_norm(50:L2),'red', 'DisplayName','Семейство 2')
plot(THETA(100:L2),PR_3_norm(100:L2),'magenta','DisplayName','Семейство 3')
plot(THETA(140:L2),PR_4_norm(140:L2),'green','DisplayName','Семейство 4')
hold off
xlim([0 6.5])
xlabel('Угловая дальность, витков')
ylabel('Норма производной базис-вектора')
title('Pr norm')
legend;
grid on;

figure(17)
plot(THETA(1:L2),PR_1_angle,'cyan', 'DisplayName','Семейство 1')
hold on
plot(THETA(50:L2),PR_2_angle(50:L2),'red', 'DisplayName','Семейство 2')
plot(THETA(100:L2),PR_3_angle(100:L2),'magenta','DisplayName','Семейство 3');
plot(THETA(140:L2),PR_4_angle(140:L2),'green','DisplayName','Семейство 4');
hold off
xlabel('Угловая дальность, витков')
ylabel('Угол производной базис-вектора\n относительно Ox, радиан')
title('Pr angle')
xlim([0 6.5])
legend;
grid on;

%% 3д график для угла базис-вектора
figure(21);
PVevery_fix = PVevery;
PVevery_fix(PVevery_fix==0)=nan;

PV_1 = squeeze(PVevery_fix(:,:,1:2,1));
PV_2 = squeeze(PVevery_fix(:,:,1:2,2));
PV_3 = squeeze(PVevery_fix(:,:,1:2,3));
PV_4 = squeeze(PVevery_fix(:,:,1:2,4));
PV_5 = squeeze(PVevery_fix(:,:,1:2,5));

z = repmat(ds(1:L2), L1,1);
s = surf(PV_1(:,:,1),PV_1(:,:,2),z,'DisplayName','Семейство 1', 'FaceColor','cyan','FaceAlpha',0.2);
s.EdgeColor = 'none';
hold on;
s = surf(PV_2(:,:,1),PV_2(:,:,2),z,'DisplayName','Семейство 2', 'FaceColor','red','FaceAlpha',0.2);
s.EdgeColor = 'none';
s = surf(PV_3(:,:,1),PV_3(:,:,2),z,'DisplayName','Семейство 3', 'FaceColor','magenta','FaceAlpha',0.2);
s.EdgeColor = 'none';
s = surf(PV_4(:,:,1),PV_4(:,:,2),z,'DisplayName','Семейство 4', 'FaceColor','green','FaceAlpha',0.2);
s.EdgeColor = 'none';
s = surf(PV_5(:,:,1),PV_5(:,:,2),z,'DisplayName','Семейство 0', 'FaceColor','green','FaceAlpha',0.2);
s.EdgeColor = 'none';
%фиксированная разность фаз
% for i = 1:36:L1
%     plot3(PV_1(i,:,1), PV_1(i,:,2), ds(1:L2), 'black', 'HandleVisibility','off')
% end
% for i = 1:36:L1
%     plot3(PV_2(i,:,1), PV_2(i,:,2), ds(1:L2), 'black', 'HandleVisibility','off')
% end
% for i = 1:36:L1
%     plot3(PV_3(i,:,1), PV_3(i,:,2), ds(1:L2), 'black', 'HandleVisibility','off')
% end
% for i = 1:36:L1
%     plot3(PV_4(i,:,1), PV_4(i,:,2), ds(1:L2), 'black', 'HandleVisibility','off')
% end
% 


hold off;
xlabel('Базис-вектор, X')
ylabel('Базис-вектор, Y')
zlabel('Угловая дальность, витков')
%%
%График для среза семейств
Mevery_fix = Mevery;
Mevery_fix(Mevery_fix==0)=nan;
ax = figure(18);
M1 = Mevery_fix(1,1:L2,1)';
%cmap1 = colormap(ax, summer);
%cmap2 = colormap(ax, cool);
%C1 = rescale(M1,0,255);
%mixed_cmap = cat(1,cmap1, cmap2);
%colormap(ax,mixed_cmap);

Mscaled1 = 100*M1/m0;
plot(THETA(1:L2),Mscaled1,'cyan', 'DisplayName','Семейство 1');
hold on;
%C2 = rescale(M1,256,511);
M2 = Mevery_fix(1, 1:L2,2)';
Mscaled2 = 100*M2/m0;
plot(THETA(1:L2),Mscaled2,'red', 'DisplayName','Семейство 2');

M3 = Mevery_fix(1, 1:L2,3)';
Mscaled3 = 100*M3/m0;
plot(THETA(1:L2),Mscaled3,'magenta','DisplayName','Семейство 3');

M4 = Mevery_fix(1, 1:L2,4)';
Mscaled4 = 100*M4/m0;
plot(THETA(1:L2),Mscaled4,'green','DisplayName','Семейство 4');
%colormap(s,spring)
xlabel('Угловая дальность, витков')
ylabel('Затраты топлива, %')
hold off;
grid on;
legend;
%%
%График для среза семейств по времени
Mevery_fix = Jevery;
Mevery_fix(Mevery_fix==0)=nan;
ax = figure(18);
M1 = Mevery_fix(1,1:L2,1)';
%cmap1 = colormap(ax, summer);
%cmap2 = colormap(ax, cool);
%C1 = rescale(M1,0,255);
%mixed_cmap = cat(1,cmap1, cmap2);
%colormap(ax,mixed_cmap);
FaceAlpha=0.6;
Mscaled1 = 100*M1/m0;
cyan_color = [0 1 1 FaceAlpha];
plot(T1(1:L2,1)/365.25,Mscaled1,'Color',cyan_color, 'DisplayName','Family 2');
hold on;
%C2 = rescale(M1,256,511);
M2 = Mevery_fix(1, 1:L2,2)';
Mscaled2 = 100*M2/m0;
red_color = [1 0 0 FaceAlpha];
plot(T2(1:L2,1)/365.25,Mscaled2,'Color', red_color, 'DisplayName','Family 3');

M3 = Mevery_fix(1, 1:L2,3)';
Mscaled3 = 100*M3/m0;
magenta_color = [1 0 1 FaceAlpha];
plot(T3(1:L2,1)/365.25,Mscaled3,'Color',magenta_color,'DisplayName','Family 4');

M4 = Mevery_fix(1, 1:L2,4)';
Mscaled4 = 100*M4/m0;
green_color = [0 1 0 FaceAlpha];
plot(T4(1:L2,1)/365.25,Mscaled4,'Color',green_color,'DisplayName','Family 5');

% не попадает, поэтому закомментировал
% M5 = Mevery_fix(361, 1:L2,5)';
% Mscaled5 = 100*M5/m0;
% plot(T5(1:L2,361)/365.25,Mscaled5,'blue','DisplayName','Семейство 1');
%plot(T1(1:70,1)/365.25,Mscaled1(1:70),'Color','blue', 'DisplayName','Pareto front with a fixed start date');
%plot(T2(102:end,1)/365.25,Mscaled2(102:end),'Color','blue', 'HandleVisibility','off');
%plot(Tevery_a(:,113)/365.25,100*Mevery_a(:,113)/m0,'black','DisplayName','Global Pareto front');
%colormap(s,spring)
xlabel('Time of flight, years')
ylabel('Propellant consumption, %')
hold off;
grid on;
legend;
%% время
Tevery_fix = Tevery;
Tevery_fix(Tevery_fix==0)=nan;
ax = figure(19);
T1 = Tevery_fix(1:L1, 1:L2,1)';

s = surf(180*omega_space(1:L1)/pi,ds(1:L2),T1,'DisplayName','Семейство 1', 'FaceColor','cyan');
s.EdgeColor = 'none';
hold on;
s_contour = contour3(180*omega_space(1:L1)/pi,ds(1:L2),T1,(1:10)*365);
K_contours = 1;
while K_contours<length(s_contour)
    H_contours = s_contour(1,K_contours);
    N_contours = s_contour(2,K_contours);
    X_contours = s_contour(1,K_contours+1:K_contours+N_contours);
    Y_contours = s_contour(2,K_contours+1:K_contours+N_contours);
    Z_contours = interp2(X_axis,Y_axis,Mscaled1, X_contours, Y_contours);
    K_contours = K_contours+N_contours+1;
    disp(K_contours)
end
%C2 = rescale(M1,256,511);
% M2 = Mevery_fix(1:L1, 1:L2,2)';
% Mscaled2 = 100*M2/m0;
% s = surf(180*omega_space(1:L1)/pi,ds(1:L2),Mscaled2,'DisplayName','Семейство 2', 'FaceColor','red');
% s.EdgeColor = 'none';
% 
% M3 = Mevery_fix(1:L1, 1:L2,3)';
% Mscaled3 = 100*M3/m0;
% s = surf(180*omega_space(1:L1)/pi,ds(1:L2),Mscaled3,'DisplayName','Семейство 3', 'FaceColor','magenta');
% s.EdgeColor = 'none';
% 
% M4 = Mevery_fix(1:L1, 1:L2,4)';
% Mscaled4 = 100*M4/m0;
% s = surf(180*omega_space(1:L1)/pi,ds(1:L2),Mscaled4,'DisplayName','Семейство 4', 'FaceColor','green');
% s.EdgeColor = 'none';
%colormap(s,spring)
xlabel('Разность фаз, градусов')
ylabel('Угловая дальность, витков')
zlabel('Время перелёта, дни')
hold off;
xlim([-180, 180])
legend;
%%%
function txt = displayCoordinates(~,info)
    x = info.Position(1);
    y = info.Position(2);
    txt = ['(' num2str(x) ', ' num2str(y) ')'];
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
output_txt = {['X: ',num2str(pos(1),4)],...
              ['Y: ',num2str(pos(2),4)]};
% If there is a Z-coordinate in the position, display it as well. Also get
% the size of the object.
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
    dataSize  = size(get(get(event_obj,'Target'),'ZData'));
else
    dataSize  = size(get(get(event_obj,'Target'),'YData'));
end
% Get the linear index
dataIndex = get(event_obj,'DataIndex');
if sum(dataSize>1)>1
    % display subscripts for data with more than 1 dimension
   
    S = cell(1,numel(dataSize));
    %global i_one j_one;
    [i_one, j_one] = cell2mat(S);
    F_one = 1;
    %plotTrajectoryCallback(i_one, j_one, F_one, PXevery, PHIevery)
    [S{:}] = ind2sub(dataSize,dataIndex);
    
    subsString = sprintf('%i,',cell2mat(S));

    %plotOneTrajectoryFlat;
    output_txt{end+1} = ['@ (', subsString(1:end-1), ')'];
else
    % display linear index for vector data
    output_txt{end+1} = ['@ ',num2str(dataIndex)];
end
end

function [skip, i_nearest, j_nearest] = checkNear(arr, i, j, F)
    % checkNear
    s = size(arr);
    if arr(i,j,F)==0
        %уже есть решение
        skip = 0;
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
        skip = -1;
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