%% выбираем функционал
Mevery_fix = Jevery;
Mevery_fix(Mevery_fix==0)=nan;

ae = 149597870.700; %км
mug = 132712.43994*(10^6); %км3с−2
Tunit_earth = (2*pi/sqrt(mug))*((ae).^(3/2))/(3600*24);
Tunit_mars = (2*pi/sqrt(mug))*((1.52*ae).^(3/2))/(3600*24);

Tevery_fix = Tevery/Tunit_earth;
Tevery_fix(Tevery_fix==0)=nan;
T1 = Tevery_fix(1:L1, 1:L2,1)';
T2 = Tevery_fix(1:L1, 1:L2,2)';
T3 = Tevery_fix(1:L1, 1:L2,3)';
T4 = Tevery_fix(1:L1, 1:L2,4)';
T5 = Tevery_fix(1:L1, 1:L2,5)';


ANevery_fix = ANevery/(2*pi);
ANevery_fix(ANevery_fix==0)=nan;
AN1 = ANevery_fix(1:L1, 1:L2,1)';
AN2 = ANevery_fix(1:L1, 1:L2,2)';
AN3 = ANevery_fix(1:L1, 1:L2,3)';
AN4 = ANevery_fix(1:L1, 1:L2,4)';
AN5 = ANevery_fix(1:L1, 1:L2,5)';

CONDevery_fix = CONDevery;
CONDevery_fix(CONDevery_fix==1)=nan;
COND1 = CONDevery_fix(1:L1, 1:L2,1)';
COND2 = CONDevery_fix(1:L1, 1:L2,2)';
COND3 = CONDevery_fix(1:L1, 1:L2,3)';
COND4 = CONDevery_fix(1:L1, 1:L2,4)';
COND5 = CONDevery_fix(1:L1, 1:L2,5)';

X_axis = 180*omega_space(1:L1)/pi;
Y_axis = ds(1:L2);

M1 = Mevery_fix(1:L1, 1:L2,1)';
M2 = Mevery_fix(1:L1, 1:L2,2)';
M3 = Mevery_fix(1:L1, 1:L2,3)';
M4 = Mevery_fix(1:L1, 1:L2,4)';
M5 = Mevery_fix(1:L1, 1:L2,5)';


drop1Fam_j = 1:10;
drop1Fam_i = 181:L1;
M1(drop1Fam_j,drop1Fam_i) = nan;
T1(drop1Fam_j,drop1Fam_i) = nan;
AN1(drop1Fam_j,drop1Fam_i) = nan;
COND1(drop1Fam_j,drop1Fam_i) = nan;
for i=181:L1
    for j=10:20
        if j <= 10+(i-181)*7/110
            M1(j,i) = nan;
            T1(j,i) = nan;
            AN1(j,i) = nan;
        end
    end
end

drop2Fam_j = 10:39;
drop2Fam_i = 161:L1;
M2(drop2Fam_j,drop2Fam_i) = nan;
T2(drop2Fam_j,drop2Fam_i) = nan;
AN2(drop2Fam_j,drop2Fam_i) = nan;
COND2(drop2Fam_j,drop2Fam_i) = nan;

drop3Fam_j = 10:71;
drop3Fam_i = 160:L1;
M3(drop3Fam_j,drop3Fam_i) = nan;
T3(drop3Fam_j,drop3Fam_i) = nan;
AN3(drop3Fam_j,drop3Fam_i) = nan;
COND3(drop3Fam_j,drop3Fam_i) = nan;

% drop3Fam_j = 10:87;
% drop3Fam_i = 332:L1;
% M3(drop3Fam_j,drop3Fam_i) = nan;
% T3(drop3Fam_j,drop3Fam_i) = nan;
% AN3(drop3Fam_j,drop3Fam_i) = nan;

T_united = [T1,T2,T3,T4,T5];
Y_axis_united = repmat(Y_axis,size(T_united,2),1)';
M_united = [M1,M2,M3,M4,M5];
AN_united = [AN1,AN2,AN3,AN4,AN5];
COND_united = [COND1,COND2,COND3,COND4,COND5];

PXevery_a_MARS = zeros(L2,8); %a = 1.52, Mars
PHIevery_a_MARS = zeros(L2,1);
OMevery_a_MARS = zeros(L2,1);

global_minimum_line = zeros(5,size(T_united,1));

PX_united = [PXevery(:,:,:,1);PXevery(:,:,:,2);PXevery(:,:,:,3);PXevery(:,:,:,4);PXevery(:,:,:,5);];
PHI_united = [PHIevery(:,:,1);PHIevery(:,:,2);PHIevery(:,:,3);PHIevery(:,:,4);PHIevery(:,:,5);];
OM_united = repmat(omega_space,1,5)';
for i = 1:size(T_united,1)
    m_local = M_united(i,:);
    m_min = min(m_local,[],'all');
    j_min = find(m_local==m_min);

    global_minimum_line(:,i) = [T_united(i, j_min), AN_united(i, j_min), M_united(i, j_min), Y_axis_united(i, j_min), COND_united(i, j_min)];

    PXevery_a_MARS(i,:) = PX_united(j_min, i, :); %a = 1.52, Mars
    PHIevery_a_MARS(i) = PHI_united(j_min, i);
    OMevery_a_MARS(i) = OM_united(j_min, 1);
end


%plot everything
figure(22);
FaceAlpha=0.4;

% s = surf(T5,Y_axis,M5,'DisplayName','Семейство 1', 'FaceColor','blue','FaceAlpha',FaceAlpha);
% s.EdgeColor = 'none';
s = surf(T5,AN5,M5,'DisplayName','Семейство 1', 'FaceColor','blue','FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';


hold on;

s = surf(T1,AN1,M1,'DisplayName','Семейство 2', 'FaceColor','cyan','FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';

s = surf(T2,AN2,M2,'DisplayName','Семейство 3', 'FaceColor','red','FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';

s = surf(T3,AN3,M3,'DisplayName','Семейство 4', 'FaceColor','magenta','FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';

s = surf(T4,AN4,M4,'DisplayName','Семейство 5', 'FaceColor','green','FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';


glob_min_line = plot3(global_minimum_line(1,:),global_minimum_line(2,:),global_minimum_line(3,:),'black', 'LineWidth', 1.5,'DisplayName','Глобальный минимум');

s_contour = contour3(T_united,AN_united,M_united,round(logspace(-2,3,11),2),'ShowText','on', 'HandleVisibility','off');
hold off;
zlabel('Значение функционала, безразм.')
xlabel('Длительность перелёта, безразм.')
ylabel('Угловая дальность, безразм.')


set(gca,'zscale','log')
legend;
xlim([0, 12])
ylim([0, 6])
set(gca,'fontsize',12)

figure(23);
plot(global_minimum_line(2,:), global_minimum_line(2,:)-global_minimum_line(4,:))
ylabel('Разность фиктивного времени, безразм.')
xlabel('Угловая дальность, безразм.')
set(gca,'fontsize',12)
%% числа обусловленности
figure(32);
s = surf(T5,AN5,COND5,'DisplayName','Семейство 1', 'FaceColor','blue','FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';


hold on;

s = surf(T1,AN1,COND1,'DisplayName','Семейство 2', 'FaceColor','cyan','FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';

s = surf(T2,AN2,COND2,'DisplayName','Семейство 3', 'FaceColor','red','FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';

s = surf(T3,AN3,COND3,'DisplayName','Семейство 4', 'FaceColor','magenta','FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';

s = surf(T4,AN4,COND4,'DisplayName','Семейство 5', 'FaceColor','green','FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';

glob_min_line = plot3(global_minimum_line(1,:),global_minimum_line(2,:),global_minimum_line(5,:),'black', 'LineWidth', 1.5,'DisplayName','Глобальный минимум');

s_contour = contour3(T_united,AN_united,COND_united,round(logspace(1,9,9),2),'ShowText','on', 'HandleVisibility','off');
hold off;
grid on;
xlim([0, 12])
ylim([0, 6])
set(gca,'zscale','log')
legend;
zlabel('Значение числа обусловленности, безразм.')
xlabel('Длительность перелёта, безразм.')
ylabel('Угловая дальность, безразм.')
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
PVevery_mars=zeros(L2, 3);
PRevery_mars=zeros(L2, 3);
for j = 1:L2
    vecPX1 = reshape(PXevery_a_MARS(j,:), [1,8]);
    vecPV_cartesian = a_reactive(u0,w0,vecPX1(1:4)',vecPX1(5:8)');
    b = vecPX1'-dFdX' * [0; 0; 0; vecPV_cartesian(1:3)]; %TODO переписать в символьных
    vecPR_cartesian = dFdX(1:3,:)'\b;

    vecPV_cartesian = rotmZYX^(-1)*vecPV_cartesian(1:3)/(ae/sqrt(mug_0)).^2;
    vecPR_cartesian = rotmZYX^(-1)*vecPR_cartesian/(ae/sqrt(mug_0)).^2;
    PVevery_mars(j,:) = vecPV_cartesian';
    PRevery_mars(j,:) = vecPR_cartesian';
end
%%

PVevery_fix = PVevery;
PVevery_fix(PVevery_fix==0)=nan;

PV_1 = squeeze(PVevery_fix(:,:,1:2,1));
PV_2 = squeeze(PVevery_fix(:,:,1:2,2));
PV_3 = squeeze(PVevery_fix(:,:,1:2,3));
PV_4 = squeeze(PVevery_fix(:,:,1:2,4));
PV_5 = squeeze(PVevery_fix(:,:,1:2,5));

PV_1(drop1Fam_i,drop1Fam_j,:) = nan;
for i=181:L1
    for j=10:20
        if j <= 10+(i-181)*7/110
            PV_1(i,j,:) = nan;
        end
    end
end

PV_2(drop2Fam_i,drop2Fam_j,:) = nan;

%PV_3(160:end,10:71,:) = nan;

PV_3(drop3Fam_i,drop3Fam_j,:) = nan;

z = repmat(ds(1:L2), L1,1);
figure(28);
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
s = plot3(PVevery_mars(:,1), PVevery_mars(:,2), ds(1:L2),'black', 'LineWidth', 1.5,'DisplayName', 'Решения линии минимума');
legend;
hold off
zlabel('Угловая дальность, безразм.')
xlabel('Pv x')
ylabel('Pv y')

figure(29);
plot3(PVevery_mars(:,1), PVevery_mars(:,2), ds(1:L2), 'black', 'HandleVisibility','off')
grid on;
zlabel('Угловая дальность, безразм.')
xlabel('Pv x')
ylabel('Pv y')


%%
PRevery_fix = PRevery;
PRevery_fix(PRevery_fix==0)=nan;

PR_1 = squeeze(PRevery_fix(:,:,1:2,1));
PR_2 = squeeze(PRevery_fix(:,:,1:2,2));
PR_3 = squeeze(PRevery_fix(:,:,1:2,3));
PR_4 = squeeze(PRevery_fix(:,:,1:2,4));
PR_5 = squeeze(PRevery_fix(:,:,1:2,5));

PR_1(181:end,1:10,:) = nan;
for i=181:L1
    for j=10:20
        if j <= 10+(i-181)*7/110
            PR_1(i,j,:) = nan;
        end
    end
end

PR_2(161:end,10:39,:) = nan;

PR_3(160:end,10:71,:) = nan;

PR_3(332:end,10:87,:) = nan;

z = repmat(ds(1:L2), L1,1);
figure(30);
s = surf(PR_1(:,:,1),PR_1(:,:,2),z,'DisplayName','Семейство 1', 'FaceColor','cyan','FaceAlpha',0.2);
s.EdgeColor = 'none';
hold on;
s = surf(PR_2(:,:,1),PR_2(:,:,2),z,'DisplayName','Семейство 2', 'FaceColor','red','FaceAlpha',0.2);
s.EdgeColor = 'none';
s = surf(PR_3(:,:,1),PR_3(:,:,2),z,'DisplayName','Семейство 3', 'FaceColor','magenta','FaceAlpha',0.2);
s.EdgeColor = 'none';
s = surf(PR_4(:,:,1),PR_4(:,:,2),z,'DisplayName','Семейство 4', 'FaceColor','green','FaceAlpha',0.2);
s.EdgeColor = 'none';
s = surf(PR_5(:,:,1),PR_5(:,:,2),z,'DisplayName','Семейство 0', 'FaceColor','green','FaceAlpha',0.2);
s.EdgeColor = 'none';
s = plot3(PRevery_mars(:,1), PRevery_mars(:,2), ds(1:L2),'black', 'LineWidth', 1.5,'DisplayName', 'Решения линии минимума');
legend;
hold off
zlabel('Угловая дальность, безразм.')
xlabel('Pr x')
ylabel('Pr y')

figure(31);
plot3(PRevery_mars(:,1), PRevery_mars(:,2), ds(1:L2), 'black', 'HandleVisibility','off')
grid on;
zlabel('Угловая дальность, безразм.')
xlabel('Pr x')
ylabel('Pr y')
