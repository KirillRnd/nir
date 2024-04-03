%% выбираем функционал
Mevery_fix = Jevery;
Mevery_fix(Mevery_fix==0)=nan;

Tevery_fix = Tevery/Tunit_mars;
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

X_axis = 180*omega_space(1:L1)/pi;
Y_axis = ds(1:L2);

M1 = Mevery_fix(1:L1, 1:L2,1)';
M2 = Mevery_fix(1:L1, 1:L2,2)';
M3 = Mevery_fix(1:L1, 1:L2,3)';
M4 = Mevery_fix(1:L1, 1:L2,4)';
M5 = Mevery_fix(1:L1, 1:L2,5)';

M1(1:10,181:end) = nan;
T1(1:10,181:end) = nan;
for i=181:L1
    for j=10:20
        if j <= 10+(i-181)*7/110
            M1(j,i) = nan;
            T1(j,i) = nan;
            AN1(j,i) = nan;
        end
    end
end

M2(10:39,161:end) = nan;
T2(10:39,161:end) = nan;
AN2(10:39,161:end) = nan;

M3(10:71,160:end) = nan;
T3(10:71,160:end) = nan;
AN3(10:71,160:end) = nan;

M3(10:87,332:end) = nan;
T3(10:87,332:end) = nan;
AN3(10:87,332:end) = nan;

T_united = [T1,T2,T3,T4,T5];
Y_axis_united = repmat(Y_axis,size(T_united,2),1)';
M_united = [M1,M2,M3,M4,M5];
AN_united = [AN1,AN2,AN3,AN4,AN5];

global_minimum_line = zeros(4,size(T_united,1));
for i = 1:size(T_united,1)
    m_local = M_united(i,:);
    m_min = min(m_local,[],'all');
    j_min = find(m_local==m_min);
    global_minimum_line(:,i) = [T_united(i, j_min), AN_united(i, j_min), M_united(i, j_min), Y_axis_united(i, j_min)];
end

%plot everything
figure(22);
FaceAlpha=0.4;

% s = surf(T5,Y_axis,M5,'DisplayName','Семейство 1', 'FaceColor','blue','FaceAlpha',FaceAlpha);
% s.EdgeColor = 'none';
s = surf(T5,AN5,M5,'DisplayName','Семейство 1', 'FaceColor','blue','FaceAlpha',FaceAlpha);
s.EdgeColor = 'none';


hold on;

% s = surf(T1,Y_axis,M1,'DisplayName','Семейство 2', 'FaceColor','cyan','FaceAlpha',FaceAlpha);
% s.EdgeColor = 'none';
% 
% s = surf(T2,Y_axis,M2,'DisplayName','Семейство 3', 'FaceColor','red','FaceAlpha',FaceAlpha);
% s.EdgeColor = 'none';
% 
% s = surf(T3,Y_axis,M3,'DisplayName','Семейство 4', 'FaceColor','magenta','FaceAlpha',FaceAlpha);
% s.EdgeColor = 'none';
% 
% s = surf(T4,Y_axis,M4,'DisplayName','Семейство 5', 'FaceColor','green','FaceAlpha',FaceAlpha);
% s.EdgeColor = 'none';

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
xlim([0, 6])
ylim([0, 6])
set(gca,'fontsize',12)

figure(23);
plot(global_minimum_line(2,:), global_minimum_line(2,:)-global_minimum_line(4,:))
ylabel('Разность фиктивного времени, безразм.')
xlabel('Угловая дальность, безразм.')
set(gca,'fontsize',12)