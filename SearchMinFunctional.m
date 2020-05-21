clc;
clear;
warning('off');
angle = 6*pi/6;
case_traj=2;
[functional, dis, s, y] = trajectorySearch(0,angle,case_traj);
n = 0;
for i = 1:10
    [functional_tmp, dis_tmp, s_tmp, y_tmp] = trajectorySearch(i,angle,case_traj);
    if (functional_tmp < functional) && (dis_tmp <= dis*100)
        functional = functional_tmp
        dis = dis_tmp;
        s = s_tmp;
        y = y_tmp;
        n = i;
    end
end

u = y(:, 1:2);
r=zeros(length(u),2);
for i = 1:length(u)
    rr = u(i,:);
    L = [[rr(1) -rr(2)];
    [rr(2) rr(1)]];
    r(i,:)=L*rr';
end

ae = 149597870700;
mug = 132712.43994*(10^6)*(10^(3*3));

%Проверка "на глаз"
plot(0, 0,'y--o')
hold on;
th = 0:pi/50:2*pi;
plot(ae*cos(th),ae*sin(th),'k');
plot(1.52*ae*cos(th),1.52*ae*sin(th),'r');
plot(r(:, 1), r(:, 2),'b')
plot(r(end, 1), r(end, 2),'bO')
plot(1.52*ae*cos(angle), 1.52*ae*sin(angle),'rO')
axis equal
hold off;