mug = 132712.43994*(10^6)*(10^(3*3));
ae = 149597870700;
global alpha_vsc;
global v_opt;
Vsc=norm(v_opt)-(mug/(1*ae))^(1/2);
vsc=[Vsc*cos(alpha_vsc) Vsc*sin(alpha_vsc) 0];
 y0 = cat(2,[1*ae 0 0],[0 (mug/(1*ae))^(1/2) 0] + vsc,[0 0 0],[0 0 0],...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(eye(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(eye(3),[1,9])...
    )';
%Определяем tf
T=2*pi*sqrt((1*ae)^3/mug);
%tf=6*T/6;
tf=(4+124)*3600*24-t_opt;
optionsInn = odeset('AbsTol',1e-12);

[t,y] = ode45(@(t,y) integrateTraectory(t,y,mug),[0 tf],y0,optionsInn);
plot(y(:,1),y(:,2));
axis equal
hold on;
drdzdt=reshape(y(end,49:66),[3,6]);
drdz=reshape(y(end,13:30),[3,6]);

dfdz = cat(1,drdz,drdzdt);
global th_mars;  
%th_mars=3*pi/6+pi/15;
rf = [1.52*ae*cos(th_mars) 1.52*ae*sin(th_mars) 0];
vf = [(mug/(1.52*ae))^(1/2)*cos(th_mars+pi/2) (mug/(1.52*ae))^(1/2)*sin(th_mars+pi/2) 0];


b=y(end,1:6)-cat(2,rf,vf);
%оптимизируем траекторию
z0=[0 0 0 0 0 0];

optionsExt = odeset('AbsTol',1e-16);

[tau,z] = ode45(@(t,y) optimiseToMars(t,y,b,tf),[0 1],z0, optionsExt);

y0 = cat(2,[1*ae 0 0],[0 (mug/(1*ae))^(1/2) 0] + vsc,z(end,:),...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(eye(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(eye(3),[1,9])...
    )';
[t,y] = ode45(@(t,y) integrateTraectory(t,y,mug),[0 tf],y0,optionsInn);
plot(y(:,1),y(:,2));
axis equal
hold on
th = 0:pi/50:2*pi;
%орбита Земли
plot(ae*cos(th),ae*sin(th));
%орбита Марса
plot(1.52*ae*cos(th),1.52*ae*sin(th));
%конечная точка траектории
plot(rf(1),rf(2),'b--o')
%Земля
plot(ae, 0,'b--o')
%Солнце
plot(0, 0,'y--o')
%положение земли в момент времени tf
%plot(rf(1)+b(1),rf(2)+b(2),'r--o')
hold off;
norm([rf(1)-y(end,1),rf(2)-y(end,2)])