
mug = 132712.43994*(10^6)*(10^(3*3));
ae = 149597870700;
 y0 = cat(2,[1*ae 0 0],[0 (mug/(1*ae))^(1/2) 0],[0 0 0],[0 0 0],...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(eye(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(eye(3),[1,9])...
    )';
%Определяем tf
T=2*pi*sqrt((1*ae)^3/mug);
tf=3*T/12;
optionsInn = odeset('AbsTol',1e-12);

[t,y] = ode45(@(t,y) integrateTraectory(t,y,mug),[0 tf],y0,optionsInn);
plot(y(:,1),y(:,2));
axis equal
drdzdt=reshape(y(end,49:66),[3,6]);
drdz=reshape(y(end,13:30),[3,6]);

dfdz = cat(1,drdz,drdzdt);

rf = [0 1.5*ae 0];
vf = [-(mug/(1.5*ae))^(1/2) 0 0];

b=y(end,1:6)-cat(2,rf,vf);
%оптимизируем траекторию
z0=[0 0 0 0 0 0];

optionsExt = odeset('AbsTol',1e-16);

[tau,z] = ode45(@(t,y) optimiseToMars(t,y,b,tf),[0 1],z0, optionsExt);

y0 = cat(2,[1*ae 0 0],[0 (mug/(1*ae))^(1/2) 0],z(end,:),...
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
plot(ae*cos(th),ae*sin(th));
plot(1.5*ae*cos(th),1.5*ae*sin(th));
plot(rf(1),rf(2),'b--o')
plot(ae, 0,'b--o')
plot(0, 0,'y--o')
plot(rf(1)+b(1),rf(2)+b(2),'r--o')
hold off;