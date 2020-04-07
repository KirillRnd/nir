function d_b = find_zad_opt(th_mars,gamma,v_opt,t_opt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here



mug = 132712.43994*(10^6)*(10^(3*3));
ae = 149597870700;
Vsc=norm(v_opt)-(mug/(1*ae))^(1/2);
%оценка без ограничений
%z0=[0.0026 0.0005 0 -0.5149*1.0e-09 -0.0444*1.0e-09 0];
z0 = [0.1 0.1 0 0 0 0];
th_earth = +6*pi/4;
r_earth=[ae*cos(th_earth) ae*sin(th_earth) 0];
v_earth=[((mug/(1*ae))^(1/2)+norm(Vsc))*cos(th_earth+pi/2) ((mug/(1*ae))^(1/2)+norm(Vsc))*sin(th_earth+pi/2)  0];


u_dv = 3000*9.81;
m0 = 1e+5;
pv0_alpha = gamma*u_dv/m0;

y0 = cat(2, r_earth, v_earth , z0,...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(eye(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(eye(3),[1,9])...
    )';



%Определяем tf
%T=2*pi*sqrt((1*ae)^3/mug);
%tf=6*T/6;
tf=124*3600*24-t_opt;
optionsInn = odeset('AbsTol',1e-12);

[t,y] = ode45(@(t,y) integrateTraectory(t,y,mug,gamma,th_mars,v_opt,t_opt,...
    norm(y0(7:9))),[0 tf],y0,optionsInn);
plot(y(:,1),y(:,2));
hold on;
axis equal
plot(y(end,1),y(end,2),'r--o');

drdzdt=reshape(y(end,49:66),[3,6]);
drdz=reshape(y(end,13:30),[3,6]);

dfdz = cat(1,drdz,drdzdt);
%th_mars=3*pi/6+pi/15;

%кеплеровы параметры Марс
ex=0.093394;
a=1.523662*ae;
i=0;
w=286.46230/360*2*pi;
omega=49.57/360*2*pi;

E=th_mars;
b = a*sqrt(1-ex^2);
r00=[a*(cos(E)-ex) b*sin(E) 0]';
tmpq=my_eul_to_quat(omega, i, w,"ZXZs");

Atr=quat2dcm(tmpq);
rf=Atr*r00;
rn=rf/norm(rf);

sinTh=b*sin(E)/norm(rf);
p=a*(1-ex^2);
cInt=sqrt(mug * p);

k=Atr*[0 0 1]';

Vr=sqrt(mug/p)*ex*sinTh*rn;
Vn=cInt/norm(rf)*cross(k,rn);
vf=Vn+Vr;


%rf = [1.52*ae*cos(th_mars) 1.52*ae*sin(th_mars) 0];
%vf = [(mug/(1.52*ae))^(1/2)*cos(th_mars+pi/2) (mug/(1.52*ae))^(1/2)*sin(th_mars+pi/2) 0];


b=y(end,1:6)-cat(2,rf',vf');
%оптимизируем траекторию


optionsExt = odeset('AbsTol',1e-16);

[tau,z] = ode45(@(t,y) optimiseToMars(t,y,b,tf,gamma,th_mars,v_opt,...
    t_opt, y0),[0 1],z0, optionsExt);

y0(7:12)=z(end,:);
[t,y] = ode45(@(t,y) integrateTraectory(t,y,mug,gamma,th_mars,v_opt,t_opt,...
    norm(y0(7:9))),[0 tf],y0,optionsInn);
plot(y(:,1),y(:,2),'k');
axis equal
hold on
th = 0:pi/50:2*pi;
%орбита Земли
plot(ae*cos(th),ae*sin(th),'b');
%орбита Марса
plot((p./(1+ex*cos(th))).*cos(th),(p./(1+ex*cos(th))).*sin(th),'r');
%конечная точка траектории
plot(rf(1),rf(2),'b--o')
%Земля
plot(r_earth(1), r_earth(2),'b--o')
%Солнце
plot(0, 0,'y--o')
hold off;
d_r = norm([rf(1)-y(end,1),rf(2)-y(end,2)]);
d_v = norm([vf(1)-y(end,4),vf(2)-y(end,5)]);
d_b = (d_r + d_v*1e+06)*1e-10;
end

