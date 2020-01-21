mug = 132712.43994*(10^6)*(10^(3*3));
ae = 149597870700;
%Ќачальные услови€
r0=[1*ae 0]';
V0=[0 (mug/(1*ae))^(1/2)]';
u0=KS_straight(r0,2);
L0=KS_L(u0,2);
h0=V0'*V0/2-mug/norm(r0);

v0=L0'*V0/(2*(-2*h0)^(1/2));
%получил из метода продолжени€
pr0=[-0.0017 0.0078];
pV0=1.0e-08 * [0.0316 -0.1839];

y0 = cat(1,u0,v0,h0,0,0,0,0,0);

rf = [0 1.5*ae]';
Vf = [-(mug/(1.5*ae))^(1/2) 0]';
uf=KS_straight(rf,2);
Lf=KS_L(uf,2);
hf=norm(Vf)^2/2-mug/norm(rf);
vf=Lf'*Vf/(2*sqrt(-2*hf));

yf = cat(1,uf,vf,hf,0,0,0,0,0);
%оптимизируем траекторию


T=2*pi*sqrt((1*ae)^3/mug);
tf=5*T/6;
[tau,y] = ode45(@(tau,y) odeTrajectoryKS_2d(tau,y,mug),[0 1],y0);
u=y(:,1:2);
r=zeros(length(u),2);
for c = 1:length(u)
    r(c,:)=KS_back(u(c,1:2)',2);
end 
axis equal
hold on
th = 0:pi/50:2*pi;
plot(ae*cos(th),ae*sin(th));
plot(r(:,1),r(:,2));


plot(ae, 0,'b--o');
plot(0, 0,'y--o');
%[tau,z] = ode45(@(tau,y) odeTrajectoryKS_2d(tau,y,mug));
hold off;