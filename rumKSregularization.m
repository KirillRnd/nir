mug = 132712.43994*(10^6)*(10^(3*3));
ae = 149597870700;
%Начальные условия
r0=[1*ae 0];
V0=[0 (mug/(1*ae))^(1/2)];
u0=KS_straight(r0);
L0=KS_L(u0,2);
h0=V0^2/2-mug/norm(r0);
v0=L0'*V0/(2*sqrt(-2*h0));

y0 = cat(2,u0,v0,h0);

rf = [0 1.5*ae 0];
vf = [-(mug/(1.5*ae))^(1/2) 0 0];

%оптимизируем траекторию

[tau,z] = ode45(@(tau,y) odeTrajectoryKS_2d(tau,y,mug));
