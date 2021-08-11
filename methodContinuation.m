mug = 132712.43994*(10^6)*(10^(3*3));
ae = 149597870700;
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;
days2sec=24*3600;
t0=juliandate(2022,1,1);
dt = 570;
N=1350;
m0=367;
eta=0.45;
n=1;
case_traj = 1;
tf=t0+dt;
planet_end='Mars';
[r0, V0] = planetEphemeris(t0,'SolarSystem','Earth','430');
[rf, Vf] = planetEphemeris(tf,'SolarSystem',planet_end,'430');
eul = [0 pi/4 0];
rotmZYX = eul2rotm(eul);
r0=r0'*1e+03;
V0=V0'*1e+03;
rf=rf'*1e+03;
Vf=Vf'*1e+03;
mu0 = defineMu0(n,r0,V0,rf,dt*days2sec,mug);
mu_tau=@(tau)(mu0+(mug-mu0)*tau);
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);   


r = sym('r', [3 1],'real');

r_norm=sqrt(r(1)^2+r(2)^2+r(3)^2);
assume(r_norm,'real');
sym_dUdr=-mug*r/r_norm^3;

sym_ddUdrdr=jacobian(sym_dUdr, r);
pV = sym('pV', [3 1],'real');
assume(sym_ddUdrdr,'real');
ddUdrdrpV=sym_ddUdrdr*pV;
assume(ddUdrdrpV,'real');
jac_sym_ddUdrdr=jacobian(ddUdrdrpV',r);

dUdr = matlabFunction(sym_dUdr, 'Vars', {r});
ddUdrdr = matlabFunction(sym_ddUdrdr, 'Vars', {r});
jac_ddUdrdr = matlabFunction(jac_sym_ddUdrdr, 'Vars', {r,pV});

%y0=[1 0 0 0 -mug^(1/2) 0 0 0 0 0 0 0]';
y0 = cat(2,r0',V0'*sqrt(mu_tau(0)/mu_tau(1)),zeros([1, 3]),zeros([1, 3]),...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(eye(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(eye(3),[1,9]),...
    zeros([1, 3]), 0.5*sqrt(mu_tau(1)/mu_tau(0))*(1-mu_tau(0)/mu_tau(1))*V0',...
    zeros([1, 3]), zeros([1, 3])...
    )';
tspan=[0 dt*days2sec];
%tspan = t_nonlinear;
[t,y_initial] = ode113(@(t,y) internalIntegration(t,y,dUdr,ddUdrdr,jac_ddUdrdr,mu_tau,0),tspan,y0,options);
%plot(y(:,1),y(:,2));
rf_0 = y_initial(end,1:3)';
Vf_0 = y_initial(end,4:6)';
pvf_0 = y_initial(end,7:9)';
if case_traj == 1
    b=cat(1,rf_0,pvf_0)-cat(1,rf,[0 0 0]');
elseif case_traj == 2
    b=cat(1,rf_0,Vf_0)-cat(1,rf,Vf*sqrt(mu_tau(0)/mu_tau(1)));
end

%РѕРїС‚РёРјРёР·РёСЂСѓРµРј С‚СЂР°РµРєС‚РѕСЂРёСЋ
z0=zeros([1, 6]);
tic;
tauspan=linspace(0, 1,11);
[tau,z] = ode113(@(t,z) externalIntegration(t,z,b,dUdr,ddUdrdr,jac_ddUdrdr,y0,tspan,mu_tau,V0,Vf, case_traj),tauspan,z0,options);
evaluation_time = toc;
y0_final=y0;
y0_final(4:6)=V0;
y0_final(7:12)=z(end,:);
[t,y_final] = ode113(@(t,y) internalIntegration(t,y,dUdr,ddUdrdr,jac_ddUdrdr,mu_tau,1),tspan,y0_final,options);
%РљРѕРѕСЂРґРёРЅР°С‚С‹
rr_cont = y_final(:, 1:3);
%РЈСЃРєРѕСЂРµРЅРёРµ (pV)
a=y_final(:, 7:9);
a_vec=vecnorm(a, 2, 2).^2;
%Р¤СѓРЅРєС†РёРѕРЅР°Р»
Jt = cumtrapz(t, a_vec)/(2*eta);
%РњР°С‚СЂРёС†Р° С‡СѓРІСЃС‚РІРёС‚РµР»СЊРЅРѕСЃС‚Рё

drdpr=reshape(y_final(end,13:21),[3,3]);
drdpv=reshape(y_final(end,22:30),[3,3]);
drdz=cat(2,drdpr,drdpv);

ddrdprdt=reshape(y_final(end,31:39),[3,3]);
ddrdpvdt=reshape(y_final(end,40:48),[3,3]);
ddrdzdt=cat(2,ddrdprdt,ddrdpvdt);

dpvdpr=reshape(y_final(49:57),[3,3]);
dpvdpv=reshape(y_final(58:66),[3,3]);
dpvdz=cat(2,dpvdpr,dpvdpv);

drdtau = y_final(end,85:87)';
dvdtau = y_final(end,88:90)';
dpvdtau = y_final(end,91:93)';

if case_traj == 1
    dfdz = cat(1,drdz,dpvdz);
elseif case_traj == 2
    dfdz = cat(1,drdz,ddrdzdt);
end
C = cond(dfdz);
rf_cont=zeros(11,3);
for i = 1:11
    y0_final(4:6)=V0'*sqrt(mu_tau(tauspan(i))/mu_tau(1));
    y0_final(7:12)=z(i,:);
    [t,y_final] = ode113(@(t,y) internalIntegration(t,y,dUdr,ddUdrdr,jac_ddUdrdr,mu_tau,tauspan(i)),tspan,y0_final,options);
    rf_cont(i,:)=y_final(end,1:3)/ae;
end
[mars_r_f, mars_v_f]=planetEphemeris(tf,'SolarSystem',planet_end,'430');
mars_r_f=mars_r_f'*1e+03;
mars_v_f=mars_v_f'*1e+03;
dr = norm(rr_cont(end, :)'-mars_r_f);
dV = norm(y_final(end, 4:6)'-mars_v_f);



%Проверка "на глаз"
figure(1);
plot3(0, 0, 0, 'y--o')
set(gca,'FontSize',14)
hold on;

t_orbit = linspace(t0,t0+T_earth/(24*3600), 1000);
earth_traj = planetEphemeris(t_orbit','SolarSystem','Earth','430');
earth_traj=earth_traj*1e+03/ae;

t_orbit = linspace(t0,t0+T_mars/(24*3600), 1000);
mars_traj = planetEphemeris(t_orbit','SolarSystem',planet_end,'430');
mars_traj=mars_traj*1e+03/ae;
plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k')
plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'r')
plot3(rr_cont(:, 1)/ae, rr_cont(:, 2)/ae, rr_cont(:, 3)/ae, 'g', 'LineWidth', 2.5);

[mars_r_f, mars_v_f]=planetEphemeris([t0, dt],'SolarSystem',planet_end,'430');
mars_r_f=mars_r_f'*1e+03;
mars_v_f=mars_v_f'*1e+03;
plot3(mars_r_f(1)/ae, mars_r_f(2)/ae,mars_r_f(3)/ae,'rO')
plot3(rf_0(1)/ae, rf_0(2)/ae,rf_0(3)/ae,'bO')
plot3([rf_0(1),rf_0(1)-b(1)]/ae, [rf_0(2),rf_0(2)-b(2)]/ae,[rf_0(3),rf_0(3)-b(3)]/ae,'-k')
plot3(rf_cont(:,1),rf_cont(:,2),rf_cont(:,3),'b+');
%эти две точки должны находиться рядом
%plot3(rr_cont(500, 1)/ae, rr_cont(500, 2)/ae, rr_cont(500, 3)/ae, 'gO', 'LineWidth', 2.5);
%plot3(rr_old(500, 1), rr_old(500, 2), rr_old(500, 3), 'bO', 'LineWidth', 2.5);
axis equal

title('Траектория КА')
xlabel('x, a.e.')
ylabel('y, a.e.')
zlabel('z, a.e.')
view(0,90)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
hold off;