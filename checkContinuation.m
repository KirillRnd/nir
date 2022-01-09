function [rr_cont, Jt, C, evaluation_time, dr, dV,PR,PV, t_cont] = checkContinuation(t0, dt, t_nonlinear, case_traj,planet_end, eta,n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mug_0 = 132712.43994*(10^6)*(10^(3*3));
ae = 149597870700;
V_unit=sqrt(mug_0/ae);
a_unit=(ae/sqrt(mug_0)).^2;
T_earth = 365.256363004*3600*24;
T_unit = T_earth/(2*pi);
days2sec=24*3600;

mug=1;

tf=t0+dt;

eul = [0 pi/4 0];
rotmZYX = eul2rotm(eul);

[r0, V0] = planetEphemeris(t0,'SolarSystem','Earth','430');
[rf, Vf] = planetEphemeris(tf,'SolarSystem',planet_end,'430');
r0=rotmZYX*r0'/ae*1e+03;
V0=rotmZYX*V0'/V_unit*1e+03;
rf=rotmZYX*rf'/ae*1e+03;
Vf=rotmZYX*Vf'/V_unit*1e+03;
% 
% phi = getAngleRoRf(r0,rf,V0);
% 
% n = floor(dt/(T_earth/days2sec) - phi/(2*pi));
% r0 = [ae 0 0]';
% V0 = [0 sqrt(mug/ae) 0]'; 
% a_m = 5*pi/3;
% rf = 1.52*ae*[cos(a_m) sin(a_m) 0]';
% Vf = sqrt(mug/(1.52*ae))*[cos(a_m+pi/2) sin(a_m+pi/2) 0]'; 
mu0 = defineMu0(n,r0,V0,rf,dt*days2sec/T_unit,mug);
mu_tau=@(tau)(mu0+(mug-mu0)*tau);
%mu_tau=@(tau)(mug+0*tau);
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);   


r = sym('r', [3 1],'real');
assume(norm(r),'positive');
sym_dUdr=-mug*r/norm(r)^3;

sym_ddUdrdr=jacobian(sym_dUdr, r);
pV = sym('pV', [3 1],'real');

ddUdrdrpV=sym_ddUdrdr*pV;
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
%tspan=[0 dt*days2sec];
tspan = t_nonlinear/T_unit;
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

%оптимизируем траекторию
z0=zeros([1, 6]);
tic;
[tau,z] = ode113(@(t,z) externalIntegration(t,z,b,dUdr,ddUdrdr,jac_ddUdrdr,y0,tspan,mu_tau,V0,Vf, case_traj),[0 1],z0,options);
evaluation_time = toc;
y0_final=y0;
y0_final(4:6)=V0;
y0_final(7:12)=z(end,:);
% pr0=-z(end,4:6)';
% pv0=z(end,1:3)';
%Отключим ненадолго весь метод
%y0_final(7:12)=zeros(1,6);
[t,y_final] = ode113(@(t,y) internalIntegration(t,y,dUdr,ddUdrdr,jac_ddUdrdr,mu_tau,1),tspan,y0_final,options);
%Координаты
PV=y_final(:, 7:9);
PR=-y_final(:, 10:12); %он же dpvdt
rr_cont_rot = ae*y_final(:, 1:3);
rr_cont = arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', rr_cont_rot(:, 1),rr_cont_rot(:, 2),rr_cont_rot(:, 3),'UniformOutput',false);
rr_cont = cell2mat(rr_cont')';
%Ускорение (pV)
a=y_final(:, 7:9)/a_unit;
a_vec=vecnorm(a, 2, 2).^2;
%Функционал
Jt = cumtrapz(t*T_unit, a_vec)/(2*eta);
%Матрица чувствительности
drdpv=reshape(y_final(end,13:21),[3,3]);
drddpvdt=reshape(y_final(end,22:30),[3,3]);
drdz=cat(2,drdpv,drddpvdt);

ddrdpvdt=reshape(y_final(end,49:57),[3,3]);
ddVdpvdt=reshape(y_final(end,58:66),[3,3]);
ddrdzdt=cat(2,ddrdpvdt,ddVdpvdt);

dpvdpv=reshape(y_final(31:39),[3,3]);
dPvdpv=reshape(y_final(40:48),[3,3]);
dpvdz=cat(2,dpvdpv,dPvdpv);

if case_traj == 1
    dfdz = cat(1,drdz,dpvdz);
elseif case_traj == 2
    dfdz = cat(1,drdz,ddrdzdt);
end
C = cond(dfdz);

[mars_r_f, mars_v_f]=planetEphemeris(tf,'SolarSystem',planet_end,'430');
mars_r_f=rotmZYX*mars_r_f'*1e+03;
mars_v_f=rotmZYX*mars_v_f'*1e+03;
dr = norm(ae*rr_cont(end, :)'-mars_r_f);
dV = norm(V_unit*y_final(end, 4:6)'-mars_v_f);

t_cont=t*T_unit;
end

