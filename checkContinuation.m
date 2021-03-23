function [rr_cont] = checkContinuation(t0, dt, t_nonlinear, case_traj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mug = 132712.43994*(10^6)*(10^(3*3));
ae = 149597870700;
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;
days2sec=24*3600;
n = floor(dt/(T_earth/days2sec));
tf=t0+dt;

[r0, V0] = planetEphemeris(t0,'SolarSystem','Earth','430');
[rf, Vf] = planetEphemeris(tf,'SolarSystem','Mars','430');
r0=r0'*1e+03;
V0=V0'*1e+03;
rf=rf'*1e+03;
Vf=Vf'*1e+03;

% r0 = [ae 0 0]';
% V0 = [0 sqrt(mug/ae) 0]'; 
% a_m = 5*pi/3;
% rf = 1.52*ae*[cos(a_m) sin(a_m) 0]';
% Vf = sqrt(mug/(1.52*ae))*[cos(a_m+pi/2) sin(a_m+pi/2) 0]'; 
mu0 = defineMu0(n,r0,V0,rf,dt*days2sec,mug);
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
y0 = cat(2,r0',V0'*sqrt(mu_tau(0)/mu_tau(1)),[0 0 0],[0 0 0],...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(eye(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(eye(3),[1,9]),...
    zeros([1, 3]), 0.5*sqrt(mu_tau(1)/mu_tau(0))*(1-mu_tau(0)/mu_tau(1))*V0',...
    zeros([1, 3]), zeros([1, 3])...
    )';
%tspan=[0 dt*days2sec];
tspan = t_nonlinear;
[t,y_initial] = ode113(@(t,y) internalIntegration(t,y,dUdr,ddUdrdr,jac_ddUdrdr,mu_tau,0, case_traj),tspan,y0,options);
%plot(y(:,1),y(:,2));
rf_0 = y_initial(end,1:3)';
Vf_0 = y_initial(end,4:6)';
b=cat(1,rf_0,Vf_0)-cat(1,rf,Vf*sqrt(mu_tau(0)/mu_tau(1)));
%оптимизируем траекторию
z0=zeros([1, 6]);
tic;
[tau,z] = ode113(@(t,z) externalIntegration(t,z,b,dUdr,ddUdrdr,jac_ddUdrdr,y0,tspan,mu_tau,V0,Vf, case_traj),[0 1],z0,options);
toc
y0_final=y0;
y0_final(4:6)=V0;
y0_final(7:12)=z(end,:);
[t,y_final] = ode113(@(t,y) internalIntegration(t,y,dUdr,ddUdrdr,jac_ddUdrdr,mu_tau,1, case_traj),tspan,y0_final,options);
rr_cont = y_final(:, 1:3);
end

