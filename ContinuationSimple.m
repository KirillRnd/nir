function [rr_cont, Jt, C, evaluation_time, dr, dV,PR,PV, t_cont] = ContinuationSimple(t0, dt, t_nonlinear,z0, case_traj,planet_end, eta,n, orbits, omega)
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

eul = [pi/12 pi/4 pi/12];
rotmZYX = eul2rotm(eul);

st.t = t0;
st.planet = 'Earth';
st.mode = orbits;
st.delta_omega = omega;

[r0, V0] = planetModel(st);

st.t = tf;
st.planet = planet_end;
st.mode = orbits;
st.delta_omega = omega;

[rf, Vf] = planetModel(st);
r0=rotmZYX*r0'/ae*1e+03;
V0=rotmZYX*V0'/V_unit*1e+03;
rf=rotmZYX*rf'/ae*1e+03;
Vf=rotmZYX*Vf'/V_unit*1e+03;

options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);   


y0 = cat(2,r0',V0',zeros([1, 3]),zeros([1, 3]),...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(eye(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(eye(3),[1,9]))';
y0(7:12)=z0;
tspan = t_nonlinear/T_unit;
[t,y_initial] = ode113(@(t,y) internalIntegrationSimple(t,y),tspan,y0,options);
%plot(y(:,1),y(:,2));
rf_0 = y_initial(end,1:3)';
Vf_0 = y_initial(end,4:6)';
pvf_0 = y_initial(end,7:9)';
if case_traj == 1
    b=cat(1,rf_0,pvf_0)-cat(1,rf,[0 0 0]');
elseif case_traj == 2
    b=cat(1,rf_0,Vf_0)-cat(1,rf,Vf);
end

%оптимизируем траекторию
%z0=zeros([1, 6]);
tic;
[tau,z] = ode113(@(tau,z) externalIntegrationSimple(tau,z,b,y0,tspan,case_traj),[0 1],z0,options);
evaluation_time = toc;
y0_final=y0;
y0_final(4:6)=V0;
y0_final(7:12)=z(end,:);
% pr0=-z(end,4:6)';
% pv0=z(end,1:3)';
%Отключим ненадолго весь метод
%y0_final(7:12)=zeros(1,6);
[t,y_final] = ode113(@(t,y) internalIntegrationSimple(t,y),tspan,y0_final,options);
%Координаты
PV=y_final(:, 7:9);
PR=-y_final(:, 10:12); %он же dpvdt


PV = arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', PV(:, 1),PV(:, 2),PV(:, 3),'UniformOutput',false);
PV = cell2mat(PV')';
PR = arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', PR(:, 1),PR(:, 2),PR(:, 3),'UniformOutput',false);
PR = cell2mat(PR')';

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

st.t = tf;
st.planet = planet_end;
st.mode = orbits;
st.delta_omega = omega;

[mars_r_f, mars_v_f]=planetModel(st);
mars_r_f=rotmZYX*mars_r_f'*1e+03;
mars_v_f=rotmZYX*mars_v_f'*1e+03;
dr = norm(ae*y_final(end, 1:3)'-mars_r_f);
dV = norm(V_unit*y_final(end, 4:6)'-mars_v_f);

t_cont=t*T_unit;
end

