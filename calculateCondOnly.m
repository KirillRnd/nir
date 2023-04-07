function C = calculateCondOnly(t_start, px, s_f, phi)
%UNTITLED9 Summary of this function goes here
%   Вычисляет невязку в зависимости от входных параметров
%условия на fmincon
%ЗАДАЧА ПРОЛЁТА case_traj=1; ЗАДАЧА сопровождения case_traj=2;

%Начальные условия
ae = 149597870700;
mug_0 = 132712.43994*(10^6)*(10^(3*3));

V_unit=sqrt(mug_0/ae);
planet_start = 'Earth';
[r0, V0] = planetEphemeris(t_start,'SolarSystem',planet_start,'430');

eul = [0 pi/4 0];
rotmZYX = eul2rotm(eul);

r0 = [rotmZYX*r0'/ae; 0]*1e+03;
V0 = [rotmZYX*V0'/V_unit; 0]*1e+03;

mug=1;


%домножаем на коэффициент 1е-12, чтобы fmincon работал с более крупными
%величинами и не выдавал лишних ворнингов

%задаем начальные условия
%options = optimoptions(options,'OutputFcn',@myoutput);
%options = optimoptions(options, 'Algorithm', 'sqp');
phi0=phi;
h0=(norm(V0)^2)/2-mug/norm(r0);

u0 = rToU(r0, phi0);
w0 = vFromV(V0,r0,mug,phi0);

%tau0=getEccentricAnomaly(r0(1:3),V0(1:3),mug);
%tau0=0;
tau0=2*u0'*w0/sqrt(-2*h0);
y0 = cat(1, u0, w0, px', tau0)';

%t_start_fix=T_unit*(y0(10)-2*(y0(1:4)*y0(5:8)')/sqrt(-2*(y0(9)')))/(24*60*60);
int_s0sf = linspace(0, s_f, 1e+3);

options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);  

ddeltady0=zeros([8,8]);
step_h=1e-6;
for i=9:16
    y0_delta=zeros([1,17]);
    y0_delta(i)=step_h;
    [s,y] = ode113(@(s,y) integrateTraectory(s,y),int_s0sf,y0+y0_delta, options);
    u_end=y(end,1:4)';
    w_end=y(end,5:8)';
    px_end=y(end,9:16)';
    g=get_target_g(u_end,w_end);
    ortdgduv=get_ortdgduv(u_end,w_end);
    tr=[px_end'*ortdgduv]';
    p_plus=[g;tr];

    [s,y] = ode113(@(s,y) integrateTraectory(s,y),int_s0sf,y0-y0_delta, options);
    u_end=y(end,1:4)';
    w_end=y(end,5:8)';
    px_end=y(end,9:16)';
    g=get_target_g(u_end,w_end);
    ortdgduv=get_ortdgduv(u_end,w_end);
    tr=[px_end'*ortdgduv]';
    p_minus=[g;tr];
    partial=(p_plus-p_minus)/(2*step_h);
    ddeltady0(:,i-8)=partial;
end

C=cond(ddeltady0);

end

