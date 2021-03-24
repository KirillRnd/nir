%��� ������ ���������� ������� ��������� � �������� �������� ������
t_start = juliandate(2022,0,0);
eta=0.45;
UorR='u';
step = 1/4;
ds = 2/2:step:3/2;
rad = step/2;
L=length(ds);
DR=zeros([1,L]);
DV=zeros([1,L]);
CONV=zeros([1,L]);
SF=zeros([1,L]);
PX=zeros([10,L]);
case_traj = 2;
x0=zeros([1, 12]);
warning('off');
planet_start = 'Earth';
planet_end = 'Mars';
%������������� ��� ������������� ���������
direction = -1;

for i=1:L
    %���� ���������� ������� �� ���������� ����� 4-� ������� ��� �������
    %������
    i
    ds(i)
    UorR = 'u';
    modifier_p=1e-06;
    modifier_f=1e+08;
    [dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks] = checkMethod(t_start,ds(i),rad,UorR,direction,modifier_p,modifier_f,x0,eta, case_traj,planet_end);
    DR(i)=dr;
    DV(i)=dV;
    CONV(i)=C;
    PX(:,i)=px;
    SF(i)=s_f;
%     %���������� �����������
%     if DR(i)>1e+07 && UorR == 'u'
%         direction = -1*direction
%         [dr,dV, C, px, sf] = checkMethod(t_start,ds(i),rad,UorR,direction,modifier_p,modifier_f,x0);
%         if dr<DR(i)
%             DR(i)=dr;
%             DV(i)=dV;
%             CONV(i) =C;
%             PX(:,i)=px;
%             SF(i)=sf;
%         else
%             direction = -1*direction;
%         end
%     end
% 
%     %������� ��������� � ��������� �����������
%     if DR(i)>1e+07 && UorR == 'u'
%         UorR = 'r';
%         [dr,dV, C, px,sf] = checkMethod(t_start,ds(i),rad,UorR,direction,modifier_p,modifier_f, x0);
%         if dr<DR(i)
%             DR(i)=dr;
%             DV(i)=dV;
%             CONV(i) =C;
%             PX(:,i)=px;
%             SF(i)=sf;
%         else
%              UorR = 'u';
%         end
%     end
    %������� ��������� �������
%     if DR(i)>1e+07 && UorR == 'u'
%         modifier_f=1e+06;
%         [dr,dV, C, px,sf] = checkMethod(t_start,ds(i),rad,UorR,direction,modifier_p,modifier_f,x0);
%         if dr<DR(i)
%             DR(i)=dr;
%             DV(i)=dV;
%             CONV(i) =C;
%             PX(:,i)=px;
%             SF(i)=sf;
%         else
%             modifier_f=1e+08;
%         end
%     end
    
end

%������ ������������ ����������


%t_end=t(end);
eul = [0 pi/4 0];
rotmZYX = eul2rotm(eul);
%�������� "�� ����"
figure(1);
plot3(0, 0, 0, 'y--o')
set(gca,'FontSize',14)
hold on;

% mars_traj = 1.52*[cos(th), sin(th), zeros(100,1)];
% earth_traj  = [cos(th), sin(th), zeros(100,1)];
t0 = t_start;
r_unit=ae;

mug_0 = 132712.43994*(10^6)*(10^(3*3));
V_unit=sqrt(mug_0/ae);
T_unit = T_earth/(2*pi);
ae = 149597870700;
T_mars=T_earth*1.8808476;
T_earth = 365.256363004*3600*24;
t_orbit = linspace(t0,t0+T_earth/(24*3600), 1000);
earth_traj = planetEphemeris(t_orbit','SolarSystem','Earth','430');
earth_traj=earth_traj*1e+03/ae;

earth_traj_New = arrayfun(@(x,y,z)rotmZYX*[x, y, z]', earth_traj(:, 1),earth_traj(:, 2),earth_traj(:, 3),'UniformOutput',false);
earth_traj_New = cell2mat(earth_traj_New')';

t_orbit = linspace(t0,t0+T_mars/(24*3600), 1000);
mars_traj = planetEphemeris(t_orbit','SolarSystem','Mars','430');
mars_traj=mars_traj*1e+03/ae;

mars_traj_New = arrayfun(@(x,y,z)rotmZYX*[x, y, z]', mars_traj(:, 1),mars_traj(:, 2),mars_traj(:, 3),'UniformOutput',false);
mars_traj_New = cell2mat(mars_traj_New')';

plot3(earth_traj(:, 1), earth_traj(:, 2), earth_traj(:, 3), 'k')
plot3(mars_traj(:, 1), mars_traj(:, 2), mars_traj(:, 3), 'r')





axis equal

%title('���������� ��')
xlabel('x, a.e.')
ylabel('y, a.e.')
zlabel('z, a.e.')

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
hold off;

%������� ���������� � ��������������� ����������"
figure(5);
plot3(0, 0, 0, 'y--o')
set(gca,'FontSize',14)
hold on;

th = linspace(0 ,4*pi,1000)';

mars_traj_ks = arrayfun(@(r1, r2, r3) rToU([r1,r2,r3], 0), mars_traj_New(:, 1),mars_traj_New(:, 2),mars_traj_New(:, 3),'UniformOutput',false);
mars_traj_ks = cell2mat(mars_traj_ks')';
mars_traj_ks=-mars_traj_ks*direction;
earth_traj_ks = arrayfun(@(r1, r2, r3) rToU([r1,r2,r3], 0), earth_traj_New(:, 1),earth_traj_New(:, 2),earth_traj_New(:, 3),'UniformOutput',false);
earth_traj_ks = cell2mat(earth_traj_ks')';
plot3(earth_traj_ks(:, 1), earth_traj_ks(:, 2), earth_traj_ks(:, 3), 'k')
plot3(mars_traj_ks(:, 1), mars_traj_ks(:, 2), mars_traj_ks(:, 3), 'r--')
plot3(-mars_traj_ks(:, 1), -mars_traj_ks(:, 2), -mars_traj_ks(:, 3), 'r--')

axis equal

title('���������� �� KS')
xlabel('u1')
ylabel('u2')
zlabel('u3')

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
hold off;

[r0, V0] = planetEphemeris(t_start,'SolarSystem',planet_start,'430');

eul = [0 pi/4 0];
rotmZYX = eul2rotm(eul);

r0 = [rotmZYX*r0'/ae; 0]*1e+03;
V0 = [rotmZYX*V0'/V_unit; 0]*1e+03;
mug=1;
h0=(norm(V0)^2)/2-mug/norm(r0);
t0=0;


u0 = rToU(r0, 0);

v0 = vFromV(V0,r0,mug, 0);
for i=1:L
    px = PX(:,i);
    s_f =SF(i);

    tau0= getEccentricAnomaly(r0(1:3),V0(1:3),mug);
    y0 = cat(1, u0, v0, 0, tau0,  px)';

    int_s0sf = linspace(0, s_f, 1e+03);
    %options = odeset('Events', @(s, y) eventIntegrationTraj(s, y, tf));
    options = odeset('AbsTol',1e-10);
    options = odeset(options,'RelTol',1e-10);
    %�����������, ��������� ����������� ���������� �� fmincon

    [s,y] = ode113(@(s,y) integrateTraectory(s,y,h0),int_s0sf, y0, options);


    uu = y(:, 1:4);
    vv = y(:, 5:8);
    rr=zeros(length(uu),4);
    a=zeros(length(uu),4);
    a_ks=zeros(length(uu),4);
    t=zeros(length(uu),1);
    VV=zeros(length(uu),4);
    t_start_fix=T_unit*(y(1, 10)-2*(y(1, 1:4)*y(1, 5:8)')/sqrt(-2*(y(1, 9)'+h0)))/(24*60*60);

    for j = 1:length(uu)
        u = uu(j,:)';
        r=KS(u);
        rr(j,:)=r;
        L_tmp=L_KS(u);
        u2=norm(u)^2;
        v=y(j, 5:8)';
        h=y(j, 9)'+h0;
        tau=y(j ,10)';
        pu=y(j, 11:14)';
        pv=y(j, 15:18)';
        ph=y(j, 19)';
        ptau=y(i, 20)';
        dtds=u2/sqrt(-2*h);
        aa_ks=L_tmp*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*(u2)*u/((-2*h)^(3/2)))/dtds;
        a_ks(j, :)=aa_ks/(ae/sqrt(mug_0)).^2;
        res=symF(u,v,h,pu,pv,ph,ptau);
        dvds=res(5:8);
        dhds=res(9);
        V = 2*sqrt(-2*h)*L_tmp*v/(u2);
        VV(j, :)=V;
        %a(i, :)=((-2*h/(norm(r)^2))*(2*(L_KS(v)*v+L_KS(u)*dvds)-(2*u'*v/(sqrt(-2*h)) + norm(r)*dhds/((-2*h)^(3/2)))*V)+mug*r/(norm(r)^3))/(ae/sqrt(mug_0)).^2;

        %a(i, :)=KS(aa);
        t(j) = T_unit*(tau-2*(u'*v)/sqrt(-2*h));
    end

    t_end = T_unit*(tau-2*(u'*v)/sqrt(-2*h))/(24*60*60)-t_start_fix;

    t = t - t(1);
    [rr_cont, Jt_cont] = checkContinuation(t_start, t_end, t, case_traj,planet_end,eta, floor(ds(i)));
    functional_cont = Jt_cont(end);

    figure(1);
    hold on;
    [mars_r_f, mars_v_f]=planetEphemeris([t_start, t_end],'SolarSystem',planet_end,'430');
    mars_r_f=mars_r_f'*1e+03;
    mars_v_f=mars_v_f'*1e+03;
    
    rr_old = arrayfun(@(x,y,z)rotmZYX^(-1)*[x, y, z]', rr(:, 1),rr(:, 2),rr(:, 3),'UniformOutput',false);
    rr_old = cell2mat(rr_old')';
    plot3(rr_old(:, 1), rr_old(:, 2), rr_old(:, 3), 'b', 'LineWidth', 2.5);

    plot3(rr_old(end, 1), rr_old(end, 2), rr_old(end, 3),'bO')
    plot3(mars_r_f(1)/ae, mars_r_f(2)/ae,mars_r_f(3)/ae,'rO')
    
    hold off;
    
    
    figure(5)
    hold on;
    plot3(uu(:, 1), uu(:, 2), uu(:, 3), 'b', 'LineWidth', 2.5);
    hold off;
    
     d = rr_old(:, 1:3)*ae-rr_cont;
     d_norm = vecnorm(d, 2, 2);
     disp(['������� ������� �  ����������� ', num2str(mean(d_norm),'%10.2e\n'), '�'])
end