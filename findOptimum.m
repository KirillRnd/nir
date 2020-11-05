%symbolic_Jacob;
mug = 132712.43994*(10^6)*(10^(3*3));
ae = 149597870700;

r0 = [0 1*ae 0 0]';
V0 = [-(mug/(1*ae))^(1/2) 0 0 0]';

u0 = [0 0 0 0]';
h0=(norm(V0)^2)/2-mug/norm(r0);

u0(4) = 0;
u0(1) = sqrt((norm(r0)+r0(1))/2);
u0(2) = r0(2)/(2*u0(1));
u0(3) = r0(3)/(2*u0(1));

L = L_KS(u0); 

v0 = L'*V0/(2*sqrt(-2*h0));

pu0=[0 0 0 0]'*1e-13;
pv0=[0 0 0 0]'*1e-12;
ph0=0';
pt0=0;
t0 = 0;
y0 = cat(1, u0, v0, h0, t0, pu0, pv0, ph0, pt0)';
%Определяем tf
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);
s_f=pi;
[s,y] = ode113(@(s,y) integrateTraectory(s, y, symF),[0 s_f],y0, options);

uu = y(:, 1:4);
rr=zeros(length(uu),4);
a=zeros(length(uu),4);

t=zeros(length(uu),1);
VV=zeros(length(uu),4);
for i = 1:length(uu)
    u = uu(i,:)';
    r=KS(u);
    rr(i,:)=r;
    L=L_KS(u);
    u2=norm(u)^2;
    v=y(i, 5:8)';
    h=y(i, 9)';
    tau=y(i ,10)';
    pu=y(i, 11:14)';
    pv=y(i, 15:18)';
    ph=y(i, 19)';
    ptau=y(i, 20)';
    %aa=L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*(rr'*rr)*rr/(-2*h)^(3/2));
    res=symF(h,ph,ptau,pu(1),pu(2),pu(3),pu(4),pv(1),pv(2),pv(3),pv(4),u(1),u(2),u(3),u(4),v(1),v(2),v(3),v(4));
    dvds=res(5:8);
    dhds=res(9);
    V = 2*sqrt(-2*h)*L*v/(u2);
    VV(i, :)=V;
    a(i, :)=(-2*h/(norm(r)^2))*(2*(L_KS(v)*v+L_KS(u)*dvds)-(2*u'*v/(sqrt(-2*h)) + norm(r)*dhds/((-2*h)^(3/2)))*V)+mug*r/(norm(r)^3);
    
    %a(i, :)=KS(aa);
    t(i) = tau-2*(u'*v)/(-2*h);
end

t_end=t(end);

%Проверка "на глаз"
figure(1);
plot(0, 0,'y--o')
set(gca,'FontSize',14)
hold on;
th = 0:pi/50:2*pi;
plot(cos(th),sin(th),'k');
plot(1.52*cos(th),1.52*sin(th),'r');
plot(rr(:, 1)./ae, rr(:, 2)./ae,'b', 'LineWidth', 1.5)
%a_scale=3e+10/mean(vecnorm(a, 2, 2));
a_scale=0;
d = 24*3600;
idxes=1;
for i=1:ceil(t(end)/d)
    ix = find(t>d*i*10, 1);
    idxes=[idxes, ix];
end    
for i = idxes
    plot([rr(i, 1), rr(i, 1)+a_scale*a(i, 1)]./ae, [rr(i, 2), rr(i, 2)+a_scale*a(i, 2)]./ae,'k')
end
plot(rr(end, 1)./ae, rr(end, 2)./ae,'bO')
plot(1.52*cos(angle_M), 1.52*sin(angle_M),'rO')
axis equal

%title('Траектория КА')
xlabel('x, a.e.')
ylabel('y, a.e.')

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
hold off;