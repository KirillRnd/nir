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



r0 = [1*ae 0 0 0]';
V0 = [0 (mug/(1*ae))^(1/2) 0 0]';

r=r0;
V=V0;
%Проверка "на глаз"
figure(1);
plot(0, 0,'y--o')
hold on;
th = 0:pi/50:2*pi;
%plot(ae*cos(th),ae*sin(th),'k');
plot(1.52*ae*cos(th),1.52*ae*sin(th),'r');
%plot(r(:, 1), r(:, 2),'b')
acc = (VV(2:end, :)-VV(1:end-1, :))./(t(2:end)-t(1:end-1));
for i=1:length(a)-1
    plot(r(1), r(2),'b.')
    dt=t(i+1)-t(i);
    
    r_1=V;
    V_1=-mug*r/norm(r)^3+a(i, :)';
    
    r=r+r_1*dt;
    V=V+V_1*dt;
end
%plot(r(end, 1), r(end, 2),'bO')
plot(1.52*ae*cos(angle_M), 1.52*ae*sin(angle_M),'rO')
axis equal

title('Траектория КА')
xlabel('x, м')
ylabel('y, м')

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

hold off;

figure(2)
plot(t(2:end), vecnorm(acc-a(2:end, :), 2, 2));
sum(vecnorm(acc-a(2:end, :), 2, 2))