function Jt = integrateFunctional(s, y, symF, eta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%summ_J = zeros(length(s),1);
%int_J = zeros(length(s),1);
mug_0 = 132712.43994*(10^6)*(10^(3*3));
ae = 149597870700;
T_earth = 365.256363004*3600*24;
T_norm = T_earth/(2*pi);
mug=1;
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
    a(i, :)=((-2*h/(norm(r)^2))*(2*(L_KS(v)*v+L_KS(u)*dvds)-(2*u'*v/(sqrt(-2*h)) + norm(r)*dhds/((-2*h)^(3/2)))*V)+mug*r/(norm(r)^3))/(ae/sqrt(mug_0)).^2;
    %a(i, :)=KS(aa);
    t(i) = T_norm*(tau-2*(u'*v)/(-2*h));
end
%eta=0.45;
Jt = cumtrapz(t, vecnorm(a, 2, 2).^2)/(2*eta);

