function Jt = integrateFunctional(s, y, eta,h0)
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
a_ks=zeros(length(uu),4);
t=zeros(length(uu),1);
VV=zeros(length(uu),4);
for i = 1:length(uu)
    u = uu(i,:)';
    r=KS(u);
    rr(i,:)=r;
    L=L_KS(u);
    u2=norm(u)^2;
    v=y(i, 5:8)';
    h=y(i, 9)'+h0;
    tau=y(i ,10)';
    pu=y(i, 11:14)';
    pv=y(i, 15:18)';
    ph=y(i, 19)';
    ptau=y(i, 20)';
    %aa=L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*(rr'*rr)*rr/(-2*h)^(3/2));
    res=symF(u,v,h,pu,pv,ph,ptau);
    dvds=res(5:8);
    dhds=res(9);
    dtds=u2/sqrt(-2*h);
    V = 2*sqrt(-2*h)*L*v/(u2);
    VV(i, :)=V;
    aa_ks=L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*(u2)*u/((-2*h)^(3/2)));
    a_ks(i, :)=aa_ks/(ae/sqrt(mug_0)).^2;
    t(i) = T_norm*tau-2*((ae/sqrt(mug_0)).^2)*(u'*v)/(-2*h);
end
%eta=0.45;
Jt = cumtrapz(t, vecnorm(a_ks, 2, 2).^2)/(2*eta);

