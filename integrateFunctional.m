function Jt = integrateFunctional(t, y, eta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%summ_J = zeros(length(s),1);
%int_J = zeros(length(s),1);
mug_0 = 132712.43994*(10^6)*(10^(3*3));
ae = 149597870700;
T_earth = 365.256363004*3600*24;
T_unit = T_earth/(2*pi);
a_unit=(sqrt(mug_0)/ae).^2;
mug=1;
uu = y(:, 1:4);
rr=zeros(length(uu),4);
a_ks=zeros(length(uu),4);
VV=zeros(length(uu),4);
for i = 1:length(uu)
    u = uu(i,:)';
    r=KS(u);
    rr(i,:)=r;
    L=L_KS(u);
    u2=norm(u)^2;
    w=y(i, 5:8)';
    h=-mug/(u'*u+4*w'*w);
    tau=y(i ,17)';
    pu=y(i, 9:12)';
    pw=y(i, 13:16)';
    V = 2*sqrt(-2*h)*L*w/(u2);
    VV(i, :)=V;
    dtds=u2/sqrt(-2*h);
    %aa_ks=L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*(u2)*u/((-2*h)^(3/2)))/dtds;
    aa_ks=a_reactive(u,w,pu,pw);
    a_ks(i, :)=aa_ks*a_unit;
end
%eta=0.45;
a_vec=vecnorm(a_ks(:,1:3), 2, 2).^2;
%nr = vecnorm(rr(:,1:3),2,2);
Jt = cumtrapz(t, a_vec)/(2*eta);
%Jt = cumtrapz(t, a_vec)/(2*eta);
end
