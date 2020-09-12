function Jt = integrateFunctional(s, y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%summ_J = zeros(length(s),1);
%int_J = zeros(length(s),1);
u = y(:, 1:4);
r=zeros(length(u),4);
a=zeros(length(u),4);
t=zeros(length(u),1);
for i = 1:length(u)
    rr = u(i,:)';
    L = [[rr(1) -rr(2) -rr(3) rr(4)];
    [rr(2) rr(1) -rr(4) -rr(3)];
    [rr(3) rr(4) rr(1) rr(2)];
    [rr(4) -rr(3) rr(2) -rr(1)]];
    r(i,:)=L*rr;
    u2=norm(u)^2;
    v=y(i, 5:8)';
    h=y(i, 9)';
    tau=y(i ,10)';
    pu=y(i, 11:14)';
    pv=y(i, 15:18)';
    ph=y(i, 19)';
    ptau=y(i, 20)';
    aa=L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*(rr'*rr)*rr/(-2*h)^(3/2));
    
    La = [[aa(1) -aa(2) -aa(3) aa(4)];
    [aa(2) aa(1) -aa(4) -aa(3)];
    [aa(3) aa(4) aa(1) aa(2)];
    [aa(4) -aa(3) aa(2) -aa(1)]];
    a(i, :)=La*aa;
    t(i) = tau-2*(rr'*v)/(-2*h);
end
eta=0.45;
Jt = cumtrapz(t, vecnorm(a, 2, 2))/eta;

