function res = odeTrajectoryKS_2d(tau, y, mug)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tau

u=y(1:2);
v=y(3:4);
h=y(5);

pu=y(6:7);
pv=y(8:9);
ph=y(10);

L=KS_L(u,2);

dudt=v;
dvdt=-u/4-v*pv'*v/h-(pv*norm(u)^2)/(4*h);
dhdt=2*v'*pv;

dpudt=u*(-2/(norm(u)^5)*norm(L*pv)^2+norm(pv)^2/(norm(u)^4))+...
    pv*(-1/4-u'*pv/(2*h));
dpvdt=-pu+2*pv*v'*pv/h-2*ph*pv;
dphdt=-pv'/(h^2)*((u'*u)*pv/4+v'*pv*v);

% dpudt=zeros(2,1);
% dpvdt=zeros(2,1);
% dphdt=0;

res=zeros(10,1);

res(1:2)=dudt;
res(3:4)=dvdt;
res(5)=dhdt;

res(6:7)=dpudt;
res(8:9)=dpvdt;
res(10)=dphdt;
res(10)
end

