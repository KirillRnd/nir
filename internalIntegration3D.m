function res = internalIntegration3D(~,y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
mug=1;

res=zeros(12,1);
r=y(1:3);
v=y(4:6);
pr=y(7:9);
pv=y(10:12);

nr=norm(r);
res(1:3)=v;
res(4:6)=-mug*r/nr^3+pv;
res(7:9)=mug*pv/nr^3 - mug*(3*(r*r')/norm(r)^(5))*pv;
%res(7:9)=mug*ddUdrdr(r)*pv;
res(10:12)=-pr;
end