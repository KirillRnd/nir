function res = internalIntegration3D_withJ_TH(~,y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
mug=1;

res=zeros(13,1);
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

res(13)=(norm(pv)^2)/2; %dJdt
res(14)=(r(1)*v(2) - r(2) * v(1)) / (norm(r)^2); %dTHdt
end