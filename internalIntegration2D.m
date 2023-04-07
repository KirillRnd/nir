function res = internalIntegration2D(~,y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
mug=1;

res=zeros(8,1);
r=y(1:2);
v=y(3:4);
pv=y(5:6);
dpvdt=y(7:8);

nr=norm(r);
res(1:2)=v;
res(3:4)=-mug*r/nr^3+pv;
res(5:6)=dpvdt;
res(7:8)=(r*(3*r'*pv/nr^2)-pv)*mug/nr^3;
end