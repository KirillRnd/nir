function u = KS_straight(r,d)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if d==3
    u=zeros(4,1);
    u(3)=0;
    x=r(1);
    y=r(2);
    z=r(3);

    u(2)=sqrt((norm(r)-x)/2);
    u(1)=y/(2*u(2));
    u(4)=z/(2*u(2));
end
if d==2
    u=zeros(2,1);
    x=r(1);
    
    u(2)=sqrt((norm(r)-x)/2);
    u(1)=sqrt((norm(r)+x)/2);
end
end

