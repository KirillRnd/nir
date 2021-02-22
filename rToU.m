function u = rToU(r,th)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Определим r = [1 0 0] как точку старта
x=r(1);
y=r(2);
z=r(3);

phi = th/2;
if (phi > pi/2) && (phi < 3*pi/2)
    j = -1;
else
    j = 1;
end

if (phi > pi) && (phi < 2*pi)
    i = -1;
else
    i = 1;
end

u = [0 0 0 0]';
u(4) = 0;
u(1) = j*sqrt((norm(r)+x)/2);
if abs(u(1)) < 1e-02
    u(2)=i*y*sqrt((norm(r)-x)/2/(y^2+z^2));
    u(3)=u(2)*z/y;
else
u(2) = y/(2*u(1));
u(3) = z/(2*u(1));
end
end

