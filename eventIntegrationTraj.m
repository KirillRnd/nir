function [value, isterminal, direction] = eventIntegrationTraj(s, y, rf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
u=y(1:2);
v=y(3:4);
h=y(5);
L = [[u(1) -u(2)];
    [u(2) u(1)]]; 

r=[(L*u)' 0];
%b=[cos(angle) sin(angle) 0];
b=[rf 0];
xa = r(1);
ya = r(2);
xb = b(1);
yb = b(2);
an = pi/2*((1+sign(xa))*(1-sign(ya^2))-(1+sign(xb))*(1-sign(yb^2)))...
    +pi/4*((2+sign(xa))*sign(ya)-(2+sign(xb))*sign(yb))...
    +sign(xa*ya)*atan((abs(xa)-abs(ya))/(abs(xa)+abs(ya)))...
    -sign(xb*yb)*atan((abs(xb)-abs(yb))/(abs(xb)+abs(yb)));
%an = atan2(norm(cross(r,b)), dot(r,b));

value = an;
isterminal = 1;
direction = 0;
end

