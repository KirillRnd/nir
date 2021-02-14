function u = rToU(r)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
u = [0 0 0 0]';
u(4) = 0;
u(1) = sqrt((norm(r)+r(1))/2);
u(2) = r(2)/(2*u(1));
u(3) = r(3)/(2*u(1));
end

