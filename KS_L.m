function L = KS_L(u,d)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if d==2
    L=[u(1) -u(2); u(2) u(1)];
end    
if d==3
    L=[u(1) -u(2) -u(3) u(4);
        u(2) u(1) -u(4) -u(3);
        u(3) u(4) u(1) u(2);
        u(4) -u(3) u(2) -u(1)];
end  
end

