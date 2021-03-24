function [c,ceq] = ubOrtPv(x, u0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Из билинейного соотношения
u_b = [u0(4) -u0(3) u0(2) -u0(1)];
pv = x(5:8)';
c = [];
ceq = u_b*pv;
end

