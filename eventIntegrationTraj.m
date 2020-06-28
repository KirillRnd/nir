function [value, isterminal, direction] = eventIntegrationTraj(s, y, angle, n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
u=y(1:4);
v=y(5:8);
h=y(9);
%ÇÀÄÀ×À ÏĞÎË¨ÒÀ
tau=y(19);
t = tau-2*(u'*v)/(-2*h);
T_earth = 365.256363004*3600*24;
value = T_earth*(n*2*pi + angle)/(2*pi)-t;
isterminal = 1;
direction = 0;
end

