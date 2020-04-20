function [value, isterminal, direction] = eventIntegrationTraj(s, y, angle, n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
u=y(1:2);
v=y(3:4);
h=y(5);
%ÇÀÄÀ×À ÏĞÎË¨ÒÀ
tau=y(11);
t = tau-2*(u'*v)/(-2*h);
T_earth = 365.256363004*3600*24;
value = T_earth*(n*2*pi + angle)/(2*pi)-t;
isterminal = 1;
direction = 0;
end

