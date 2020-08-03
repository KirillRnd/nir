function [value, isterminal, direction] = eventIntegrationTraj(s, y, tf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
u=y(1:4);
v=y(5:8);
h=y(9);
%ÇÀÄÀ×À ÏĞÎË¨ÒÀ
tau=y(10);
t = tau-2*(u'*v)/(-2*h);
value = tf-t;
isterminal = 1;
direction = 0;
end

