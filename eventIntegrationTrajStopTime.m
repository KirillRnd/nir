function [value, isterminal, direction] = eventIntegrationTrajStopTime(s, y, time0, t_end, h0, t_start_fix)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

T_earth = 365.256363004*3600*24;
T_unit = T_earth/(2*pi);
u=y(1:4);
v=y(5:8);
h_end=y(9)+h0;
tau=y(10);
t_end_tmp = T_unit*(tau-2*(u'*v)/sqrt(-2*h_end))/(24*60*60)-t_start_fix;

value(1) = 10-toc(time0);
if value < 0
    value=0;
end

value(2) = t_end - t_end_tmp;
isterminal = [1 1];
direction = [0 0];
end