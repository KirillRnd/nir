function [value, isterminal, direction] = eventIntegrationTrajStopTimeSimple(~, y, time0,maxtime, t_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
T_earth = 365.256363004*3600*24;
T_unit = T_earth/(2*pi);
t_end_tmp=T_unit*y(13)/(24*60*60);

value(1) = maxtime-toc(time0);
if value < 0
    value=0;
end

value(2) = t_end - t_end_tmp;
isterminal = [1 1];
direction = [0 0];
end