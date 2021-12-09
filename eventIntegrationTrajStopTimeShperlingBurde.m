function [value, isterminal, direction] = eventIntegrationTrajStopTimeShperlingBurde(~, y, time0,maxtime, t_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mug=1;
T_earth = 365.256363004*3600*24;
T_unit = T_earth/(2*pi);
r=y(1:3);
Vs=y(4:6);
h=y(17);
%h=-mug*norm(r)/(r'*r+Vs'*Vs);
tau=y(13);
t_end_tmp = T_unit*(tau-(r'*Vs)/norm(r)/sqrt(-2*h))/(24*60*60);

value(1) = maxtime-toc(time0);
if value < 0
    value=0;
end

value(2) = t_end - t_end_tmp;
isterminal = [1 1];
direction = [0 0];
end