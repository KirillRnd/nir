function [value, isterminal, direction] = eventIntegrationTraj(s, y, time0,maxtime)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

value = maxtime-toc(time0);
if value < 0
    value=0;
end
isterminal = 1;
direction = 0;
end

