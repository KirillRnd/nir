function [value, isterminal, direction] = eventIntegrationTraj(s, y, time0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

value = 10-toc(time0);
if value < 0
    value=0;
end
isterminal = 1;
direction = 0;
end

