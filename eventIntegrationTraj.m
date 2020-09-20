function [value, isterminal, direction] = eventIntegrationTraj(s, y, StartTime)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
time=clock-StartTime;
time=time(end)
value = time < 1.0;
isterminal = 1;
direction = 0;
end

