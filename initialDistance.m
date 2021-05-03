function dr = initialDistance(t_start,t_end,planet_start,planet_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

earth_r_f=planetEphemeris([t_start, t_end],'SolarSystem',planet_start,'430');
earth_r_f=earth_r_f'*1e+03;

mars_r_f=planetEphemeris([t_start, t_end],'SolarSystem',planet_end,'430');
mars_r_f=mars_r_f'*1e+03;

dr = norm(earth_r_f-mars_r_f);
end

