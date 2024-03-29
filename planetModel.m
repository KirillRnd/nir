function [r,V] = planetModel(structure)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if strcmp(structure.mode,'Ephemeris')
    [r, V] = planetEphemeris(structure.t,'SolarSystem',structure.planet,'430');
elseif strcmp(structure.mode,'Flat')
    [r, V] = planetFlat(structure.t, structure.planet, structure.delta_omega);
end