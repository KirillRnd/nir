function [r,V] = planetModel(structure)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if strcmp(structure.mode,'Ephemeris')

    if strcmp(structure.planet,'Asteroid2015RE36')
        [r, V] = planetFlat(structure.t, structure.planet, structure.delta_omega);
    else
        [r, V] = planetEphemeris(structure.t,'SolarSystem',structure.planet,'430');
    end
elseif strcmp(structure.mode,'Flat')
    [r, V] = planetFlat(structure.t, structure.planet, structure.delta_omega);
elseif strcmp(structure.mode,'Flat_a')
    [r, V] = planetFlat_a(structure.t, structure.planet, structure.delta_omega, structure.a_rel);
elseif strcmp(structure.mode,'Simple')
    [r, V] = planetModelSimple(structure.t, structure.planet);
end