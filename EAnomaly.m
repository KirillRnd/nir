function E = EAnomaly(e, M)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
E=M;
while E-e*sin(E)-M > 1e-4
    E=e*sin(E)+M;
end
end

