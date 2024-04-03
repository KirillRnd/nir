function delta_rad = calculate_angular_distance(rr)
%calculate_angular_distance вычисляет угловую длину траектории 
% в формате Nx3. Возвращает радианы

angle_function = @(a,v)atan2(norm(cross(a, v)), dot(a,v));
an_distances = zeros(length(rr)-1,1);
for i=1:length(rr)-1
    r1 = rr(i,:);
    r2 = rr(i+1,:);
    an_distances(i) = angle_function(r1, r2);
end
delta_rad = sum(an_distances);
end