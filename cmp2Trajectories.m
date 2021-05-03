function m = cmp2Trajectories(rr1,rr2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
DR=zeros(length(rr1),1);
for i = 1:length(rr1)
    r = rr1(i,:);
    dr = rr2-r;
    DR(i)=min(vecnorm(dr,2,2));
end
m = max(DR);
end

