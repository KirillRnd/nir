function m = massLP(Jt, m0, N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
m = zeros(length(Jt), 1);

for i=1:length(Jt)
    m(i)=m0/(1+Jt(i)*m0/N);
end    
end

