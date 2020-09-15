function t = integrateFunctional(s, y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
uu = y(:, 1:4);
t=zeros(length(uu),1);
for i = 1:length(uu)
    u = uu(i,:)';
    r=KS(u);
    v=y(i, 5:8)';
    h=y(i, 9)';
    tau=y(i ,10)';
    t(i) = tau-2*(u'*v)/(-2*h);
end
