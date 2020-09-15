function res = integrateTraectory(s, y, symF, symF_a0, amax)
%integrateTraectory интегрирует от начальной 
%точки до времени tf
%s - фиктивное время, tau - временной элемент
%Известные значения
%s
u=y(1:4);
v=y(5:8);
h=y(9);
tau=y(10);
pu=y(11:14);
pv=y(15:18);
ph=y(19);
ptau=y(20);
%Сохрняем провизводные
L=L_KS(u);
u2=u'*u;
lambda = L*(-pv*u2/(4*h)+(2*ph-pv'*v/h)*v+ptau*u*u2/((-2*h)^(3/2)));
if norm(lambda) == 0
    res=symF_a0(h,ptau,pu(1),pu(2),pu(3),pu(4),pv(1),pv(2),pv(3),pv(4),u(1),u(2),u(3),u(4),v(1),v(2),v(3),v(4));
else
    res=symF(amax,h,ph,ptau,pu(1),pu(2),pu(3),pu(4),pv(1),pv(2),pv(3),pv(4),u(1),u(2),u(3),u(4),v(1),v(2),v(3),v(4));
end
end

