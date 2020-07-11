function res = integrateTraectory(s, y, symF)
%integrateTraectory интегрирует от начальной 
%точки до времени tf
%s - фиктивное время, tau - временной элемент
%Известные значения
%s
u=y(1:4);
v=y(5:8);
h=y(9);
pu=y(10:13);
pv=y(14:17);
ph=y(18);
tau=y(19);
ptau=y(20);
%Сохрняем провизводные
res=symF(h,ph,pu(1),pu(2),pu(3),pu(4),pv(1),pv(2),pv(3),pv(4),u(1),u(2),u(3),u(4),v(1),v(2),v(3),v(4));
end

