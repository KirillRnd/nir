function res = integrateTraectory(y)
%integrateTraectory интегрирует от начальной 
%точки до времени tf
%s - фиктивное время, tau - временной элемент
%Известные значения
%s
r=y(1:3);
V=y(4:6);
pV=y(7:9);
dpVdt=y(10:12);

drdz=reshape(y(13:30),[3,6]);
ddrdzdt=reshape(y(49:66),[3,6]);

dpVdz=reshape(y(31:48),[3,6]);
ddpVdtdz=reshape(y(67:84),[3,6]);

res = symF(r,V,pV,dpVdt,drdz,ddrdzdt,dpVdz,ddpVdtdz);
end

