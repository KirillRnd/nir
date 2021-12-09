function res = integrateTraectoryShperlingBurde(s, y)
%integrateTraectory интегрирует от начальной 
%точки до времени tf
%s - фиктивное время, tau - временной элемент
%Известные значения
%s
r=y(1:3);
Vs=y(4:6);
pX=y(7:9);
pVX=y(10:12);
tau=y(13);
A=y(14:16);
h=y(17);
dtds=norm(r)/sqrt(-2*h);
a=(pVX*(r'*r) + (Vs*Vs')*pVX)/dtds/(-2*h);
%Сохрняем провизводные
res=symFShperlingBurde(r,Vs,pX,pVX,A,h);
end

