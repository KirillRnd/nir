function res = integrateTraectoryCylindric(~, y)
%integrateTraectory интегрирует от начальной 
%точки до времени tf
%s - фиктивное время, tau - временной элемент
%Известные значения
%s
X=y(1:3);
VX=y(4:6);
pX=y(7:9);
pVX=y(10:12);
%Сохрняем провизводные
res=symFCylindric(X,VX,pX,pVX);
end

