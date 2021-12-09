function res = integrateTraectorySimpleTimeTransformation(~, y)
%integrateTraectory интегрирует от начальной 
%точки до времени tf
%s - фиктивное время, tau - временной элемент
%Известные значения
%s
r=y(1:3);
Vs=y(4:6);
prs=y(7:9);
pVs=y(10:12);
t=y(13);
%Сохрняем провизводные
res=symFSimpleTimeTransformation(r,Vs,prs,pVs);
end

