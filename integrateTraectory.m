function res = integrateTraectory(tau, y, mug)
%integrateTraectory интегрирует от начальной 
%точки до времени tf

res=zeros(5,1);
u=y(1:2);
v=y(3:4);
h=y(5);

L = [[u(1) -u(2)];
    [u(2) u(1)]]; 

%Вычисляем производные
res(1:2)=v;
res(3:4)=-u/4-v*0/(2*h);
res(5)=0;
end

