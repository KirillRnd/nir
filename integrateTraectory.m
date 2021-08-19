function res = integrateTraectory(~, y)
%integrateTraectory интегрирует от начальной 
%точки до времени tf
%s - фиктивное время, tau - временной элемент
%Известные значения
%s
u=y(1:4);
w=y(5:8);

%u_alt = [u(4) -u(3) u(2) -u(1)]';
%v_til = u_alt*(u_alt'*v)/(u_alt'*u_alt);
%v_alt=norm(v)*(v-v_til)/norm(v-v_til);

%h=y(9);
tau=y(9);
pu=y(10:13);
pv=y(14:17);
%ph=y(19);
%ptau=y(20);
%ptau=0;
%Сохрняем провизводные
res=symF(u,w,pu,pv);
end

