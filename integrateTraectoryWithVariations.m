function res = integrateTraectoryWithVariations(~, Y, ~)
%integrateTraectory интегрирует от начальной 
%точки до времени tf
%s - фиктивное время, tau - временной элемент
%Известные значения
%s
res=zeros(420,1);
y=Y(1:20);
dydy0=reshape(Y(21:420), [20 20]);

u=y(1:4);
v=y(5:8);

%u_alt = [u(4) -u(3) u(2) -u(1)]';
%v_til = u_alt*(u_alt'*v)/(u_alt'*u_alt);
%v_alt=norm(v)*(v-v_til)/norm(v-v_til);

h=y(9);
tau=y(10);
pu=y(11:14);
pv=y(15:18);
ph=y(19);
ptau=y(20);

dfdy=symJ(u,v,h,pu,pv,ph,ptau);
ddudy0dt=dfdy*dydy0;
%Сохраняем провизводные
res(1:20)=symF(u,v,h,pu,pv,ph,ptau);
res(21:420)=reshape(ddudy0dt, [400 1]);
end

