function res = integrateTraectoryWithVariations(~, Y)
%integrateTraectory интегрирует от начальной 
%точки до времени tf
%s - фиктивное время, tau - временной элемент
%Известные значения
%s
res=zeros(306,1);
y=Y(1:17);
dydy0=reshape(Y(18:306), [17 17]);

%u_alt = [u(4) -u(3) u(2) -u(1)]';
%v_til = u_alt*(u_alt'*v)/(u_alt'*u_alt);
%v_alt=norm(v)*(v-v_til)/norm(v-v_til);

u=y(1:4);
w=y(5:8);
pu=y(9:12);
pw=y(13:16);
tau=y(17);

dfdy=symJ(u,w,pu,pw);
ddudy0dt=dfdy*dydy0;
%Сохраняем провизводные
res(1:17)=symF(u,w,pu,pw);
res(18:306)=reshape(ddudy0dt, [289 1]);
end

