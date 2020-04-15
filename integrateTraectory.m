function res = integrateTraectory(s, y, mug)
%integrateTraectory интегрирует от начальной 
%точки до времени tf
%s - фиктивное время, tau - временной элемент
%Известные значения
%s
u=y(1:2);
v=y(3:4);
h=y(5);
pu=y(6:7);
pv=y(8:9);
ph=y(10);
tau=y(11);
%Вспомогательные величины
L = [[u(1) -u(2)];
    [u(2) u(1)]]; 
u2=norm(u)^2;
alpha=2*ph-pv'*v/h;
a =L*(-(u2)*pv/(4*h)+v*alpha);
LTa=L'*a;

dhds=2*v'*LTa;
du2Lpvdu=2*u(1)*u(2)*[[-pv(2) pv(1)];[pv(1) pv(2)]]+...
    (3*u(1)^2+u(2)^2)*[[pv(1) pv(2)];[0 0]]+...
    (3*u(2)^2+u(1)^2)*[[0 0];[-pv(2) pv(1)]];
dLvdu=[[v(1) v(2)];[-v(2) v(1)]];
dLTadu=u2*(-pv*u'/(2*h))+2*(-u2*pv/(4*h)+v*alpha)*u';
%Вспомогательные производные управления
dadh=((u2*L*pv)/4+L*v*(pv'*v))/(h^2);
dadv=L*alpha+(L*v)*pv'/h;
dadu=-du2Lpvdu/(4*h)+dLvdu*alpha;
%Производные сопряженных переменных
dpuds=-(-a'*dadu+pv'*(-1/4-u2*dLTadu/(4*h)-LTa'*u/(2*h))+alpha*v'*dLTadu);
dpvds=-pu'+pv'*(u2*L'*dadv/(4*h)+(v'*LTa)/h)-alpha*(v'*L'*dadv+LTa');
dphds=-pv'*(u2*LTa/(4*(h^2))-u2*L'*dadh/(4*h))...
    -(pv'*v/(h^2))*v'*LTa-alpha*v'*(L'*dadh);
%Производная временнОго элемента
dtauds=(mug+4*u2*dhds+u2*u'*LTa)/((-2*h)^(3/2));
%Сохрняем провизводные
res=zeros(11,1);
res(1:2)=v;
res(3:4)=-u/4-u2*LTa/(4*h)-v*dhds/(2*h);
res(5)=dhds;

res(6:7)=dpuds;
res(8:9)=dpvds;
res(10)=dphds;
res(11)=dtauds;
end

