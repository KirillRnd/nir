function res = integrateTraectory(s, y, mug)
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
%Вспомогательные величины
L = [[u(1) -u(2) -u(3) u(4)];
    [u(2) u(1) -u(4) -u(3)];
    [u(3) u(4) u(1) u(2)];
    [u(4) -u(3) u(2) -u(1)]]; 
u2=norm(u)^2;
alpha=2*ph-pv'*v/h;
a =L*(-(u2)*pv/(4*h)+v*alpha);
LTa=L'*a;

dhds=2*v'*LTa;
du2Lpvdu=u2*[[pv(1) pv(2) pv(3) pv(4)];
        [-pv(2) pv(1) pv(4) -pv(3)];
        [-pv(3) -pv(4) pv(1) pv(2)];
        [pv(4) -pv(3) pv(2) -pv(1)]] + 2 * u *(L*pv)';
dLvdu=[[v(1) v(2) v(3) v(4)];
        [-v(2) v(1) v(4) -v(3)];
        [-v(3) -v(4) v(1) v(2)];
        [v(4) -v(3) v(2) -v(1)]];
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
res=zeros(19,1);
res(1:4)=v;
res(5:8)=-u/4-u2*LTa/(4*h)-v*dhds/(2*h);
res(9)=dhds;

res(10:13)=dpuds;
res(14:17)=dpvds;
res(18)=dphds;
res(19)=dtauds;
end

