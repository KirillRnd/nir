function res = integrateTraectory(t,y,mug, gamma,th_mars,v_opt,t_opt,pv0)
%integrateTraectory интегрирует от начальной 
%точки до времени tf

res=zeros(84,1);
r=y(1:3);
V=y(4:6);
pv=y(7:9);
Pv=y(10:12);

drdz=reshape(y(13:30),[3,6]);
drdzdt=reshape(y(49:66),[3,6]);

dpvdz=reshape(y(31:48),[3,6]);
dpvdzdt=reshape(y(67:84),[3,6]);

%Разбираем r - радиус вектор
x=r(1);
y=r(2);
z=r(3);
%Разбираем pv - сопряженную скорость
pvx=pv(1);
pvy=pv(2);
pvz=pv(3);

%Вспомогательные матрицы
ddUdrdr = mug*reshape([[ (3*x*x)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - 1/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2),(3*x*y)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2), (3*x*z)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2)]...
    [(3*y*x)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2), (3*y*y)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - 1/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2), (3*y*z)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2)]...
    [(3*z*x)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2), (3*z*y)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2), (3*z*z)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - 1/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2)]],[3,3]);
dddUdrdrpvdr=mug*reshape([[ conj(pvx)*((6*x)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*x^2*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)) + (3*y*conj(pvy))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*z*conj(pvz))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*x*y*conj(pvy)*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2) - (15*x*z*conj(pvz)*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2),                                                                                                         conj(pvx)*((3*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*x^2*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)) + (3*x*conj(pvy))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*x*y*conj(pvy)*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2) - (15*x*z*conj(pvz)*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2),                                                                                                         conj(pvx)*((3*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*x^2*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)) + (3*x*conj(pvz))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*x*y*conj(pvy)*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2) - (15*x*z*conj(pvz)*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)]...
[                                                                                                         conj(pvy)*((3*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*y^2*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)) + (3*y*conj(pvx))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*x*y*conj(pvx)*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2) - (15*y*z*conj(pvz)*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2), conj(pvy)*((6*y)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*y^2*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)) + (3*x*conj(pvx))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*z*conj(pvz))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*x*y*conj(pvx)*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2) - (15*y*z*conj(pvz)*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2),                                                                                                         conj(pvy)*((3*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*y^2*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)) + (3*y*conj(pvz))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*x*y*conj(pvx)*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2) - (15*y*z*conj(pvz)*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)]...
[                                                                                                         conj(pvz)*((3*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*z^2*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)) + (3*z*conj(pvx))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*x*z*conj(pvx)*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2) - (15*y*z*conj(pvy)*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2),                                                                                                         conj(pvz)*((3*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*z^2*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)) + (3*z*conj(pvy))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*x*z*conj(pvx)*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2) - (15*y*z*conj(pvy)*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2), conj(pvz)*((6*z)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*z^2*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)) + (3*x*conj(pvx))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*y*conj(pvy))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (15*x*z*conj(pvx)*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2) - (15*y*z*conj(pvy)*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)]],[3,3]);
 
%Вычисляем производные
res(1:3)=V;
nap=pv/norm(pv);

u_dv = 3000*9.81;
m0 = 1e+5;
flag_upr = 0;
if flag_upr == 1
    if m0-gamma*(t_opt+t)<=0
        res(4:6)=-mug*r/((norm(r))^3);
    else
        res(4:6)=-mug*r/((norm(r))^3)+(nap)*u_dv*gamma/(m0-gamma*(t_opt+t));
    end
    res(49:66)=reshape(ddUdrdr*drdz+dpvdz*(gamma*u_dv/(m0-gamma*(t_opt+t)))*(2+1/(norm(pv))),[18,1]);
else
    res(4:6)=-mug*r/((norm(r))^3)+pv;
    res(49:66)=reshape(ddUdrdr*drdz+dpvdz,[18,1]);
end

res(7:9)=Pv;
res(10:12)=ddUdrdr*pv;

res(13:30)=reshape(drdzdt,[18,1]);

res(31:48)=reshape(dpvdzdt,[18,1]);
res(67:84)=reshape(dddUdrdrpvdr*drdz+ddUdrdr*dpvdz,[18,1]);
end

