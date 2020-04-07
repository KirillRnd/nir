function res = optimiseToMars(tau, Z, b, tf, gamma,th_mars,v_opt,t_opt, y0)
%optimiseToMars Оптимизирует траекторию 
%методом продолжения по параметру
%   

%гравитационный параметр
mug = 132712.43994*(10^6)*(10^(3*3));

%задаём начальные параметры для интегрирования к точке tf
y0(7:12)= Z';

%устанавливаем точность
optionsInn = odeset('AbsTol',1e-12);
pv0=norm(Z(1:3));
[t,y] = ode45(@(t,y) integrateTraectory(t,y,mug, gamma,th_mars,v_opt,t_opt,pv0),[0 tf],y0,optionsInn);

%формируем матрицу чувствительности
drdz=reshape(y(end,13:30),[3,6]);
drdzdt=reshape(y(end,49:66),[3,6]);

dfdz = cat(1,drdz,drdzdt);

tau
%вычисляем производную
res=-(dfdz^-1)*b';
end

