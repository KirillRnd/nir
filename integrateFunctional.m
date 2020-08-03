function summ = integrateFunctional(s, y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
summ = 0;
for i = 2:length(s)
    step = s(i)-s(i-1);
    
    y_local=y(i,:)';
    u=y_local(1:4);
    v=y_local(5:8);
    h=y_local(9);
    pv=y_local(15:18);
    ph=y_local(19);
    ptau=y_local(20);
    %Вспомогательные величины
    L = [[u(1) -u(2) -u(3) u(4)];
    [u(2) u(1) -u(4) -u(3)];
    [u(3) u(4) u(1) u(2)];
    [u(4) -u(3) u(2) -u(1)]]; 
    u2=norm(u)^2;
    a_i =L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*(u'*u)*u/(-2*h)^(3/2));
    
    y_local=y(i-1,:)';
    u=y_local(1:4);
    v=y_local(5:8);
    h=y_local(9);
    pv=y_local(15:18);
    ph=y_local(19);
    ptau=y_local(20);
    %Вспомогательные величины
    L = [[u(1) -u(2) -u(3) u(4)];
    [u(2) u(1) -u(4) -u(3)];
    [u(3) u(4) u(1) u(2)];
    [u(4) -u(3) u(2) -u(1)]]; 
    u2=norm(u)^2;
    a_i_1 =L*(-(u2)*pv/(4*h) + v*(2*ph-(1/h)*pv'*v)+ptau*(u'*u)*u/(-2*h)^(3/2));
    
    trapeze = step*(0.5*norm(a_i_1)^2+0.5*norm(a_i)^2)/2;
    summ = summ + trapeze;
end

