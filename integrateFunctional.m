function summ = integrateFunctional(s, y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
summ = 0;
for i = 2:length(s)
    step = s(i)-s(i-1);
    
    y_local=y(i,:)';
    u=y_local(1:2);
    v=y_local(3:4);
    h=y_local(5);
    pv=y_local(8:9);
    ph=y_local(10);
    %Вспомогательные величины
    L = [[u(1) -u(2)];
        [u(2) u(1)]]; 
    u2=norm(u)^2;
    alpha=2*ph-pv'*v/h;
    a_i =L*(-(u2)*pv/(4*h)+v*alpha);
    
    y_local=y(i-1,:)';
    u=y_local(1:2);
    v=y_local(3:4);
    h=y_local(5);
    pv=y_local(8:9);
    ph=y_local(10);
    %Вспомогательные величины
    L = [[u(1) -u(2)];
        [u(2) u(1)]]; 
    u2=norm(u)^2;
    alpha=2*ph-pv'*v/h;
    a_i_1 =L*(-(u2)*pv/(4*h)+v*alpha);
    
    trapeze = step*(0.5*norm(a_i_1)^2+0.5*norm(a_i)^2)/2;
    summ = summ + trapeze;
end

