function u = rToU(r)
%UNTITLED Summary of this function goes here
%   Переводим физические координаты в параметрические 
x=r(1);
y=r(2);
z=r(3);

R1=sqrt(0.5*(norm(r)+x));
R2=sqrt(0.5*(norm(r)-x));
%Находим синус S и косинус C угла theta, где theta=phi+gamma
C=y/(2*R1*R2);
S=z/(2*R1*R2);
theta=asin(S);
if C<0
    if S>0
        theta=pi-theta;
    elseif S<0
        theta=-pi-theta;
    else
        theta=theta+pi;
    end
end
%phi выбран нулевым для обращения u4 в ноль
phi=0;
gamma = theta-phi;
u = [0 0 0 0]';
% находим параметрические координаты
u(1)=R1*cos(phi);
u(2)=R2*cos(gamma);
u(3)=R2*sin(gamma);
u(4)=R1*sin(phi);
end

