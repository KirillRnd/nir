function u = rToU(r,phi)
%UNTITLED Summary of this function goes here
%   Переводим физические координаты в параметрические 
x=r(1);
y=r(2);
z=r(3);

R1=sqrt(0.5*(norm(r)+x));
R2=sqrt(0.5*(norm(r)-x));
%Находим синус S и косинус C угла theta, где theta=phi+gamma
% C=y/(2*R1*R2);
% S=z/(2*R1*R2);
if z == 0 && y ==0
    theta=0;
else
    theta=atan2(z,y);
end

gamma = theta-phi;
u = [0 0 0 0]';
% находим параметрические координаты
u(1)=R1*cos(phi);
u(2)=R2*cos(gamma);
u(3)=R2*sin(gamma);
u(4)=R1*sin(phi);
end

