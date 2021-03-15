function u = rToU(r,th)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Определим r = [1 0 0] как точку старта
x=r(1);
y=r(2);
z=r(3);

R1=sqrt(0.5*(norm(r)+x));
R2=sqrt(0.5*(norm(r)-x));
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

phi = th/2;
if (phi > pi/2) && (phi < 3*pi/2)
    j = -1;
else
    j = 1;
end

if (phi > pi) && (phi < 2*pi)
    i = -1;
else
    i = 1;
end


PHI=0;
gamma = theta-PHI;
u = [0 0 0 0]';

u(1)=R1*cos(PHI);
u(2)=R2*cos(gamma);
u(3)=R2*cos(gamma);
u(4)=R1*sin(PHI);


% u(4) = 0;
% u(1) = j*sqrt((norm(r)+x)/2);
% if abs(u(1)) < 1e+03
%     u(2)=i*y*sqrt((norm(r)-x)/2/(y^2+z^2));
%     u(3)=u(2)*z/y;
% else
% u(2) = i*y/(2*u(1));
% u(3) = z/(2*u(1));
% end
end

