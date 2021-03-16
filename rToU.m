function u = rToU(r)
%UNTITLED Summary of this function goes here
%   ��������� ���������� ���������� � ��������������� 
x=r(1);
y=r(2);
z=r(3);

R1=sqrt(0.5*(norm(r)+x));
R2=sqrt(0.5*(norm(r)-x));
%������� ����� S � ������� C ���� theta, ��� theta=phi+gamma
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
%phi ������ ������� ��� ��������� u4 � ����
phi=0;
gamma = theta-phi;
u = [0 0 0 0]';
% ������� ��������������� ����������
u(1)=R1*cos(phi);
u(2)=R2*cos(gamma);
u(3)=R2*sin(gamma);
u(4)=R1*sin(phi);
end

