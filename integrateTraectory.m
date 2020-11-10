function res = integrateTraectory(s, y)
%integrateTraectory ����������� �� ��������� 
%����� �� ������� tf
%s - ��������� �����, tau - ��������� �������
%��������� ��������
%s
u=y(1:4);
v=y(5:8);
h=y(9);
tau=y(10);
pu=y(11:14);
pv=y(15:18);
ph=y(19);
ptau=y(20);
%�������� ������������
res=symF(u,v,h,pu,pv,ph,ptau);
end

