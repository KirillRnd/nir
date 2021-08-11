function res = integrateTraectory(~, y)
%integrateTraectory ����������� �� ��������� 
%����� �� ������� tf
%s - ��������� �����, tau - ��������� �������
%��������� ��������
%s
u=y(1:4);
w=y(5:8);

%u_alt = [u(4) -u(3) u(2) -u(1)]';
%v_til = u_alt*(u_alt'*v)/(u_alt'*u_alt);
%v_alt=norm(v)*(v-v_til)/norm(v-v_til);

h=y(9);
tau=y(10);
pu=y(11:14);
pv=y(15:18);
ph=y(19);
%ptau=y(20);
%ptau=0;
%�������� ������������
res=symF(u,w,h,pu,pv,ph);
end

