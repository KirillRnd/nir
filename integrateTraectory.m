function res = integrateTraectory(~, y)
%integrateTraectory ����������� �� ��������� 
%����� �� ������� tf
%s - ��������� �����, tau - ��������� �������
%��������� ��������
%s
u=y(1:4);
w=y(5:8);
pu=y(9:12);
pw=y(13:16);
tau=y(17);
%��������� ������������
res=symF(u,w,pu,pw);
%res(9:12)=-res(9:12);
%res(13:16)=-res(13:16);
end

