rad = 1/32;
phi = 11/8;

%������ ������� ����������. ��� ������ ��� ���������
a_check = phi*pi - rad*pi;
b_check = phi*pi + rad*pi;

%������ ������� ����������. ��� �� �������� ��� ����
a = (phi - rad)*pi;
b = (phi + rad)*pi;

%��� ��� ��������� ��������
a_check-a
b_check-b