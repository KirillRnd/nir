function res = optimiseToMars(tau, Z, b, tf)
%optimiseToMars ������������ ���������� 
%������� ����������� �� ���������
%   

%�������������� ��������
mug = 132712.43994*(10^6)*(10^(3*3));
%��������������� �������
ae = 149597870700;

%����� ��������� ��������� ��� �������������� � ����� tf
y0 = cat(2,[1*ae 0 0],[0 (mug/(1*ae))^(1/2) 0],Z',...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(eye(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(zeros(3),[1,9]),...
    reshape(zeros(3),[1,9]), reshape(eye(3),[1,9])...
    )';

%������������� ��������
optionsInn = odeset('AbsTol',1e-12);
[t,y] = ode45(@(t,y) integrateTraectory(t,y,mug),[0 tf],y0,optionsInn);

%��������� ������� ����������������
drdz=reshape(y(end,13:30),[3,6]);
drdzdt=reshape(y(end,49:66),[3,6]);

dfdz = cat(1,drdz,drdzdt);

tau
%��������� �����������
res=-(dfdz^-1)*b';
end

