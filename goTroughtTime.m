%��� ������ ���������� ������� ��������� � �������� �������� ������
t_start = juliandate(2022,0,0);
UorR='u';
step = 1/16;
ds = 1/2:step:4/2;
rad = step/2;
L=length(ds);
DR=zeros([1,L]);
DV=zeros([1,L]);
CONV=zeros([1,L]);

modifier_p=1e-04;
modifier_f=1e+04;

warning('off');
%������������� ��� ������������� ���������
direction = 1;
for i=1:L
    %���� ���������� ������� �� ���������� ����� 4-� ������� ��� �������
    %������
    i
    ds(i)
    modifier_p=1e-04;
    modifier_f=1e+04;
    UorR = 'u';
    
    [dr,dV, C] = checkMethod(t_start,ds(i),rad,UorR,direction,modifier_p,modifier_f);
    DR(i)=dr;
    DV(i)=dV;
    CONV(i)=C;
    
    
    %���������� �����������
    if DR(i)>1e+07 && UorR == 'u'
        direction = -1*direction
        [dr,dV, C] = checkMethod(t_start,ds(i),rad,UorR,direction,modifier_p,modifier_f);
        if dr<DR(i)
            DR(i)=dr;
            DV(i)=dV;
            CONV(i) =C;
        else
            direction = -1*direction;
        end
    end
    %������� ��������� � ��������� �����������
    if DR(i)>1e+07 && UorR == 'u'
        UorR = 'r';
        [dr,dV, C] = checkMethod(t_start,ds(i),rad,UorR,direction,modifier_p,modifier_f);
        if dr<DR(i)
            DR(i)=dr;
            DV(i)=dV;
            CONV(i) =C;
        else
             UorR = 'u';
        end
    end
    %������� ��������� �������
    if DR(i)>1e+07 && UorR == 'u'
        modifier_p=1e-06;
        [dr,dV, C] = checkMethod(t_start,ds(i),rad,UorR,direction,modifier_p,modifier_f);
        if dr<DR(i)
            DR(i)=dr;
            DV(i)=dV;
            CONV(i) =C;
        else
            modifier_p=1e-04;
        end
    end
    
end