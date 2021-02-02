t_start = juliandate(2001,0,0);
step = 1/16;
ds = 1/2:step:5/2;
rad = step/2;
L=length(ds);
DR=zeros([1,L]);
DV=zeros([1,L]);
CONV=zeros([1,L]);

warning('off');
for i=1:L
    i
    ds(i)
    [dr,dV, C] = checkMethod(t_start,ds(i),rad);
    DR(i)=dr;
    DV(i)=dV;
    CONV(i)=C;
end