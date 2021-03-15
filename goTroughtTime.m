%t_start = juliandate(2001,0,0);
t_Mars_0=0.25;
step = 1/16;
ds = 1/2:step:10/2;
rad = step/2;
L=length(ds);
DR=zeros([1,L]);
DV=zeros([1,L]);
CONV=zeros([1,L]);

warning('off');
for i=1:L
    i
    ds(i)
    [dr,dV, C] = checkMethod(t_Mars_0,ds(i),rad);
    DR(i)=dr;
    DV(i)=dV;
    CONV(i)=C;
end