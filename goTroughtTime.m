t0=juliandate(2001,0,0);
DR=zeros([1,49]);
DV=zeros([1,49]);
CONV=zeros([1,49]);
ds = pi/4:pi/48:9*pi/4;
warning('off');
for i=1:49
    ds(i)
    [dr,dV, C] = checkMethod(t0,ds(i));
    DR(i)=dr;
    DV(i)=dV;
    CONV(i)=C;
end