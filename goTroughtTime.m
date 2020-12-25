t0=juliandate(2001,0,0);
ds = pi/4:pi/48:9*pi/4;

L=length(ds);
DR=zeros([1,L]);
DV=zeros([1,L]);
CONV=zeros([1,L]);

warning('off');
for i=1:L
    ds(i)
    [dr,dV, C] = checkMethod(t0,ds(i));
    DR(i)=dr;
    DV(i)=dV;
    CONV(i)=C;
end