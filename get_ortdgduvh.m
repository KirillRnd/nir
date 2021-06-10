function ortdgduvh = get_ortdgduvh(in1,in2,h)
%GET_ORTDGDUVH
%    ORTDGDUVH = GET_ORTDGDUVH(IN1,IN2,H)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    10-Jun-2021 19:34:48

u1 = in1(1,:);
u2 = in1(2,:);
u3 = in1(3,:);
u4 = in1(4,:);
v1 = in2(1,:);
v2 = in2(2,:);
v3 = in2(3,:);
v4 = in2(4,:);
t2 = u1.*v2;
t5 = u2.*v1;
t3 = t2-t5;
t4 = 1.0./t3;
t6 = t4.*u1.*u2;
ortdgduvh = reshape([t4.*u1.*u4,-t4.*u1.*u3,t6,-t4.*u1.^2,t4.*(u1.*v4-u4.*v1),-t4.*(u1.*v3-u3.*v1),1.0,0.0,0.0,t4.*u2.*u4,-t4.*u2.*u3,t4.*u2.^2,-t6,t4.*(u2.*v4-u4.*v2),-t4.*(u2.*v3-u3.*v2),0.0,1.0,0.0],[9,2]);
