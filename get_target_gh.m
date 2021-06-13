function g = get_target_gh(in1,in2,h)
%GET_TARGET_GH
%    G = GET_TARGET_GH(IN1,IN2,H)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    10-Jun-2021 19:34:47

u1 = in1(1,:);
u2 = in1(2,:);
u3 = in1(3,:);
u4 = in1(4,:);
v1 = in2(1,:);
v2 = in2(2,:);
v3 = in2(3,:);
v4 = in2(4,:);
t2 = sqrt(2.0);
t3 = sqrt(-h);
t4 = u1.^2;
t5 = u2.^2;
t6 = u3.^2;
t7 = u4.^2;
t8 = t4+t5+t6+t7;
t9 = 1.0./t8;
g = [t4-t5-t6+t7;u1.*u2.*2.0-u3.*u4.*2.0;u1.*u3.*2.0+u2.*u4.*2.0;t9.*(t2.*t3.*u1.*v1.*2.0-t2.*t3.*u2.*v2.*2.0-t2.*t3.*u3.*v3.*2.0+t2.*t3.*u4.*v4.*2.0);t9.*(t2.*t3.*u1.*v2.*2.0+t2.*t3.*u2.*v1.*2.0-t2.*t3.*u3.*v4.*2.0-t2.*t3.*u4.*v3.*2.0);t9.*(t2.*t3.*u1.*v3.*2.0+t2.*t3.*u3.*v1.*2.0+t2.*t3.*u2.*v4.*2.0+t2.*t3.*u4.*v2.*2.0);h];