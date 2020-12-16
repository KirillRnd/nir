function f = symF(in1,in2,in3,in4,in5,in6,in7,in8)
%SYMF
%    F = SYMF(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    16-Dec-2020 18:25:00

V1 = in2(1,:);
V2 = in2(2,:);
V3 = in2(3,:);
pV1 = in3(1,:);
pV2 = in3(2,:);
pV3 = in3(3,:);
r1 = in1(1,:);
r2 = in1(2,:);
r3 = in1(3,:);
t2 = abs(r1);
t3 = abs(r2);
t4 = abs(r3);
t5 = t2.^2;
t6 = t3.^2;
t7 = t4.^2;
t8 = t5+t6+t7;
t9 = 1.0./t8.^(3.0./2.0);
f = [V1;V2;V3;pV1-r1.*t9.*1.3271243994e20;pV2-r2.*t9.*1.3271243994e20;pV3-r3.*t9.*1.3271243994e20;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0];
