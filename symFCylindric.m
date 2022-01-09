function f = symFCylindric(in1,in2,in3,in4)
%SYMFCYLINDRIC
%    F = SYMFCYLINDRIC(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    01-Sep-2021 21:04:14

V1 = in2(1,:);
V2 = in2(2,:);
V3 = in2(3,:);
pVX1 = in4(1,:);
pVX2 = in4(2,:);
pVX3 = in4(3,:);
pX1 = in3(1,:);
pX2 = in3(2,:);
pX3 = in3(3,:);
phi = in1(2,:);
rho = in1(1,:);
z = in1(3,:);
t2 = cos(phi);
t3 = sin(phi);
t4 = V2.^2;
t5 = rho.^2;
t6 = rho.^3;
t7 = z.^2;
t10 = 1.0./rho;
t8 = pVX1.*t2;
t9 = pVX1.*t3;
t11 = t5+t7;
t12 = pVX2.*t2.*t10;
t13 = pVX2.*t3.*t10;
t14 = -t13;
t15 = t11.^(5.0./2.0);
t16 = 1.0./t11.^(3.0./2.0);
t18 = t9+t12;
t17 = 1.0./t15;
t19 = t8+t14;
f = [V1;V2;V3;rho.*t4-rho.*t16+t2.*t19+t3.*t18;-t10.*(V1.*V2.*2.0-t2.*t18+t3.*t19);pVX3-t16.*z;-(t17.*(pVX1.*1.0./t10.^5.*2.0-pVX2.^2.*t15+pVX3.*t5.^2.*z.*3.0-pVX1.*t6.*t7+pVX1.*t4.*t6.*t15+V1.*V2.*pVX2.*rho.*t15.*2.0))./t6;0.0;-t17.*(-pVX3.*t5+pVX3.*t7.*2.0+pVX1.*rho.*z.*3.0);-pX1+V2.*pVX2.*t10.*2.0;-pX2-V2.*pVX1.*rho.*2.0+V1.*pVX2.*t10.*2.0;-pX3];