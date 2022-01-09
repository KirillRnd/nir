function H_opt = calculateHamiltonian(in1,in2,in3,in4)
%calculateHamiltonian
%    H_opt = calculateHamiltonian(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    07-Jan-2022 17:23:45

pu1 = in3(1,:);
pu2 = in3(2,:);
pu3 = in3(3,:);
pu4 = in3(4,:);
pw1 = in4(1,:);
pw2 = in4(2,:);
pw3 = in4(3,:);
pw4 = in4(4,:);
u1 = in1(1,:);
u2 = in1(2,:);
u3 = in1(3,:);
u4 = in1(4,:);
w1 = in2(1,:);
w2 = in2(2,:);
w3 = in2(3,:);
w4 = in2(4,:);
t2 = u1.^2;
t3 = u1.^3;
t4 = u2.^2;
t5 = u2.^3;
t6 = u3.^2;
t7 = u3.^3;
t8 = u4.^2;
t9 = u4.^3;
t10 = w1.^2;
t11 = w2.^2;
t12 = w3.^2;
t13 = w4.^2;
t18 = sqrt(2.0);
t14 = t10.*4.0;
t15 = t11.*4.0;
t16 = t12.*4.0;
t17 = t13.*4.0;
t19 = t2+t4+t6+t8;
t20 = 1.0./t19;
t21 = t14+t15+t16+t17+t19;
t22 = t21.^2;
t23 = 1.0./t21;
t24 = sqrt(t23);
et1 = (t22.*(pw3.*t3+pw1.*t7+pw4.*t5+pw2.*t9+pw1.*t2.*u3+pw1.*t4.*u3+pw2.*t2.*u4+pw3.*t4.*u1+pw4.*t2.*u2+pw2.*t4.*u4+pw3.*t6.*u1+pw1.*t8.*u3+pw2.*t6.*u4+pw3.*t8.*u1+pw4.*t6.*u2+pw4.*t8.*u2+pw1.*t14.*u3+pw3.*t16.*u1+pw2.*t15.*u4+pw4.*t17.*u2+pw1.*u1.*w1.*w3.*4.0+pw1.*u2.*w1.*w4.*4.0+pw1.*u4.*w1.*w2.*4.0+pw2.*u1.*w2.*w3.*4.0+pw2.*u3.*w1.*w2.*4.0+pw2.*u2.*w2.*w4.*4.0+pw3.*u3.*w1.*w3.*4.0+pw3.*u2.*w3.*w4.*4.0+pw3.*u4.*w2.*w3.*4.0+pw4.*u1.*w3.*w4.*4.0+pw4.*u3.*w1.*w4.*4.0+pw4.*u4.*w2.*w4.*4.0).^2)./1.6e+1;
et2 = (t22.*(pw1.*t3-pw2.*t5-pw3.*t7+pw4.*t9+pw1.*t4.*u1-pw2.*t2.*u2+pw1.*t6.*u1-pw3.*t2.*u3+pw1.*t8.*u1-pw2.*t6.*u2-pw3.*t4.*u3+pw4.*t2.*u4-pw2.*t8.*u2+pw4.*t4.*u4-pw3.*t8.*u3+pw4.*t6.*u4-pw2.*t11.*u2.*4.0+pw1.*t14.*u1-pw3.*t12.*u3.*4.0+pw4.*t17.*u4-pw1.*u2.*w1.*w2.*4.0+pw2.*u1.*w1.*w2.*4.0-pw1.*u3.*w1.*w3.*4.0+pw3.*u1.*w1.*w3.*4.0+pw1.*u4.*w1.*w4.*4.0-pw2.*u3.*w2.*w3.*4.0-pw3.*u2.*w2.*w3.*4.0+pw4.*u1.*w1.*w4.*4.0+pw2.*u4.*w2.*w4.*4.0-pw4.*u2.*w2.*w4.*4.0+pw3.*u4.*w3.*w4.*4.0-pw4.*u3.*w3.*w4.*4.0).^2)./1.6e+1;
et3 = (t22.*(pw2.*t3+pw1.*t5-pw4.*t7-pw3.*t9+pw1.*t2.*u2+pw2.*t4.*u1+pw1.*t6.*u2+pw2.*t6.*u1-pw3.*t2.*u4-pw4.*t2.*u3+pw1.*t8.*u2+pw2.*t8.*u1-pw3.*t4.*u4-pw4.*t4.*u3-pw3.*t6.*u4-pw4.*t8.*u3+pw1.*t14.*u2+pw2.*t15.*u1-pw3.*t12.*u4.*4.0-pw4.*t13.*u3.*4.0+pw1.*u1.*w1.*w2.*4.0+pw2.*u2.*w1.*w2.*4.0-pw1.*u3.*w1.*w4.*4.0-pw1.*u4.*w1.*w3.*4.0+pw3.*u1.*w2.*w3.*4.0+pw3.*u2.*w1.*w3.*4.0-pw2.*u3.*w2.*w4.*4.0-pw2.*u4.*w2.*w3.*4.0+pw4.*u1.*w2.*w4.*4.0+pw4.*u2.*w1.*w4.*4.0-pw3.*u3.*w3.*w4.*4.0-pw4.*u4.*w3.*w4.*4.0).^2)./1.6e+1;
H_opt = t18.*t20.*t24.*(pw1.*u1.*(-1.0./4.0)-(pw2.*u2)./4.0-(pw3.*u3)./4.0-(pw4.*u4)./4.0+pu1.*w1+pu2.*w2+pu3.*w3+pu4.*w4+(t18.*t20.*t24.*(et1+et2+et3))./2.0);
