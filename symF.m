function f = symF(in1,in2,h,in4,in5,ph)
%SYMF
%    F = SYMF(IN1,IN2,H,IN4,IN5,PH)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    11-Aug-2021 12:51:04

pu1 = in4(1,:);
pu2 = in4(2,:);
pu3 = in4(3,:);
pu4 = in4(4,:);
pw1 = in5(1,:);
pw2 = in5(2,:);
pw3 = in5(3,:);
pw4 = in5(4,:);
u1 = in1(1,:);
u2 = in1(2,:);
u3 = in1(3,:);
u4 = in1(4,:);
w1 = in2(1,:);
w2 = in2(2,:);
w3 = in2(3,:);
w4 = in2(4,:);
t2 = 1.0./h;
t3 = pw1.*t2.*w1;
t4 = pw2.*t2.*w2;
t5 = pw3.*t2.*w3;
t6 = pw4.*t2.*w4;
t13 = ph.*2.0;
t7 = t3+t4+t5+t6-t13;
t8 = u1.^2;
t9 = u2.^2;
t10 = u3.^2;
t11 = u4.^2;
t12 = t8+t9+t10+t11;
t14 = sqrt(2.0);
t15 = sqrt(-h);
t16 = t7.*w1;
t17 = pw1.*t2.*t12.*(1.0./4.0);
t18 = t16+t17;
t19 = t7.*w2;
t20 = pw2.*t2.*t12.*(1.0./4.0);
t21 = t19+t20;
t22 = t7.*w3;
t23 = pw3.*t2.*t12.*(1.0./4.0);
t24 = t22+t23;
t25 = t7.*w4;
t26 = pw4.*t2.*t12.*(1.0./4.0);
t27 = t25+t26;
t28 = 1.0./t12;
t29 = t18.*u1;
t30 = t27.*u4;
t47 = t21.*u2;
t48 = t24.*u3;
t31 = t29+t30-t47-t48;
t32 = 1.0./sqrt(-h);
t33 = t18.*u3;
t34 = t24.*u1;
t35 = t21.*u4;
t36 = t27.*u2;
t37 = t33+t34+t35+t36;
t38 = t18.*u2;
t39 = t21.*u1;
t53 = t24.*u4;
t54 = t27.*u3;
t40 = t38+t39-t53-t54;
t41 = u1.*w3.*2.0;
t42 = u3.*w1.*2.0;
t43 = u2.*w4.*2.0;
t44 = u4.*w2.*2.0;
t45 = t41+t42+t43+t44;
t46 = t14.*t15.*t28.*t37.*t45;
t49 = u1.*w1.*2.0;
t50 = u4.*w4.*2.0;
t60 = u2.*w2.*2.0;
t61 = u3.*w3.*2.0;
t51 = t49+t50-t60-t61;
t52 = t14.*t15.*t28.*t31.*t51;
t55 = u1.*w2.*2.0;
t56 = u2.*w1.*2.0;
t62 = u3.*w4.*2.0;
t63 = u4.*w3.*2.0;
t57 = t55+t56-t62-t63;
t58 = t14.*t15.*t28.*t40.*t57;
t59 = t46+t52+t58;
t64 = conj(t15);
t65 = pw1.*t9;
t66 = pw1.*t10;
t67 = pw1.*t11;
t68 = w1.^2;
t69 = pw1.*t68.*4.0;
t70 = pw2.*w1.*w2.*4.0;
t71 = pw3.*w1.*w3.*4.0;
t72 = pw4.*w1.*w4.*4.0;
t73 = pw1.*t8;
t106 = h.*ph.*w1.*8.0;
t74 = t65+t66+t67+t69+t70+t71+t72+t73-t106;
t75 = pw2.*t8;
t76 = pw2.*t9;
t77 = pw2.*t10;
t78 = pw2.*t11;
t79 = w2.^2;
t80 = pw2.*t79.*4.0;
t81 = pw1.*w1.*w2.*4.0;
t82 = pw3.*w2.*w3.*4.0;
t83 = pw4.*w2.*w4.*4.0;
t105 = h.*ph.*w2.*8.0;
t84 = t75+t76+t77+t78+t80+t81+t82+t83-t105;
t85 = pw3.*t8;
t86 = pw3.*t9;
t87 = pw3.*t10;
t88 = pw3.*t11;
t89 = w3.^2;
t90 = pw3.*t89.*4.0;
t91 = pw1.*w1.*w3.*4.0;
t92 = pw2.*w2.*w3.*4.0;
t93 = pw4.*w3.*w4.*4.0;
t107 = h.*ph.*w3.*8.0;
t94 = t85+t86+t87+t88+t90+t91+t92+t93-t107;
t95 = pw4.*t8;
t96 = pw4.*t9;
t97 = pw4.*t10;
t98 = pw4.*t11;
t99 = w4.^2;
t100 = pw4.*t99.*4.0;
t101 = pw1.*w1.*w4.*4.0;
t102 = pw2.*w2.*w4.*4.0;
t103 = pw3.*w3.*w4.*4.0;
t108 = h.*ph.*w4.*8.0;
t104 = t95+t96+t97+t98+t100+t101+t102+t103-t108;
t109 = t2.*t74.*u1.*(1.0./4.0);
t110 = t2.*t104.*u4.*(1.0./4.0);
t122 = t2.*t84.*u2.*(1.0./4.0);
t123 = t2.*t94.*u3.*(1.0./4.0);
t111 = t109+t110-t122-t123;
t112 = t2.*t74.*u2.*(1.0./4.0);
t113 = t2.*t84.*u1.*(1.0./4.0);
t120 = t2.*t94.*u4.*(1.0./4.0);
t121 = t2.*t104.*u3.*(1.0./4.0);
t114 = t112+t113-t120-t121;
t115 = t2.*t74.*u3.*(1.0./4.0);
t116 = t2.*t84.*u4.*(1.0./4.0);
t117 = t2.*t94.*u1.*(1.0./4.0);
t118 = t2.*t104.*u2.*(1.0./4.0);
t119 = t115+t116+t117+t118;
t124 = pw1.*u1.*u2.*2.0;
t125 = t111.^2;
t126 = t114.^2;
t127 = t119.^2;
t128 = t125+t126+t127;
t129 = 1.0./t12.^2;
t130 = pw1.*u1.*u3.*2.0;
t131 = pw1.*u2.*u3.*2.0;
t132 = pw2.*u1.*u4.*2.0;
t133 = pw2.*u2.*u4.*2.0;
t134 = pw3.*u3.*u4.*2.0;
t135 = pw3.*u1.*w3;
t136 = pw4.*u1.*w4;
t137 = pw3.*u2.*w3;
t138 = pw4.*u2.*w4;
t139 = pw1.*u1.*w1;
t140 = pw2.*u1.*w2;
t141 = pw2.*u3.*w2;
t142 = pw4.*u3.*w4;
t143 = pw1.*u4.*w1;
t144 = pw4.*u4.*w4;
t145 = pw1.*u2.*w1;
t146 = pw2.*u2.*w2;
t147 = pw1.*u3.*w1;
t148 = pw3.*u3.*w3;
t149 = pw2.*u4.*w2;
t150 = pw3.*u4.*w3;
t151 = 1.0./h.^2;
t152 = t65+t66+t67+t69+t70+t71+t72+t73;
t153 = t75+t76+t77+t78+t80+t81+t82+t83;
t154 = t85+t86+t87+t88+t90+t91+t92+t93;
t155 = t95+t96+t97+t98+t100+t101+t102+t103;
f = [w1;w2;w3;w4;u1.*(-1.0./4.0)+t2.*t59.*w1.*(1.0./2.0)-t14.*t31.*t32.*u1.*(1.0./4.0)-t14.*t32.*t37.*u3.*(1.0./4.0)-t14.*t32.*t40.*u2.*(1.0./4.0);u2.*(-1.0./4.0)+t2.*t59.*w2.*(1.0./2.0)+t14.*t31.*t32.*u2.*(1.0./4.0)-t14.*t32.*t37.*u4.*(1.0./4.0)-t14.*t32.*t40.*u1.*(1.0./4.0);u3.*(-1.0./4.0)+t2.*t59.*w3.*(1.0./2.0)+t14.*t31.*t32.*u3.*(1.0./4.0)-t14.*t32.*t37.*u1.*(1.0./4.0)+t14.*t32.*t40.*u4.*(1.0./4.0);u4.*(-1.0./4.0)+t2.*t59.*w4.*(1.0./2.0)-t14.*t31.*t32.*u4.*(1.0./4.0)-t14.*t32.*t37.*u2.*(1.0./4.0)+t14.*t32.*t40.*u3.*(1.0./4.0);-t14.*t28.*t31.*t51.*t64-t14.*t28.*t37.*t45.*t64-t14.*t28.*t40.*t57.*t64;-((t14.*t28.*t31.*t51.*t64+t14.*t28.*t37.*t45.*t64+t14.*t28.*t40.*t57.*t64).*(u1.*w1.*4.0+u2.*w2.*4.0+u3.*w3.*4.0+u4.*w4.*4.0)+t14.*t28.*t37.*t64.*(t12.*u1.*u3.*2.0+t12.*u2.*u4.*2.0)+t14.*t28.*t40.*t64.*(t12.*u1.*u2.*2.0-t12.*u3.*u4.*2.0)+t14.*t28.*t31.*t64.*(t8.*t12-t9.*t12-t10.*t12+t11.*t12)-1.0)./conj((h.*-2.0).^(3.0./2.0));pw1.*(1.0./4.0)-t14.*t15.*t28.*(t2.*t111.*(t65+t66+t67+t69+t70+t71+t72-t106+pw1.*t8.*3.0-pw2.*u1.*u2.*2.0-pw3.*u1.*u3.*2.0+pw4.*u1.*u4.*2.0).*(1.0./2.0)+t2.*t119.*(t86+t87+t88+t90+t91+t92+t93-t107+t130+t132+pw3.*t8.*3.0+pw4.*u1.*u2.*2.0).*(1.0./2.0)+t2.*t114.*(t76+t77+t78+t80+t81+t82+t83-t105+t124+pw2.*t8.*3.0-pw3.*u1.*u4.*2.0-pw4.*u1.*u3.*2.0).*(1.0./2.0)).*(1.0./2.0)+t14.*t15.*t128.*t129.*u1;pw2.*(1.0./4.0)-t14.*t15.*t28.*(t2.*t114.*(t66+t67+t69+t70+t71+t72+t73-t106+pw1.*t9.*3.0+pw2.*u1.*u2.*2.0-pw3.*u2.*u4.*2.0-pw4.*u2.*u3.*2.0).*(1.0./2.0)-t2.*t111.*(t75+t77+t78+t80+t81+t82+t83-t105-t124+pw2.*t9.*3.0+pw3.*u2.*u3.*2.0-pw4.*u2.*u4.*2.0).*(1.0./2.0)+t2.*t119.*(t95+t97+t98+t100+t101+t102+t103-t108+t131+t133+pw4.*t9.*3.0+pw3.*u1.*u2.*2.0).*(1.0./2.0)).*(1.0./2.0)+t14.*t15.*t128.*t129.*u2;pw3.*(1.0./4.0)+t14.*t15.*t28.*(t2.*t119.*(t65+t67+t69+t70+t71+t72+t73-t106+pw1.*t10.*3.0+pw3.*u1.*u3.*2.0+pw2.*u3.*u4.*2.0+pw4.*u2.*u3.*2.0).*(-1.0./2.0)+t2.*t111.*(t85+t86+t88+t90+t91+t92+t93-t107-t130+pw3.*t10.*3.0+pw2.*u2.*u3.*2.0-pw4.*u3.*u4.*2.0).*(1.0./2.0)+t2.*t114.*(t95+t96+t98+t100+t101+t102+t103-t108-t131+t134+pw4.*t10.*3.0-pw2.*u1.*u3.*2.0).*(1.0./2.0)).*(1.0./2.0)+t14.*t15.*t128.*t129.*u3;pw4.*(1.0./4.0)-t14.*t15.*t28.*(t2.*t119.*(t75+t76+t77+t80+t81+t82+t83-t105+pw2.*t11.*3.0+pw1.*u3.*u4.*2.0+pw3.*u1.*u4.*2.0+pw4.*u2.*u4.*2.0).*(1.0./2.0)-t2.*t114.*(t85+t86+t87+t90+t91+t92+t93-t107-t132+pw3.*t11.*3.0-pw1.*u2.*u4.*2.0+pw4.*u3.*u4.*2.0).*(1.0./2.0)+t2.*t111.*(t95+t96+t97+t100+t101+t102+t103-t108-t133-t134+pw4.*t11.*3.0+pw1.*u1.*u4.*2.0).*(1.0./2.0)).*(1.0./2.0)+t14.*t15.*t128.*t129.*u4;-pu1-t14.*t15.*t28.*(t2.*t111.*(t135+t136+t140-h.*ph.*u1.*2.0+pw1.*u1.*w1.*2.0-pw1.*u2.*w2-pw1.*u3.*w3+pw1.*u4.*w4).*2.0+t2.*t114.*(t137+t138+t146-h.*ph.*u2.*2.0+pw1.*u1.*w2+pw1.*u2.*w1.*2.0-pw1.*u3.*w4-pw1.*u4.*w3).*2.0+t2.*t119.*(t141+t142+t148-h.*ph.*u3.*2.0+pw1.*u1.*w3+pw1.*u3.*w1.*2.0+pw1.*u2.*w4+pw1.*u4.*w2).*2.0).*(1.0./2.0);-pu2-t14.*t15.*t28.*(t2.*t114.*(t135+t136+t139-h.*ph.*u1.*2.0+pw2.*u1.*w2.*2.0+pw2.*u2.*w1-pw2.*u3.*w4-pw2.*u4.*w3).*2.0-t2.*t111.*(t137+t138+t145-h.*ph.*u2.*2.0-pw2.*u1.*w1+pw2.*u2.*w2.*2.0+pw2.*u3.*w3-pw2.*u4.*w4).*2.0+t2.*t119.*(t143+t144+t150-h.*ph.*u4.*2.0+pw2.*u1.*w3+pw2.*u3.*w1+pw2.*u2.*w4+pw2.*u4.*w2.*2.0).*2.0).*(1.0./2.0);-pu3+t14.*t15.*t28.*(t2.*t119.*(t136+t139+t140-h.*ph.*u1.*2.0+pw3.*u1.*w3.*2.0+pw3.*u3.*w1+pw3.*u2.*w4+pw3.*u4.*w2).*-2.0+t2.*t111.*(t141+t142+t147-h.*ph.*u3.*2.0-pw3.*u1.*w1+pw3.*u2.*w2+pw3.*u3.*w3.*2.0-pw3.*u4.*w4).*2.0+t2.*t114.*(t143+t144+t149-h.*ph.*u4.*2.0-pw3.*u1.*w2-pw3.*u2.*w1+pw3.*u3.*w4+pw3.*u4.*w3.*2.0).*2.0).*(1.0./2.0);-pu4-t14.*t15.*t28.*(t2.*t119.*(t137+t145+t146-h.*ph.*u2.*2.0+pw4.*u1.*w3+pw4.*u3.*w1+pw4.*u2.*w4.*2.0+pw4.*u4.*w2).*2.0-t2.*t114.*(t141+t147+t148-h.*ph.*u3.*2.0-pw4.*u1.*w2-pw4.*u2.*w1+pw4.*u3.*w4.*2.0+pw4.*u4.*w3).*2.0+t2.*t111.*(t143+t149+t150-h.*ph.*u4.*2.0+pw4.*u1.*w1-pw4.*u2.*w2-pw4.*u3.*w3+pw4.*u4.*w4.*2.0).*2.0).*(1.0./2.0);t14.*t28.*t64.*(t111.*(t151.*t152.*u1.*(1.0./4.0)-t151.*t153.*u2.*(1.0./4.0)-t151.*t154.*u3.*(1.0./4.0)+t151.*t155.*u4.*(1.0./4.0)).*2.0+t114.*(t151.*t152.*u2.*(1.0./4.0)+t151.*t153.*u1.*(1.0./4.0)-t151.*t154.*u4.*(1.0./4.0)-t151.*t155.*u3.*(1.0./4.0)).*2.0+t119.*(t151.*t152.*u3.*(1.0./4.0)+t151.*t154.*u1.*(1.0./4.0)+t151.*t153.*u4.*(1.0./4.0)+t151.*t155.*u2.*(1.0./4.0)).*2.0).*(1.0./2.0)+(t14.*t28.*t128.*(1.0./4.0))./t64];
