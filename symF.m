function f = symF(in1,in2,h,in4,in5,ph,ptau)
%SYMF
%    F = SYMF(IN1,IN2,H,IN4,IN5,PH,PTAU)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    12-Jun-2021 20:04:26

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
t3 = u1.^2;
t4 = u2.^2;
t5 = u3.^2;
t6 = u4.^2;
t7 = t3+t4+t5+t6;
t8 = pw1.*t2.*w1;
t9 = pw2.*t2.*w2;
t10 = pw3.*t2.*w3;
t11 = pw4.*t2.*w4;
t19 = ph.*2.0;
t12 = t8+t9+t10+t11-t19;
t13 = u1.*w1.*8.0;
t14 = u2.*w2.*8.0;
t15 = u3.*w3.*8.0;
t16 = u4.*w4.*8.0;
t17 = t13+t14+t15+t16;
t20 = h.*2.0;
t18 = 1.0./(-t20).^(3.0./2.0);
t21 = sqrt(2.0);
t22 = sqrt(-h);
t23 = t12.*w1;
t24 = t17.*w1;
t25 = t7.*u1;
t26 = t24+t25;
t27 = pw1.*t2.*t7.*(1.0./4.0);
t48 = ptau.*t18.*t26;
t28 = t23+t27-t48;
t29 = t12.*w2;
t30 = t17.*w2;
t31 = t7.*u2;
t32 = t30+t31;
t33 = pw2.*t2.*t7.*(1.0./4.0);
t49 = ptau.*t18.*t32;
t34 = t29+t33-t49;
t35 = t12.*w3;
t36 = t17.*w3;
t37 = t7.*u3;
t38 = t36+t37;
t39 = pw3.*t2.*t7.*(1.0./4.0);
t50 = ptau.*t18.*t38;
t40 = t35+t39-t50;
t41 = t12.*w4;
t42 = t17.*w4;
t43 = t7.*u4;
t44 = t42+t43;
t45 = pw4.*t2.*t7.*(1.0./4.0);
t51 = ptau.*t18.*t44;
t46 = t41+t45-t51;
t47 = 1.0./t7;
t52 = t28.*u1;
t53 = t46.*u4;
t70 = t34.*u2;
t71 = t40.*u3;
t54 = t52+t53-t70-t71;
t55 = 1.0./sqrt(-h);
t56 = t28.*u3;
t57 = t34.*u4;
t58 = t40.*u1;
t59 = t46.*u2;
t60 = t56+t57+t58+t59;
t61 = t28.*u2;
t62 = t34.*u1;
t76 = t40.*u4;
t77 = t46.*u3;
t63 = t61+t62-t76-t77;
t64 = u1.*w3.*2.0;
t65 = u3.*w1.*2.0;
t66 = u2.*w4.*2.0;
t67 = u4.*w2.*2.0;
t68 = t64+t65+t66+t67;
t69 = t21.*t22.*t47.*t60.*t68;
t72 = u1.*w1.*2.0;
t73 = u4.*w4.*2.0;
t83 = u2.*w2.*2.0;
t84 = u3.*w3.*2.0;
t74 = t72+t73-t83-t84;
t75 = t21.*t22.*t47.*t54.*t74;
t78 = u1.*w2.*2.0;
t79 = u2.*w1.*2.0;
t85 = u3.*w4.*2.0;
t86 = u4.*w3.*2.0;
t80 = t78+t79-t85-t86;
t81 = t21.*t22.*t47.*t63.*t80;
t82 = t69+t75+t81;
t87 = (-t20).^(3.0./2.0);
t88 = conj(t87);
t89 = 1.0./t88;
t90 = conj(t22);
t95 = ptau.*t26.*t89;
t91 = t23+t27-t95;
t96 = ptau.*t32.*t89;
t92 = t29+t33-t96;
t97 = ptau.*t38.*t89;
t93 = t35+t39-t97;
t98 = ptau.*t44.*t89;
t94 = t41+t45-t98;
t99 = t91.*u3;
t100 = t92.*u4;
t101 = t93.*u1;
t102 = t94.*u2;
t103 = t99+t100+t101+t102;
t104 = t91.*u1;
t105 = t94.*u4;
t112 = t92.*u2;
t113 = t93.*u3;
t106 = t104+t105-t112-t113;
t107 = t91.*u2;
t108 = t92.*u1;
t110 = t93.*u4;
t111 = t94.*u3;
t109 = t107+t108-t110-t111;
t114 = pw1.*w1;
t115 = pw2.*w2;
t116 = pw3.*w3;
t117 = pw4.*w4;
t122 = h.*ph.*2.0;
t118 = t114+t115+t116+t117-t122;
t119 = (-h).^(3.0./2.0);
t120 = abs(h);
t121 = 1.0./t120.^3;
t123 = w1.^2;
t124 = t2.*t118.*w2;
t125 = 1.0./(-h).^(5.0./2.0);
t126 = 1.0./(-h).^(3.0./2.0);
t127 = t3.*u2;
t128 = t5.*u2;
t129 = t6.*u2;
t130 = w2.^2;
t131 = t130.*u2.*8.0;
t132 = u1.*w1.*w2.*8.0;
t133 = u3.*w2.*w3.*8.0;
t134 = u4.*w2.*w4.*8.0;
t135 = t127+t128+t129+t131+t132+t133+t134+t4.*u2;
t136 = t2.*t118.*w1;
t137 = t4.*u1;
t138 = t5.*u1;
t139 = t6.*u1;
t140 = t123.*u1.*8.0;
t141 = u2.*w1.*w2.*8.0;
t142 = u3.*w1.*w3.*8.0;
t143 = u4.*w1.*w4.*8.0;
t144 = t137+t138+t139+t140+t141+t142+t143+t3.*u1;
t204 = ptau.*t21.*t119.*t121.*t144.*(1.0./4.0);
t145 = t27+t136-t204;
t185 = ptau.*t21.*t119.*t121.*t135.*(1.0./4.0);
t146 = t33+t124-t185;
t147 = t2.*t118.*w3;
t148 = t3.*u3;
t149 = t4.*u3;
t150 = t6.*u3;
t151 = w3.^2;
t152 = t151.*u3.*8.0;
t153 = u1.*w1.*w3.*8.0;
t154 = u2.*w2.*w3.*8.0;
t155 = u4.*w3.*w4.*8.0;
t156 = t148+t149+t150+t152+t153+t154+t155+t5.*u3;
t203 = ptau.*t21.*t119.*t121.*t156.*(1.0./4.0);
t157 = t39+t147-t203;
t158 = t2.*t118.*w4;
t159 = t3.*u4;
t160 = t4.*u4;
t161 = t5.*u4;
t162 = w4.^2;
t163 = t162.*u4.*8.0;
t164 = u1.*w1.*w4.*8.0;
t165 = u2.*w2.*w4.*8.0;
t166 = u3.*w3.*w4.*8.0;
t167 = t159+t160+t161+t163+t164+t165+t166+t6.*u4;
t205 = ptau.*t21.*t119.*t121.*t167.*(1.0./4.0);
t168 = t45+t158-t205;
t169 = pw1.*t2.*u1.*(1.0./2.0);
t170 = t3.*3.0;
t171 = t123.*8.0;
t172 = t4+t5+t6+t170+t171;
t173 = t169-ptau.*t21.*t126.*t172.*(1.0./4.0);
t174 = h.*ptau.*t21.*u1.*u3;
t175 = h.*ptau.*t21.*w1.*w3.*4.0;
t176 = t174+t175-pw3.*t119.*u1;
t177 = h.*ptau.*t21.*u1.*u2;
t178 = h.*ptau.*t21.*w1.*w2.*4.0;
t179 = t177+t178-pw2.*t119.*u1;
t180 = h.*ptau.*t21.*u1.*u4;
t181 = h.*ptau.*t21.*w1.*w4.*4.0;
t182 = t180+t181-pw4.*t119.*u1;
t183 = t120.^2;
t184 = (-h).^(5.0./2.0);
t206 = ptau.*t21.*t126.*t144.*(1.0./4.0);
t186 = t27+t136-t206;
t207 = ptau.*t21.*t126.*t135.*(1.0./4.0);
t187 = t33+t124-t207;
t208 = ptau.*t21.*t126.*t156.*(1.0./4.0);
t188 = t39+t147-t208;
t209 = ptau.*t21.*t126.*t167.*(1.0./4.0);
t189 = t45+t158-t209;
t190 = t169-ptau.*t21.*t119.*t121.*t172.*(1.0./4.0);
t191 = pw3.*t120.*t183.*u1;
t192 = ptau.*t21.*t184.*u1.*u3;
t193 = ptau.*t21.*t184.*w1.*w3.*4.0;
t194 = t191+t192+t193;
t195 = pw2.*t120.*t183.*u1;
t196 = ptau.*t21.*t184.*u1.*u2;
t197 = ptau.*t21.*t184.*w1.*w2.*4.0;
t198 = t195+t196+t197;
t199 = pw4.*t120.*t183.*u1;
t200 = ptau.*t21.*t184.*u1.*u4;
t201 = ptau.*t21.*t184.*w1.*w4.*4.0;
t202 = t199+t200+t201;
t210 = t186.*u3;
t211 = t187.*u4;
t212 = t188.*u1;
t213 = t189.*u2;
t214 = t210+t211+t212+t213;
t215 = conj(t119);
t216 = 1.0./t215;
t217 = t186.*u1;
t218 = t189.*u4;
t267 = t187.*u2;
t268 = t188.*u3;
t219 = t217+t218-t267-t268;
t227 = ptau.*t21.*t144.*t216.*(1.0./4.0);
t220 = t27+t136-t227;
t228 = ptau.*t21.*t135.*t216.*(1.0./4.0);
t221 = t33+t124-t228;
t229 = ptau.*t21.*t156.*t216.*(1.0./4.0);
t222 = t39+t147-t229;
t230 = ptau.*t21.*t167.*t216.*(1.0./4.0);
t223 = t45+t158-t230;
t224 = t186.*u2;
t225 = t187.*u1;
t251 = t188.*u4;
t252 = t189.*u3;
t226 = t224+t225-t251-t252;
t231 = t145.*u2;
t232 = t146.*u1;
t285 = t157.*u4;
t286 = t168.*u3;
t233 = t231+t232-t285-t286;
t234 = t145.*u3;
t235 = t146.*u4;
t236 = t157.*u1;
t237 = t168.*u2;
t238 = t234+t235+t236+t237;
t239 = pw2.*t2.*u2.*(1.0./2.0);
t240 = t4.*3.0;
t241 = t130.*8.0;
t242 = t3+t5+t6+t240+t241;
t243 = t239-ptau.*t21.*t126.*t242.*(1.0./4.0);
t244 = t177+t178-pw1.*t119.*u2;
t245 = h.*ptau.*t21.*u2.*u3;
t246 = h.*ptau.*t21.*w2.*w3.*4.0;
t247 = t245+t246-pw3.*t119.*u2;
t248 = h.*ptau.*t21.*u2.*u4;
t249 = h.*ptau.*t21.*w2.*w4.*4.0;
t250 = t248+t249-pw4.*t119.*u2;
t253 = t239-ptau.*t21.*t119.*t121.*t242.*(1.0./4.0);
t254 = pw1.*t120.*t183.*u2;
t255 = t196+t197+t254;
t256 = pw3.*t120.*t183.*u2;
t257 = ptau.*t21.*t184.*u2.*u3;
t258 = ptau.*t21.*t184.*w2.*w3.*4.0;
t259 = t256+t257+t258;
t260 = pw4.*t120.*t183.*u2;
t261 = ptau.*t21.*t184.*u2.*u4;
t262 = ptau.*t21.*t184.*w2.*w4.*4.0;
t263 = t260+t261+t262;
t264 = t145.*u1;
t265 = t168.*u4;
t306 = t146.*u2;
t307 = t157.*u3;
t266 = t264+t265-t306-t307;
t269 = t220.*u3;
t270 = t221.*u4;
t271 = t222.*u1;
t272 = t223.*u2;
t273 = t269+t270+t271+t272;
t274 = t214.*t273;
t275 = t220.*u1;
t276 = t223.*u4;
t308 = t221.*u2;
t309 = t222.*u3;
t277 = t275+t276-t308-t309;
t278 = t219.*t277;
t279 = t220.*u2;
t280 = t221.*u1;
t310 = t222.*u4;
t311 = t223.*u3;
t281 = t279+t280-t310-t311;
t282 = t226.*t281;
t283 = t274+t278+t282;
t284 = 1.0./t7.^2;
t287 = pw3.*t2.*u3.*(1.0./2.0);
t288 = t5.*3.0;
t289 = t151.*8.0;
t290 = t3+t4+t6+t288+t289;
t291 = t287-ptau.*t21.*t126.*t290.*(1.0./4.0);
t292 = t174+t175-pw1.*t119.*u3;
t293 = t245+t246-pw2.*t119.*u3;
t294 = h.*ptau.*t21.*u3.*u4;
t295 = h.*ptau.*t21.*w3.*w4.*4.0;
t296 = t294+t295-pw4.*t119.*u3;
t297 = t287-ptau.*t21.*t119.*t121.*t290.*(1.0./4.0);
t298 = pw1.*t120.*t183.*u3;
t299 = t192+t193+t298;
t300 = pw2.*t120.*t183.*u3;
t301 = t257+t258+t300;
t302 = pw4.*t120.*t183.*u3;
t303 = ptau.*t21.*t184.*u3.*u4;
t304 = ptau.*t21.*t184.*w3.*w4.*4.0;
t305 = t302+t303+t304;
t312 = pw4.*t2.*u4.*(1.0./2.0);
t313 = t6.*3.0;
t314 = t162.*8.0;
t315 = t3+t4+t5+t313+t314;
t316 = t312-ptau.*t21.*t126.*t315.*(1.0./4.0);
t317 = t180+t181-pw1.*t119.*u4;
t318 = t248+t249-pw2.*t119.*u4;
t319 = t294+t295-pw3.*t119.*u4;
t320 = t312-ptau.*t21.*t119.*t121.*t315.*(1.0./4.0);
t321 = pw1.*t120.*t183.*u4;
t322 = t200+t201+t321;
t323 = pw2.*t120.*t183.*u4;
t324 = t261+t262+t323;
t325 = pw3.*t120.*t183.*u4;
t326 = t303+t304+t325;
t327 = pw1.*t119;
t329 = h.*ptau.*t21.*u1.*2.0;
t328 = t327-t329;
t330 = pw1.*t2.*w1.*2.0;
t331 = u1.*w1.*1.6e1;
t332 = t14+t15+t16+t331;
t333 = t9+t10+t11-t19+t330-ptau.*t21.*t126.*t332.*(1.0./4.0);
t334 = pw1.*t120.*t183;
t335 = ptau.*t21.*t184.*u1.*2.0;
t336 = t334+t335;
t337 = t9+t10+t11-t19+t330-ptau.*t21.*t119.*t121.*t332.*(1.0./4.0);
t338 = pw2.*t119;
t340 = h.*ptau.*t21.*u2.*2.0;
t339 = t338-t340;
t341 = pw2.*t2.*w2.*2.0;
t342 = u2.*w2.*1.6e1;
t343 = t13+t15+t16+t342;
t344 = t8+t10+t11-t19+t341-ptau.*t21.*t126.*t343.*(1.0./4.0);
t345 = pw2.*t120.*t183;
t346 = ptau.*t21.*t184.*u2.*2.0;
t347 = t345+t346;
t348 = t8+t10+t11-t19+t341-ptau.*t21.*t119.*t121.*t343.*(1.0./4.0);
t349 = pw3.*t119;
t351 = h.*ptau.*t21.*u3.*2.0;
t350 = t349-t351;
t352 = pw3.*t2.*w3.*2.0;
t353 = u3.*w3.*1.6e1;
t354 = t13+t14+t16+t353;
t355 = t8+t9+t11-t19+t352-ptau.*t21.*t126.*t354.*(1.0./4.0);
t356 = pw3.*t120.*t183;
t357 = ptau.*t21.*t184.*u3.*2.0;
t358 = t356+t357;
t359 = t8+t9+t11-t19+t352-ptau.*t21.*t119.*t121.*t354.*(1.0./4.0);
t360 = pw4.*t119;
t362 = h.*ptau.*t21.*u4.*2.0;
t361 = t360-t362;
t363 = pw4.*t2.*w4.*2.0;
t364 = u4.*w4.*1.6e1;
t365 = t13+t14+t15+t364;
t366 = t8+t9+t10-t19+t363-ptau.*t21.*t126.*t365.*(1.0./4.0);
t367 = pw4.*t120.*t183;
t368 = ptau.*t21.*t184.*u4.*2.0;
t369 = t367+t368;
t370 = t8+t9+t10-t19+t363-ptau.*t21.*t119.*t121.*t365.*(1.0./4.0);
t371 = 1.0./h.^2;
t372 = conj(t184);
t373 = 1.0./t372;
t374 = t114+t115+t116+t117;
t375 = t371.*t374.*w1;
t376 = pw1.*t7.*t371.*(1.0./4.0);
t377 = ptau.*t21.*t144.*t373.*(3.0./8.0);
t378 = t375+t376+t377;
t379 = t371.*t374.*w2;
t380 = pw2.*t7.*t371.*(1.0./4.0);
t381 = ptau.*t21.*t135.*t373.*(3.0./8.0);
t382 = t379+t380+t381;
t383 = t371.*t374.*w3;
t384 = pw3.*t7.*t371.*(1.0./4.0);
t385 = ptau.*t21.*t156.*t373.*(3.0./8.0);
t386 = t383+t384+t385;
t387 = t371.*t374.*w4;
t388 = pw4.*t7.*t371.*(1.0./4.0);
t389 = ptau.*t21.*t167.*t373.*(3.0./8.0);
t390 = t387+t388+t389;
t391 = 1.0./h.^3;
t392 = 1.0./t90;
t397 = ptau.*t21.*t120.*t144.*t391.*t392.*(3.0./8.0);
t393 = t375+t376-t397;
t398 = ptau.*t21.*t120.*t135.*t391.*t392.*(3.0./8.0);
t394 = t379+t380-t398;
t399 = ptau.*t21.*t120.*t156.*t391.*t392.*(3.0./8.0);
t395 = t383+t384-t399;
t400 = ptau.*t21.*t120.*t167.*t391.*t392.*(3.0./8.0);
t396 = t387+t388-t400;
f = [w1;w2;w3;w4;u1.*(-1.0./4.0)+t2.*t82.*w1.*(1.0./2.0)-t21.*t54.*t55.*u1.*(1.0./4.0)-t21.*t55.*t60.*u3.*(1.0./4.0)-t21.*t55.*t63.*u2.*(1.0./4.0);u2.*(-1.0./4.0)+t2.*t82.*w2.*(1.0./2.0)+t21.*t54.*t55.*u2.*(1.0./4.0)-t21.*t55.*t60.*u4.*(1.0./4.0)-t21.*t55.*t63.*u1.*(1.0./4.0);u3.*(-1.0./4.0)+t2.*t82.*w3.*(1.0./2.0)+t21.*t54.*t55.*u3.*(1.0./4.0)-t21.*t55.*t60.*u1.*(1.0./4.0)+t21.*t55.*t63.*u4.*(1.0./4.0);u4.*(-1.0./4.0)+t2.*t82.*w4.*(1.0./2.0)-t21.*t54.*t55.*u4.*(1.0./4.0)-t21.*t55.*t60.*u2.*(1.0./4.0)+t21.*t55.*t63.*u3.*(1.0./4.0);-t21.*t47.*t68.*t90.*t103-t21.*t47.*t74.*t90.*t106-t21.*t47.*t80.*t90.*t109;-t89.*((t21.*t47.*t68.*t90.*t103+t21.*t47.*t74.*t90.*t106+t21.*t47.*t80.*t90.*t109).*(u1.*w1.*4.0+u2.*w2.*4.0+u3.*w3.*4.0+u4.*w4.*4.0)+t21.*t47.*t90.*t103.*(t7.*u1.*u3.*2.0+t7.*u2.*u4.*2.0)+t21.*t47.*t90.*t109.*(t7.*u1.*u2.*2.0-t7.*u3.*u4.*2.0)+t21.*t47.*t90.*t106.*(t3.*t7-t4.*t7-t5.*t7+t6.*t7)-1.0);pw1.*(1.0./4.0)-t21.*t22.*t47.*(t233.*(t33+t124-t207+t173.*u2-t125.*t176.*u4.*(1.0./2.0)+t125.*t179.*u1.*(1.0./2.0)-t125.*t182.*u3.*(1.0./2.0))+t238.*(t39+t147-t208+t173.*u3+t125.*t176.*u1.*(1.0./2.0)+t125.*t179.*u4.*(1.0./2.0)+t125.*t182.*u2.*(1.0./2.0))+t226.*(t33+t124-t185+t190.*u2-t2.*t121.*t194.*u4.*(1.0./2.0)+t2.*t121.*t198.*u1.*(1.0./2.0)-t2.*t121.*t202.*u3.*(1.0./2.0))+t214.*(t39+t147-t203+t190.*u3+t2.*t121.*t194.*u1.*(1.0./2.0)+t2.*t121.*t198.*u4.*(1.0./2.0)+t2.*t121.*t202.*u2.*(1.0./2.0))-t125.*t266.*(pw1.*t3.*t119.*3.0+pw1.*t4.*t119+pw1.*t5.*t119+pw1.*t6.*t119+pw1.*t119.*t123.*4.0+ph.*t184.*w1.*8.0-pw2.*t119.*u1.*u2.*2.0-pw3.*t119.*u1.*u3.*2.0+pw4.*t119.*u1.*u4.*2.0+pw2.*t119.*w1.*w2.*4.0+pw3.*t119.*w1.*w3.*4.0+pw4.*t119.*w1.*w4.*4.0-h.*ptau.*t3.*t21.*u1.*4.0-h.*ptau.*t6.*t21.*u1.*4.0-h.*ptau.*t21.*t123.*u1.*1.6e1-h.*ptau.*t21.*u4.*w1.*w4.*1.6e1).*(1.0./4.0)+t2.*t121.*t219.*(pw1.*t3.*t120.*t183.*3.0+pw1.*t4.*t120.*t183+pw1.*t5.*t120.*t183+pw1.*t6.*t120.*t183+pw1.*t120.*t123.*t183.*4.0-h.*ph.*t120.*t183.*w1.*8.0+ptau.*t3.*t21.*t184.*u1.*4.0+ptau.*t6.*t21.*t184.*u1.*4.0+ptau.*t21.*t123.*t184.*u1.*1.6e1-pw2.*t120.*t183.*u1.*u2.*2.0-pw3.*t120.*t183.*u1.*u3.*2.0+pw4.*t120.*t183.*u1.*u4.*2.0+pw2.*t120.*t183.*w1.*w2.*4.0+pw3.*t120.*t183.*w1.*w3.*4.0+pw4.*t120.*t183.*w1.*w4.*4.0+ptau.*t21.*t184.*u4.*w1.*w4.*1.6e1).*(1.0./4.0)).*(1.0./2.0)+t21.*t22.*t283.*t284.*u1;pw2.*(1.0./4.0)-t21.*t22.*t47.*(t233.*(t27+t136-t206+t243.*u1+t125.*t244.*u2.*(1.0./2.0)-t125.*t247.*u4.*(1.0./2.0)-t125.*t250.*u3.*(1.0./2.0))+t238.*(t45+t158-t209+t243.*u4+t125.*t244.*u3.*(1.0./2.0)+t125.*t247.*u1.*(1.0./2.0)+t125.*t250.*u2.*(1.0./2.0))+t226.*(t27+t136-t204+t253.*u1+t2.*t121.*t255.*u2.*(1.0./2.0)-t2.*t121.*t259.*u4.*(1.0./2.0)-t2.*t121.*t263.*u3.*(1.0./2.0))+t214.*(t45+t158-t205+t253.*u4+t2.*t121.*t255.*u3.*(1.0./2.0)+t2.*t121.*t259.*u1.*(1.0./2.0)+t2.*t121.*t263.*u2.*(1.0./2.0))+t125.*t266.*(pw2.*t3.*t119+pw2.*t4.*t119.*3.0+pw2.*t5.*t119+pw2.*t6.*t119+pw2.*t119.*t130.*4.0+ph.*t184.*w2.*8.0-pw1.*t119.*u1.*u2.*2.0+pw3.*t119.*u2.*u3.*2.0-pw4.*t119.*u2.*u4.*2.0+pw1.*t119.*w1.*w2.*4.0+pw3.*t119.*w2.*w3.*4.0+pw4.*t119.*w2.*w4.*4.0-h.*ptau.*t4.*t21.*u2.*4.0-h.*ptau.*t5.*t21.*u2.*4.0-h.*ptau.*t21.*t130.*u2.*1.6e1-h.*ptau.*t21.*u3.*w2.*w3.*1.6e1).*(1.0./4.0)-t2.*t121.*t219.*(pw2.*t3.*t120.*t183+pw2.*t4.*t120.*t183.*3.0+pw2.*t5.*t120.*t183+pw2.*t6.*t120.*t183+pw2.*t120.*t130.*t183.*4.0-h.*ph.*t120.*t183.*w2.*8.0+ptau.*t4.*t21.*t184.*u2.*4.0+ptau.*t5.*t21.*t184.*u2.*4.0+ptau.*t21.*t130.*t184.*u2.*1.6e1-pw1.*t120.*t183.*u1.*u2.*2.0+pw3.*t120.*t183.*u2.*u3.*2.0-pw4.*t120.*t183.*u2.*u4.*2.0+pw1.*t120.*t183.*w1.*w2.*4.0+pw3.*t120.*t183.*w2.*w3.*4.0+pw4.*t120.*t183.*w2.*w4.*4.0+ptau.*t21.*t184.*u3.*w2.*w3.*1.6e1).*(1.0./4.0)).*(1.0./2.0)+t21.*t22.*t283.*t284.*u2;pw3.*(1.0./4.0)-t21.*t22.*t47.*(t238.*(t27+t136-t206+t291.*u1+t125.*t292.*u3.*(1.0./2.0)+t125.*t293.*u4.*(1.0./2.0)+t125.*t296.*u2.*(1.0./2.0))-t233.*(t45+t158-t209+t291.*u4-t125.*t292.*u2.*(1.0./2.0)-t125.*t293.*u1.*(1.0./2.0)+t125.*t296.*u3.*(1.0./2.0))+t214.*(t27+t136-t204+t297.*u1+t2.*t121.*t299.*u3.*(1.0./2.0)+t2.*t121.*t301.*u4.*(1.0./2.0)+t2.*t121.*t305.*u2.*(1.0./2.0))-t226.*(t45+t158-t205+t297.*u4-t2.*t121.*t299.*u2.*(1.0./2.0)-t2.*t121.*t301.*u1.*(1.0./2.0)+t2.*t121.*t305.*u3.*(1.0./2.0))+t125.*t266.*(pw3.*t3.*t119+pw3.*t4.*t119+pw3.*t5.*t119.*3.0+pw3.*t6.*t119+pw3.*t119.*t151.*4.0+ph.*t184.*w3.*8.0-pw1.*t119.*u1.*u3.*2.0+pw2.*t119.*u2.*u3.*2.0-pw4.*t119.*u3.*u4.*2.0+pw1.*t119.*w1.*w3.*4.0+pw2.*t119.*w2.*w3.*4.0+pw4.*t119.*w3.*w4.*4.0-h.*ptau.*t4.*t21.*u3.*4.0-h.*ptau.*t5.*t21.*u3.*4.0-h.*ptau.*t21.*t151.*u3.*1.6e1-h.*ptau.*t21.*u2.*w2.*w3.*1.6e1).*(1.0./4.0)-t2.*t121.*t219.*(pw3.*t3.*t120.*t183+pw3.*t4.*t120.*t183+pw3.*t5.*t120.*t183.*3.0+pw3.*t6.*t120.*t183+pw3.*t120.*t151.*t183.*4.0-h.*ph.*t120.*t183.*w3.*8.0+ptau.*t4.*t21.*t184.*u3.*4.0+ptau.*t5.*t21.*t184.*u3.*4.0+ptau.*t21.*t151.*t184.*u3.*1.6e1-pw1.*t120.*t183.*u1.*u3.*2.0+pw2.*t120.*t183.*u2.*u3.*2.0-pw4.*t120.*t183.*u3.*u4.*2.0+pw1.*t120.*t183.*w1.*w3.*4.0+pw2.*t120.*t183.*w2.*w3.*4.0+pw4.*t120.*t183.*w3.*w4.*4.0+ptau.*t21.*t184.*u2.*w2.*w3.*1.6e1).*(1.0./4.0)).*(1.0./2.0)+t21.*t22.*t283.*t284.*u3;pw4.*(1.0./4.0)-t21.*t22.*t47.*(t238.*(t33+t124-t207+t316.*u2+t125.*t317.*u3.*(1.0./2.0)+t125.*t319.*u1.*(1.0./2.0)+t125.*t318.*u4.*(1.0./2.0))-t233.*(t39+t147-t208+t316.*u3-t125.*t317.*u2.*(1.0./2.0)-t125.*t318.*u1.*(1.0./2.0)+t125.*t319.*u4.*(1.0./2.0))+t214.*(t33+t124-t185+t320.*u2+t2.*t121.*t322.*u3.*(1.0./2.0)+t2.*t121.*t326.*u1.*(1.0./2.0)+t2.*t121.*t324.*u4.*(1.0./2.0))-t226.*(t39+t147-t203+t320.*u3-t2.*t121.*t322.*u2.*(1.0./2.0)-t2.*t121.*t324.*u1.*(1.0./2.0)+t2.*t121.*t326.*u4.*(1.0./2.0))-t125.*t266.*(pw4.*t3.*t119+pw4.*t4.*t119+pw4.*t5.*t119+pw4.*t6.*t119.*3.0+pw4.*t119.*t162.*4.0+ph.*t184.*w4.*8.0+pw1.*t119.*u1.*u4.*2.0-pw2.*t119.*u2.*u4.*2.0-pw3.*t119.*u3.*u4.*2.0+pw1.*t119.*w1.*w4.*4.0+pw2.*t119.*w2.*w4.*4.0+pw3.*t119.*w3.*w4.*4.0-h.*ptau.*t3.*t21.*u4.*4.0-h.*ptau.*t6.*t21.*u4.*4.0-h.*ptau.*t21.*t162.*u4.*1.6e1-h.*ptau.*t21.*u1.*w1.*w4.*1.6e1).*(1.0./4.0)+t2.*t121.*t219.*(pw4.*t3.*t120.*t183+pw4.*t4.*t120.*t183+pw4.*t5.*t120.*t183+pw4.*t6.*t120.*t183.*3.0+pw4.*t120.*t162.*t183.*4.0-h.*ph.*t120.*t183.*w4.*8.0+ptau.*t3.*t21.*t184.*u4.*4.0+ptau.*t6.*t21.*t184.*u4.*4.0+ptau.*t21.*t162.*t184.*u4.*1.6e1+pw1.*t120.*t183.*u1.*u4.*2.0-pw2.*t120.*t183.*u2.*u4.*2.0-pw3.*t120.*t183.*u3.*u4.*2.0+pw1.*t120.*t183.*w1.*w4.*4.0+pw2.*t120.*t183.*w2.*w4.*4.0+pw3.*t120.*t183.*w3.*w4.*4.0+ptau.*t21.*t184.*u1.*w1.*w4.*1.6e1).*(1.0./4.0)).*(1.0./2.0)+t21.*t22.*t283.*t284.*u4;-pu1-t21.*t22.*t47.*(t233.*(t333.*u2-t125.*t328.*u1.*w2+t125.*t328.*u3.*w4+t125.*t328.*u4.*w3)-t238.*(-t333.*u3+t125.*t328.*u1.*w3+t125.*t328.*u2.*w4+t125.*t328.*u4.*w2)+t214.*(t337.*u3+t2.*t121.*t336.*u1.*w3+t2.*t121.*t336.*u2.*w4+t2.*t121.*t336.*u4.*w2)+t226.*(t337.*u2+t2.*t121.*t336.*u1.*w2-t2.*t121.*t336.*u3.*w4-t2.*t121.*t336.*u4.*w3)-t125.*t266.*(ph.*t184.*u1.*2.0+pw1.*t119.*u1.*w1.*2.0-pw1.*t119.*u2.*w2+pw2.*t119.*u1.*w2-pw1.*t119.*u3.*w3+pw3.*t119.*u1.*w3+pw1.*t119.*u4.*w4+pw4.*t119.*u1.*w4-h.*ptau.*t3.*t21.*w1.*4.0-h.*ptau.*t21.*u1.*u4.*w4.*4.0)+t2.*t121.*t219.*(h.*ph.*t120.*t183.*u1.*-2.0+ptau.*t3.*t21.*t184.*w1.*4.0+pw1.*t120.*t183.*u1.*w1.*2.0-pw1.*t120.*t183.*u2.*w2+pw2.*t120.*t183.*u1.*w2-pw1.*t120.*t183.*u3.*w3+pw3.*t120.*t183.*u1.*w3+pw1.*t120.*t183.*u4.*w4+pw4.*t120.*t183.*u1.*w4+ptau.*t21.*t184.*u1.*u4.*w4.*4.0)).*(1.0./2.0);-pu2-t21.*t22.*t47.*(t233.*(t344.*u1-t125.*t339.*u2.*w1+t125.*t339.*u3.*w4+t125.*t339.*u4.*w3)-t238.*(-t344.*u4+t125.*t339.*u1.*w3+t125.*t339.*u3.*w1+t125.*t339.*u2.*w4)+t214.*(t348.*u4+t2.*t121.*t347.*u1.*w3+t2.*t121.*t347.*u3.*w1+t2.*t121.*t347.*u2.*w4)+t226.*(t348.*u1+t2.*t121.*t347.*u2.*w1-t2.*t121.*t347.*u3.*w4-t2.*t121.*t347.*u4.*w3)+t125.*t266.*(ph.*t184.*u2.*2.0+pw1.*t119.*u2.*w1-pw2.*t119.*u1.*w1+pw2.*t119.*u2.*w2.*2.0+pw2.*t119.*u3.*w3+pw3.*t119.*u2.*w3-pw2.*t119.*u4.*w4+pw4.*t119.*u2.*w4-h.*ptau.*t4.*t21.*w2.*4.0-h.*ptau.*t21.*u2.*u3.*w3.*4.0)-t2.*t121.*t219.*(h.*ph.*t120.*t183.*u2.*-2.0+ptau.*t4.*t21.*t184.*w2.*4.0+pw1.*t120.*t183.*u2.*w1-pw2.*t120.*t183.*u1.*w1+pw2.*t120.*t183.*u2.*w2.*2.0+pw2.*t120.*t183.*u3.*w3+pw3.*t120.*t183.*u2.*w3-pw2.*t120.*t183.*u4.*w4+pw4.*t120.*t183.*u2.*w4+ptau.*t21.*t184.*u2.*u3.*w3.*4.0)).*(1.0./2.0);-pu3+t21.*t22.*t47.*(t233.*(t355.*u4+t125.*t350.*u1.*w2+t125.*t350.*u2.*w1-t125.*t350.*u3.*w4)+t238.*(-t355.*u1+t125.*t350.*u3.*w1+t125.*t350.*u2.*w4+t125.*t350.*u4.*w2)-t214.*(t359.*u1+t2.*t121.*t358.*u3.*w1+t2.*t121.*t358.*u2.*w4+t2.*t121.*t358.*u4.*w2)+t226.*(t359.*u4-t2.*t121.*t358.*u1.*w2-t2.*t121.*t358.*u2.*w1+t2.*t121.*t358.*u3.*w4)-t125.*t266.*(ph.*t184.*u3.*2.0+pw1.*t119.*u3.*w1-pw3.*t119.*u1.*w1+pw2.*t119.*u3.*w2+pw3.*t119.*u2.*w2+pw3.*t119.*u3.*w3.*2.0-pw3.*t119.*u4.*w4+pw4.*t119.*u3.*w4-h.*ptau.*t5.*t21.*w3.*4.0-h.*ptau.*t21.*u2.*u3.*w2.*4.0)+t2.*t121.*t219.*(h.*ph.*t120.*t183.*u3.*-2.0+ptau.*t5.*t21.*t184.*w3.*4.0+pw1.*t120.*t183.*u3.*w1-pw3.*t120.*t183.*u1.*w1+pw2.*t120.*t183.*u3.*w2+pw3.*t120.*t183.*u2.*w2+pw3.*t120.*t183.*u3.*w3.*2.0-pw3.*t120.*t183.*u4.*w4+pw4.*t120.*t183.*u3.*w4+ptau.*t21.*t184.*u2.*u3.*w2.*4.0)).*(1.0./2.0);-pu4+t21.*t22.*t47.*(t233.*(t366.*u3+t125.*t361.*u1.*w2+t125.*t361.*u2.*w1-t125.*t361.*u4.*w3)+t238.*(-t366.*u2+t125.*t361.*u1.*w3+t125.*t361.*u3.*w1+t125.*t361.*u4.*w2)-t214.*(t370.*u2+t2.*t121.*t369.*u1.*w3+t2.*t121.*t369.*u3.*w1+t2.*t121.*t369.*u4.*w2)+t226.*(t370.*u3-t2.*t121.*t369.*u1.*w2-t2.*t121.*t369.*u2.*w1+t2.*t121.*t369.*u4.*w3)+t125.*t266.*(ph.*t184.*u4.*2.0+pw1.*t119.*u4.*w1+pw4.*t119.*u1.*w1+pw2.*t119.*u4.*w2-pw4.*t119.*u2.*w2+pw3.*t119.*u4.*w3-pw4.*t119.*u3.*w3+pw4.*t119.*u4.*w4.*2.0-h.*ptau.*t6.*t21.*w4.*4.0-h.*ptau.*t21.*u1.*u4.*w1.*4.0)-t2.*t121.*t219.*(h.*ph.*t120.*t183.*u4.*-2.0+ptau.*t6.*t21.*t184.*w4.*4.0+pw1.*t120.*t183.*u4.*w1+pw4.*t120.*t183.*u1.*w1+pw2.*t120.*t183.*u4.*w2-pw4.*t120.*t183.*u2.*w2+pw3.*t120.*t183.*u4.*w3-pw4.*t120.*t183.*u3.*w3+pw4.*t120.*t183.*u4.*w4.*2.0+ptau.*t21.*t184.*u1.*u4.*w1.*4.0)).*(1.0./2.0);ptau.*t21.*t373.*(-3.0./8.0)+t21.*t47.*t90.*(t214.*(t378.*u3+t382.*u4+t386.*u1+t390.*u2)+t219.*(t378.*u1-t382.*u2-t386.*u3+t390.*u4)+t226.*(t378.*u2+t382.*u1-t386.*u4-t390.*u3)+t273.*(t393.*u3+t395.*u1+t394.*u4+t396.*u2)+t277.*(t393.*u1-t394.*u2-t395.*u3+t396.*u4)+t281.*(t393.*u2+t394.*u1-t395.*u4-t396.*u3)).*(1.0./2.0)+t21.*t47.*t283.*t392.*(1.0./4.0);0.0];
