function g = get_target_g_u(in1)
%GET_TARGET_G_U
%    G = GET_TARGET_G_U(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    04-Jun-2021 21:25:58

u1 = in1(1,:);
u2 = in1(2,:);
u3 = in1(3,:);
u4 = in1(4,:);
g = [u1.^2-u2.^2-u3.^2+u4.^2;u1.*u2.*2.0-u3.*u4.*2.0;u1.*u3.*2.0+u2.*u4.*2.0];
