function dgduv = get_dgduv(in1,in2)
%GET_DGDUV
%    DGDUV = GET_DGDUV(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    13-Apr-2021 17:33:27

u1 = in1(1,:);
u2 = in1(2,:);
u3 = in1(3,:);
u4 = in1(4,:);
v1 = in2(1,:);
v2 = in2(2,:);
v3 = in2(3,:);
v4 = in2(4,:);
dgduv = reshape([u1.*2.0,u2.*2.0,u3.*2.0,v1,v2,v3,u2.*-2.0,u1.*2.0,u4.*2.0,-v2,v1,v4,u3.*-2.0,u4.*-2.0,u1.*2.0,-v3,-v4,v1,u4.*2.0,u3.*-2.0,u2.*2.0,v4,-v3,v2,0.0,0.0,0.0,u1,u2,u3,0.0,0.0,0.0,-u2,u1,u4,0.0,0.0,0.0,-u3,-u4,u1,0.0,0.0,0.0,u4,-u3,u2],[6,8]);