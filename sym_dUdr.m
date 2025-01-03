function sym_dUdr = sym_dUdr(in1)
%sym_dUdr
%    sym_dUdr = sym_dUdr(IN1)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    10-Dec-2024 10:01:59

r1 = in1(1,:);
r2 = in1(2,:);
r3 = in1(3,:);
t2 = r1.^2;
t3 = r2.^2;
t4 = r3.^2;
t5 = t2+t3+t4;
t6 = 1.0./t5.^(3.0./2.0);
sym_dUdr = [-r1.*t6;-r2.*t6;-r3.*t6];
end
