mug=1;
u = sym('u', [1 4],'real')';
w = sym('w', [1 4],'real')';
%h = sym('h','real');
h=-mug/(u'*u+4*(w'*w));
tau = sym('tau','real');

Vp = sym('Vp', [1 3],'real')';
Vhyp = sym('Vhyp', [1 3],'real')';
%билинейное соотношение
bil=@(x)[x(4); -x(3); x(2); -x(1)];

u2 = u'*u;
L = L_KS(u);
Lt=L(1:3,:);

%фиксированный гип. избыток
delta_V = 2*sqrt(-2*h)*Lt*w/(u'*u)-Vp;
jac_delta_V = jacobian(delta_V,[u;w]);
jac_delta_V = simplify(jac_delta_V);

dVhyp_duw = delta_V'*jac_delta_V;
dVhyp_duw = simplify(dVhyp_duw);
F_r = Lt*u;
F_r_duw=jacobian(F_r,[u;w]);
F_duw = [F_r_duw; dVhyp_duw];
ortdgduw_Vhyp=null(F_duw);

matlabFunction(F_duw,'File','get_dgduv_Vhyp','Optimize', true, 'Vars', {u,w,Vp});
matlabFunction(ortdgduw_Vhyp,'File','get_ortdgduv_Vhyp','Optimize', true, 'Vars', {u,w,Vp});