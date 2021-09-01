mug=1;
u = sym('u', [1 4],'real')';
w = sym('w', [1 4],'real')';
h=-mug/(u'*u+4*w'*w);

%pu = sym('pu', [1 4],'real')';
%pw = sym('pw', [1 4],'real')';
L = L_KS(u);
g=[L*u;2*sqrt(-2*h)*L*w/(u'*u)];
g=[g(1:3);g(5:7)];

dgduw=jacobian(g,[u;w]);

pr = sym('pr', [1 3],'real')';
pv = sym('pv', [1 3],'real')';



P = [pr;pv]'*dgduw;
pu=simplify(P(1:4)');
pw=simplify(P(5:8)');