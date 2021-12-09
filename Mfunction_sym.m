w = sym('w', [1 3],'real')';
u = sym('u', [1 3],'real')';
%M = Mfunction(u);
M = [[u(1) u(3) -u(2)];
    [u(2) -u(1) u(3)];
    [u(3) u(2) -u(1)]];

% M = [[-u(2) u(1) u(3)];
%     [u(3) u(2) -u(1)];
%     [u(1) u(3) -u(2)]];

%M=diag(u)
F1=M*u;
%chi = norm(F1);
%chi = u'*u;
chi=u'*u+4*w'*w;
dF1du=jacobian(F1, u);

k = sym('k','real')';
mu = sym('mu','real')';
F2=(dF1du*w)/chi;
dF2du=jacobian(F2, u);
dF2dw=jacobian(F2, w);
rightpart = (dF2du*w-k*dF2dw*u)/chi;
rightpart= simplify(rightpart);
leftpart=-mu*F1/(norm(F1)^3);

simplify(dF2du*w)

% function M = Mfunction(u)
% %UNTITLED2 Summary of this function goes here
% %   Detailed explanation goes here
% M = [[u(1) u(3) -u(2)];
%     [u(2) -u(1) u(3)];
%     [u(3) u(2) -u(1)]];
% end