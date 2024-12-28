xmesh = linspace(0,3,5);
solinit = bvpinit(xmesh,@guess);
sol = bvp4c(@bvpfun, @bcfun, solinit);
plot(sol.x(1,:),sol.y(1,:),'-o')
function dydx = bvpfun(x,y) % equation being solved
dydx = [y(2)
        y(1)];
end
%-------------------------------------------
function res = bcfun(ya,yb) % boundary conditions
res = [ya(1)
       yb(1)-1];
end
%-------------------------------------------
function y = guess(x) % guess at solution behavior
y = [exp(x)
     exp(x)];
end