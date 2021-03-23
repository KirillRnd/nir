function theta = getTrueAnomaly(r,V,mug)
%UNTITLED Summary of this function goes here
c=cross(r,V);
n=c/norm(c);
f=cross(V,c)-mug*r/norm(r);
if norm(f) == 0
    f = [1 0 0];
end

theta_cos=f'*r/(norm(f)*norm(r));
theta_sin=norm(cross(f, r))/norm(f)/norm(r);
theta = atan2(theta_sin,theta_cos);
if cross(f, r)'*n < 0
    theta=2*pi-theta;
end
end

