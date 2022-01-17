function theta = traj2TrueAnomaly(rr)
%принимает на вход массив Nx2
angles = atan2(rr(:,2),rr(:,1));
angles(angles<0)=2*pi+angles(angles<0);
diffs = abs(angles(2:end)-angles(1:end-1));
n=sum(diffs>pi);
thetaRad=n*2*pi+angles(end)-angles(1);
theta=thetaRad/(2*pi);
end