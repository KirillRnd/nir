%Выводим траекторию в параметрических переменных"
fig1 = figure(5);
plot3(0, 0, 0, 'k--o')
set(gca,'FontSize',14)
hold on;

th = linspace(0 ,4*pi,1000)';

mars_traj_ks = arrayfun(@(r1, r2, r3) rToU([r1,r2,r3], phi), mars_traj_New(:, 1),mars_traj_New(:, 2),mars_traj_New(:, 3),'UniformOutput',false);
mars_traj_ks = cell2mat(mars_traj_ks')';

mars_traj_ks_zero = arrayfun(@(r1, r2, r3) rToU([r1,r2,r3], phi0), mars_traj_New(:, 1),mars_traj_New(:, 2),mars_traj_New(:, 3),'UniformOutput',false);
mars_traj_ks_zero = cell2mat(mars_traj_ks_zero')';

%phi === 0
earth_traj_ks = arrayfun(@(r1, r2, r3) rToU([r1,r2,r3], phi0), earth_traj_New(:, 1),earth_traj_New(:, 2),earth_traj_New(:, 3),'UniformOutput',false);
earth_traj_ks = cell2mat(earth_traj_ks')';
plot3(earth_traj_ks(:, 1), earth_traj_ks(:, 2), earth_traj_ks(:, 3), 'k')
plot3(mars_traj_ks(:, 1), mars_traj_ks(:, 2), mars_traj_ks(:, 3), 'r')
%plot3(-mars_traj_ks_zero(:, 1), -mars_traj_ks_zero(:, 2), -mars_traj_ks_zero(:, 3), 'r--')
%plot3(mars_traj_ks_zero(:, 1), mars_traj_ks_zero(:, 2), mars_traj_ks_zero(:, 3), 'r--')
plot3(uu(:, 1), uu(:, 2), uu(:, 3), 'b', 'LineWidth', 2.5);
%a_scale=3e-01/mean(vecnorm(a_ks, 2, 2));
a_scale=0;
d = 24*3600;
idxes=1;
for i=1:ceil(t(end)/d)
    ix = find(t>d*i*10, 1);
    idxes=[idxes, ix];
end    
for i = idxes
    %plot3([uu(i, 1), uu(i, 1)+a_scale*a_ks(i, 1)], [uu(i, 2), uu(i, 2)+a_scale*a_ks(i, 2)],[uu(i, 3), uu(i, 3)+a_scale*a_ks(i, 3)],'k')
end

%plot3(uu(end, 1), uu(end, 2), uu(end, 3),'bO')
%plot3(mars_r_f(1), mars_r_f(2),mars_r_f(3),'rO')
axis equal

%title('Траектория КА KS')
xlabel('u1')
ylabel('u2')
zlabel('u3')
grid on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
hold off;
set(fig1,'color','w');

f = getframe(fig1);
[im,map] = rgb2ind(f.cdata,128,'nodither');
im(1,1,1,128) = 0;

for i=1:128
    an=i*2*pi/128;
    x=cos(an);
    y=sin(an);
    v = [x y 0.6];
    view(v)
    pause(1/24)
    %F(i) = getframe(fig1);
    f = getframe(fig1);
    drawnow
    im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
    %saveas(fig1,'gif/frame'+string(i)+'.png')
end
imwrite(im,map,'KSRotation.gif','DelayTime',1/24,'LoopCount',inf)
% writerObj = VideoWriter('KSRotation.avi');
% writerObj.FrameRate = 24;
%   % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);