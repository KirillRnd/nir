%% оптимальные поверхности
Mevery_fix = Jevery;
Mevery_fix(Mevery_fix==0)=nan;
figure(11);
M1 = Mevery_fix(1:L1, 1:L2,1)';
X_axis = 180*omega_space(1:L1)/pi;
Y_axis = ds(1:L2);
%Mscaled1 = 100*M1/m0;
Mscaled1 = M1;
s = surf(180*omega_space(1:L1)/pi,ds(1:L2),Mscaled1,'DisplayName','Семейство 1', 'FaceColor','cyan');
s.EdgeColor = 'none';
hold on;

Tevery_fix = Tevery;
Tevery_fix(Tevery_fix==0)=nan;
T1 = Tevery_fix(1:L1, 1:L2,1)';
T2 = Tevery_fix(1:L1, 1:L2,2)';
T3 = Tevery_fix(1:L1, 1:L2,3)';
T4 = Tevery_fix(1:L1, 1:L2,4)';
T5 = Tevery_fix(1:L1, 1:L2,5)';
figure(19);
s_contour = contour3(180*omega_space(1:L1)/pi,ds(1:L2),T1,(1:10)*365);
figure(11)
hold on;
K_contours = 1;
while K_contours<length(s_contour)
    H_contours = s_contour(1,K_contours);
    N_contours = s_contour(2,K_contours);
    X_contours = s_contour(1,K_contours+1:K_contours+N_contours);
    Y_contours = s_contour(2,K_contours+1:K_contours+N_contours);
    Z_contours = interp2(X_axis,Y_axis,Mscaled1, X_contours, Y_contours);
    K_contours = K_contours+N_contours+1;
    plot3(X_contours, Y_contours,Z_contours,'black', 'HandleVisibility','off')
    %disp(K_contours)
end
%C2 = rescale(M1,256,511);
M2 = Mevery_fix(1:L1, 1:L2,2)';
%Mscaled2 = 100*M2/m0;

Mscaled2 = M2;
s = surf(180*omega_space(1:L1)/pi,ds(1:L2),Mscaled2,'DisplayName','Семейство 2', 'FaceColor','red');
s.EdgeColor = 'none';
figure(19);
hold on;
s_contour = contour3(180*omega_space(1:L1)/pi,ds(1:L2),T2,(1:10)*365);
figure(11)
K_contours = 1;
while K_contours<length(s_contour)
    H_contours = s_contour(1,K_contours);
    N_contours = s_contour(2,K_contours);
    X_contours = s_contour(1,K_contours+1:K_contours+N_contours);
    Y_contours = s_contour(2,K_contours+1:K_contours+N_contours);
    Z_contours = interp2(X_axis,Y_axis,Mscaled2, X_contours, Y_contours);
    K_contours = K_contours+N_contours+1;
    plot3(X_contours, Y_contours,Z_contours,'black','HandleVisibility','off')
    %disp(K_contours)
end

M3 = Mevery_fix(1:L1, 1:L2,3)';
%Mscaled3 = 100*M3/m0;

Mscaled3 = M3;
s = surf(180*omega_space(1:L1)/pi,ds(1:L2),Mscaled3,'DisplayName','Семейство 3', 'FaceColor','magenta');
s.EdgeColor = 'none';

figure(19);
s_contour = contour3(180*omega_space(1:L1)/pi,ds(1:L2),T3,(1:10)*365);
figure(11)
K_contours = 1;
while K_contours<length(s_contour)
    H_contours = s_contour(1,K_contours);
    N_contours = s_contour(2,K_contours);
    X_contours = s_contour(1,K_contours+1:K_contours+N_contours);
    Y_contours = s_contour(2,K_contours+1:K_contours+N_contours);
    Z_contours = interp2(X_axis,Y_axis,Mscaled3, X_contours, Y_contours);
    K_contours = K_contours+N_contours+1;
    plot3(X_contours, Y_contours,Z_contours,'black','HandleVisibility','off')
    %disp(K_contours)
end

M4 = Mevery_fix(1:L1, 1:L2,4)';
%Mscaled4 = 100*M4/m0;

Mscaled4 = M4;
s = surf(180*omega_space(1:L1)/pi,ds(1:L2),Mscaled4,'DisplayName','Семейство 4', 'FaceColor','green');
s.EdgeColor = 'none';

figure(19);
s_contour = contour3(180*omega_space(1:L1)/pi,ds(1:L2),T4,(1:10)*365);
hold off;
figure(11)
K_contours = 1;
while K_contours<length(s_contour)
    H_contours = s_contour(1,K_contours);
    N_contours = s_contour(2,K_contours);
    X_contours = s_contour(1,K_contours+1:K_contours+N_contours);
    Y_contours = s_contour(2,K_contours+1:K_contours+N_contours);
    Z_contours = interp2(X_axis,Y_axis,Mscaled4, X_contours, Y_contours);
    K_contours = K_contours+N_contours+1;
    plot3(X_contours, Y_contours,Z_contours,'black','HandleVisibility','off')
    %disp(K_contours)
end

M5 = Mevery_fix(1:L1, 1:L2,5)';
Mscaled5 = M5;
s = surf(180*omega_space(1:L1)/pi,ds(1:L2),Mscaled5,'DisplayName','Семейство 0', 'FaceColor','green');
s.EdgeColor = 'none';

figure(19);
s_contour = contour3(180*omega_space(1:L1)/pi,ds(1:L2),T5,(1:10)*365);
hold off;
figure(11)
K_contours = 1;
while K_contours<length(s_contour)
    H_contours = s_contour(1,K_contours);
    N_contours = s_contour(2,K_contours);
    X_contours = s_contour(1,K_contours+1:K_contours+N_contours);
    Y_contours = s_contour(2,K_contours+1:K_contours+N_contours);
    Z_contours = interp2(X_axis,Y_axis,Mscaled5, X_contours, Y_contours);
    K_contours = K_contours+N_contours+1;
    plot3(X_contours, Y_contours,Z_contours,'black','HandleVisibility','off')
    %disp(K_contours)
end

%colormap(s,spring)
xlabel('Разность фаз, градусов')
ylabel('Угловая дальность, витков')
%zlabel('Затраты топлива, %')
zlabel('J, функционал')

set(gca,'zscale','log')
hold off;
xlim([-180, 180])
legend;
dcm = datacursormode;
dcm.Enable = 'on';


mega_st.PXevery = PXevery;
mega_st.PHIevery = PHIevery;
mega_st.Tevery_fix = Tevery_fix;
mega_st.ds = ds;
mega_st.omega_space = omega_space;
mega_st.rad = rad;
mega_st.UorR = UorR;
mega_st.decreaseNonPsysical = decreaseNonPsysical;
mega_st.modifier_f = modifier_f;
mega_st.eta = eta;
mega_st.case_traj = case_traj;
mega_st.planet_end = planet_end;
mega_st.display = display;
mega_st.terminal_state = terminal_state;
mega_st.orbits = orbits;
mega_st.Mevery_fix = Mevery_fix;

dcm.UpdateFcn = @(obj,event_obj)datatipWithSubscript(obj,event_obj, mega_st);

figure(19);
close;
%%
%находим точки минимума
ind_optimum = zeros(4,L2);%(i_ind, F, m)
for j = 1:L2
    m_local = Mevery_fix(1:L1, j,:);
    m_min = min(m_local,[],'all');
    [i_min, F_min] = find(m_local==m_min);
    
    
    ind_optimum(:,j) = [i_min, j, m_min, F_min];
end


figure(11);
hold on;
for F = 1:5
    F_local = ind_optimum(:,ind_optimum(4,:)==F);
    plot3(180*omega_space(F_local(1,:))/pi,ds(F_local(2,:)),F_local(3,:),'magenta', 'LineWidth', 2)
end
hold off;
figure(21);
hold on;
for F = 1:5
    F_local = ind_optimum(:,ind_optimum(4,:)==F);
    PV_local = squeeze(PVevery_fix(:,:,1:2,F));
    i_min = F_local(1,:);
    j_min = F_local(2,:);
    plot3(diag(PV_local(i_min,j_min,1)), diag(PV_local(i_min,j_min,2)), ds(j_min), 'black', 'LineWidth', 2)
end
hold off;
%%
function output_txt = datatipWithSubscript(obj,event_obj, mega_st)
% Display the position of the data cursor with index/subscript information
%
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
%
% A simple datatip text update function that additionally displays the
% linear index or the subscripts of an object's selected point.
pos = get(event_obj,'Position');
output_txt = {['X: ',num2str(pos(1),4)],...
              ['Y: ',num2str(pos(2),4)]};
% If there is a Z-coordinate in the position, display it as well. Also get
% the size of the object.
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
    dataSize  = size(get(get(event_obj,'Target'),'ZData'));
else
    dataSize  = size(get(get(event_obj,'Target'),'YData'));
end
% Get the linear index
dataIndex = get(event_obj,'DataIndex');
if sum(dataSize>1)>1
    % display subscripts for data with more than 1 dimension
   
    S = cell(1,numel(dataSize));
    %global i_one j_one;
    %[i_one, j_one] = cell2mat(S);
    [S{:}] = ind2sub(dataSize,dataIndex);
    i_one = S{2};
    j_one = S{1};
    Z_slice = mega_st.Mevery_fix(i_one,j_one,:);
    F_one = find(Z_slice==pos(3));
    subsString = sprintf('%i,',cell2mat(S));

    plotTrajectoryCallback(i_one, j_one, F_one, mega_st)
    output_txt{end+1} = ['@ (', subsString(1:end-1), ')'];
else
    % display linear index for vector data
    output_txt{end+1} = ['@ ',num2str(dataIndex)];
end
end