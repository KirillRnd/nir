%внешние константы
ae = 149597870700;
mug_0 = 132712.43994*(10^6)*(10^(3*3));
T_earth_days = 365.256363004;
T_earth = T_earth_days*3600*24;
T_mars=T_earth*1.8808476;
T_mars_days = T_earth_days*1.8808476;

%поворот кватернионом в плоскость эклиптики
n_eps = [1; 0; 0];
eps_0 = 2*pi*(23+26/60+21/3600)/360;
q_eps = [cos(eps_0/2); n_eps*sin(eps_0/2)];
%коэффициенты обезразмеривания
r_unit=ae;
V_unit=sqrt(mug_0/ae);
T_unit = T_earth/(2*pi);
%гиперпараметры
t_start = 0;
orbits='Flat';
N=1350;
m0=367;
eta=0.45;
planet_start = 'Earth';
planet_end = 'Mars';



%формулы аппроксимации
pr_x = @(a)(0.20074211987250948)/(7.780825627926492*a-0.34159993346898454+sin(2*pi*a+0.1231174475945867+0.5853586584957071/(a^2-1.5690074915306857*a+1.9819223815154108)) + (-41.26470185582318+86.35986725318767*a-51.08371071757919*a^2)*exp(-3.424188046763629*a) ); %OK
pr_y = @(a)(0.0015463527047873246 - 0.0031678016585249486 * cos(2*pi * a+(-2.7056930690628103*a^2+5.176300671261158*a-2.336378597152512)*exp(3.2003639491958467*a-3.2205106566099353*a^2)))/(a^2-0.18328217249773662*a+0.11514749689845788);%OK
pv_x = @(a)(0.0018661818762522732 - 0.0031537940448213170 * cos(2*pi * a+(-7.277086750486235*a^2+13.74021831427378*a-6.262552219285274)*exp(1.381195073394106*a-2.278956114228774*a^2)))/(a^2-0.19152227381595716*a+0.0993237695309777);
pv_y = @(a)(0.10456877009148134)/(4.050922375438426*a-0.15902042687671913+sin(2*pi*a+0.11429217276196703+0.012218372217935064/(a^2-2.8347986183365705*a+2.06664561825753)) +(-540.5981059117084*a^2+958.0032106019531*a-415.15730979365355)*exp(-5.877538308252977*a));

J_a  = @(a)(0.002117*exp(a))/(a*exp(a)+0.17856*(sin(7.44443*a)-a));
%% перебираем разные дальности
N_count_3 = 169;
N_count_4 = 101;
AN_MARS_lin = linspace(0.75,6, N_count_3);
deltaAN_MARS_lin = linspace(-1,1, N_count_4);
AN_MARS_predicted = zeros(N_count_3, N_count_4);
AN_MARS_real = zeros(N_count_3, N_count_4);
savefilename = 'ThesisMarsCheck_naive_3.mat';
%load(savefilename)
for i = 1:N_count_3 
    disp([num2str(i), '/', num2str(N_count_3)])
    for j = 1:N_count_4 
        if AN_MARS_real(i,j)>0
            continue
        end
        %интегрируем
        %определяем стартовое положение
        
        t0 = 0;
        AN_i = AN_MARS_lin(i);
        OM_i = 2*pi*(0.27941185*AN_i-0.01523959);
        dAN_i = deltaAN_MARS_lin(j);
        AN_MARS_predicted(i,j) = AN_i;
        if AN_i+dAN_i<0.75
            continue
        end
        B=0.2721831;
        a_rel=1.52;
        T_i = 2*pi*(AN_i*(a_rel-a_rel*B*log(a_rel))+a_rel*B^2*log(a_rel));
        %T_i = 2*pi*(AN_i*1.35036912+0.02855872);
        pr_0 = [pr_x(AN_i+dAN_i),pr_y(AN_i+dAN_i),0];
        pv_0 = [pv_x(AN_i+dAN_i),pv_y(AN_i+dAN_i),0];
        tspan = linspace(t0,t0+T_i, 1000);

        st.t = 0;
        st.planet = planet_start;
        st.mode = orbits;
        st.delta_omega = OM_i;
        
        [start_pos, start_vel] = planetModel(st);
        start_pos=start_pos*1e+03/ae;
        start_vel=start_vel*1e+03/V_unit;
        y0 = cat(2,start_pos,start_vel,pr_0,pv_0)';
        
        options = odeset('AbsTol',1e-10);
        options = odeset(options,'RelTol',1e-10);   
        [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
        AN_MARS_real(i,j) = calculate_angular_distance(y(:,1:3))/(2*pi);
        
        save(savefilename);
    end
end
%%
figure(1)
AN_MARS_real_fix = AN_MARS_real-AN_MARS_predicted;
AN_MARS_corrected = AN_MARS_predicted+repmat(deltaAN_MARS_lin,N_count_3,1);
AN_MARS_real_fix(AN_MARS_corrected<0.75)=nan;
[X, Y] = meshgrid(AN_MARS_lin, deltaAN_MARS_lin);
s = surf(X,AN_MARS_real_fix', Y, 'FaceColor','blue','FaceAlpha',0.4);
s.EdgeColor = 'none';
hold on;
s = contour3(X,AN_MARS_real_fix', Y,'ShowText','on', 'Color', 'k');
hold off;
Y_forapprox = reshape(Y.',1,[]);
A_forapprox = reshape(AN_MARS_real_fix'.',1,[]);
p = polyfit(A_forapprox(~isnan(A_forapprox))',Y_forapprox(~isnan(A_forapprox))',1)
grid;
xlabel('Угловая дальность')
ylabel('Начальная ошибка угловой дальности')
zlabel('Невязка на правом конце')
%set(gca, 'ZScale', 'log')
%title('Реальная угловая дальность')
%colorbar;
grid;
%set(gca,'ColorScale','log')
view(0,90)
%%
function res = fsolve_traj(z,y0,yf,t0,t_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tspan = linspace(t0,t0+t_end, 10);
y0(7:12)=z;
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);   
options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 0.4, 5));
[~,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
res = y(end,1:6)-yf(1:6);
%res = [norm(y(end,1:3)-yf(1:3));norm(y(end,4:6)-yf(4:6))];
%res = [norm(y(end,1:2)-yf(1:2));y(end,3)-yf(3); norm(y(end,4:5)-yf(4:5));y(end,6)-yf(6);];
%res = sign(res).*log(sign(res).*res+1);
end

function [c, ceq] = fmincon_traj(z,y0,yf,t0,t_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tspan = linspace(t0,t0+t_end, 10);
y0(7:12)=z;
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);   
[~,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
c =[];
ceq = (y(end,1:6)-yf(1:6));
end
function [value, isterminal, direction] = eventIntegrationTrajStopR(~, y, r_end, r_outer)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


r=y(1:3);
value = [0,0];
value(1) = norm(r) - r_end;
if value(1) < 0
    value(1)=0;
end
value(2) = r_outer-norm(r);
if value(2) < 0
    value(2)=0;
end
isterminal = [1 1];
direction = [0 0];
end
function R = calculateRotMatrix(r0,v0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
e1 = r0/norm(r0);
e2 = v0/norm(v0);
e3 = cross(e1,e2);
R = [e1;e2;e3];
end

function y_g = guess(x,y0,z) % guess at solution behavior
tspan = linspace(0,x, 2);
y0(7:12)=z;
if x>0
    options = odeset('AbsTol',1e-10);
    options = odeset(options,'RelTol',1e-10);   
    options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopR(s, y, 0.4, 5));
    [~,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
    y_g = y(end,:);
else
    y_g = y0;
end
end

function res = bcfcn(ya, yb, y0, yf) % boundary conditions
res = [ya(1:6) - y0(1:6);
       yb(1:6) - yf(1:6)];
end