%внешние константы
ae = 149597870700;
mug_0 = 132712.43994*(10^6)*(10^(3*3));
T_earth = 365.256363004*3600*24;
T_mars=T_earth*1.8808476;
T_mars_days = 365.256363004*1.8808476;

r_unit=ae;
V_unit=sqrt(mug_0/ae);
T_unit = T_earth/(2*pi);
%гиперпараметры
t_start=0;
orbits='Flat';
N=1350;
m0=367;
eta=0.45;
planet_start = 'Earth';
planet_end = 'Mars';
omega = -pi;

%определяем стартовое положение
st.t = t_start;
st.planet = planet_start;
st.mode = orbits;
st.delta_omega = omega;

[start_pos, start_vel] = planetModel(st);
start_pos=start_pos*1e+03/ae;
start_vel=start_vel*1e+03/V_unit;

% pr_x = @(a)(0.20074211987250948)/(7.780825627926492*a-0.34159993346898454+sin(2*pi*a+0.1231174475945867+0.5853586584957071/(a^2-1.5690074915306857*a+1.9819223815154108)) + (-41.26470185582318+86.35986725318767*a-51.08371071757919*a^2)*exp(-3.424188046763629*a) ); %OK
% pr_y = @(a)(0.0015463527047873246 - 0.0031678016585249486 * cos(2*pi * a+(-2.7056930690628103*a^2+5.176300671261158*a-2.336378597152512)*exp(3.2003639491958467*a-3.2205106566099353*a^2)))/(a^2-0.18328217249773662*a+0.11514749689845788);%OK
% pv_x = @(a)(0.0018661818762522732 - 0.0031537940448213170 * cos(2*pi * a+(-7.277086750486235*a^2+13.74021831427378*a-6.262552219285274)*exp(1.381195073394106*a-2.278956114228774*a^2)))/(a^2-0.19152227381595716*a+0.0993237695309777);
% pv_y = @(a)(0.10456877009148134)/(4.050922375438426*a-0.15902042687671913+sin(2*pi*a+0.11429217276196703+0.012218372217935064/(a^2-2.8347986183365705*a+2.06664561825753)-0.0019877228343783706*a^2) +(-540.5981059117084*a^2+958.0032106019531*a-415.15730979365355)*exp(-5.877538308252977*a));


pr_x = @(a,d)(0.20074211987250948)/(7.780825627926492*a-0.34159993346898454+sin(2*pi*a+0.1231174475945867+d*0.5853586584957071/(a^2-1.5690074915306857*a+1.9819223815154108)) + d*(-41.26470185582318+86.35986725318767*a-51.08371071757919*a^2)*exp(-3.424188046763629*a) ); %OK
pr_y = @(a,d)(0.0015463527047873246 - 0.0031678016585249486 * cos(2*pi * a+d*(-2.7056930690628103*a^2+5.176300671261158*a-2.336378597152512)*exp(3.2003639491958467*a-3.2205106566099353*a^2)))/(a^2-0.18328217249773662*a+0.11514749689845788);%OK
pv_x = @(a,d)(0.0018661818762522732 - 0.0031537940448213170 * cos(2*pi * a+d*(-7.277086750486235*a^2+13.74021831427378*a-6.262552219285274)*exp(1.381195073394106*a-2.278956114228774*a^2)))/(a^2-0.19152227381595716*a+0.0993237695309777);
pv_y = @(a,d)(0.10456877009148134)/(4.050922375438426*a-0.15902042687671913+sin(2*pi*a+0.11429217276196703+d*0.012218372217935064/(a^2-2.8347986183365705*a+2.06664561825753)) +d*(-540.5981059117084*a^2+958.0032106019531*a-415.15730979365355)*exp(-5.877538308252977*a));

B=0.2721831;
A=0.5433279;
T_a  = @(a, a_rel)(a_rel-a_rel*B*log(a_rel))*a+(a_rel*B^2*log(a_rel));
a_OM = @(OM, a_rel,k)(OM+k+A*exp(-a_rel)*log(a_rel)*log(a_rel))/(A*(1+exp(-a_rel))*log(a_rel));
OM_a = @(a, a_rel)A*(a*(exp(a_rel)+1)-log(a_rel))*exp(-a_rel)*log(a_rel);

%% перебираем разные дальности
step = 1/32;
ds = 3/4:step:12/2; %угловая дальность

step2 = 0.01;
da = 0.4:step2:2.5; %соотношение полуосей
L2 = length(ds);
L3 = length(da);

k_true = zeros(L2,L3);
k_preds = zeros(L2,L3);
for i = 1:L2 
    disp([num2str(i), '/', num2str(L2)])
    for j = 1:L3 
        t0 = t_start;
        AN_i = ds(i);        
        a_rel = da(j);

        T_i = 2*pi*T_a(AN_i, a_rel);
        K = calculateKMatrix_ver2(a_rel,AN_i,1);
        d_coef = 1;
        [pr_0,pv_0] = get_initial_adjoint(AN_i,a_rel,d_coef);
        pr_0 = pr_0';
        pv_0 = pv_0';
        y0 = cat(2,start_pos,start_vel,pr_0,pv_0)';

        tspan = linspace(t0,t0+T_i, AN_i*400);
        options = odeset('AbsTol',1e-12);
        options = odeset(options,'RelTol',1e-12);   
        %options = odeset(options, 'Events',@(s, y) eventIntegrationTrajStopE0(s, y));
        [t,y] = ode113(@(t,y) internalIntegration3D(t,y), tspan,y0,options);
        [a,eMag,i_2,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(y(end, 1:3)',y(end, 4:6)',1);

        k_true(i,j) = a_rel;
        k_preds(i,j) = a;
    end
end


%%
figure(1)
[X, Y] = meshgrid(ds, da);
k_preds_fix = k_preds';
%k_preds_fix(abs(k_preds_fix)>4) = nan;
k_diff = k_preds_fix-Y;
s = surf(X,Y,k_diff-1, 'FaceColor','blue','FaceAlpha',0.4);
s.EdgeColor = 'none';
hold on;
k_diff_contour3 = [-0.1,-0.01,0,0.1,0.3];
s = contour3(X,Y, k_diff,k_diff_contour3,'ShowText','on', 'Color', 'k');
hold off;
%Y_forapprox = reshape(Y.',1,[]);
%A_forapprox = reshape(AN_MARS_real_fix'.',1,[]);
%p = polyfit(A_forapprox(~isnan(A_forapprox))',Y_forapprox(~isnan(A_forapprox))',1)
grid;
xlabel('Угловая дальность')
ylabel('Соотношение радиусов')
zlabel('Соотношение радиусов 2')
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
options = odeset('AbsTol',1e-12);
options = odeset(options,'RelTol',1e-12);   
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
function [value, isterminal, direction] = eventIntegrationTrajStopR(~, y, r_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

r=y(1:3);

value = norm(r) - r_end;
isterminal = 1;
direction = 0;
end
function [value, isterminal, direction] = eventIntegrationTrajStopE0(~, y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

r=y(1:3);
v=y(4:6);
e = getEccentricity(r,v,1);
value = e;
isterminal = 1;
direction = 0;
end
function R = calculateRotMatrix(r0,v0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
e1 = r0;
e2 = v0;
e3 = cross(e1,e2);
R = [e1;e2;e3];
end
function [Kr,Kv] = calculateKMatrix(k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
a1 = (k+2.087217245873957)*log(k)/k;
a2 = -k*exp(-8.563932009313163/k)+69.18759116365857*exp(-8.563932009313163*k);
a3 = -0.07470468492737094*log(k)+(0.6602914816446814*log(k)*exp(-k))/(k^3/2);
a4 = (-0.8348043664990472*k+2.922882836448715*log(k)+3.6786032392739227)*log(k)/k;
a2=0;
a3=0;
a4=a1;
Kr = [[a1,a2,0];[a3,a4,0];[0,0,1]];
Kv = [[a4,a3,0];[a2,a1,0];[0,0,1]];
end
function K = calculateKMatrix_ver2(k_target,a,d)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%коррекция kappa (k)
c0 = 1.568311392349845;
c1 = (1.52-1)/(exp(-c0/1.52)-exp(-c0));
c2 = 1 - c1*exp(-c0);
k = c2+c1*exp(-c0/k_target);
%вычисляем gamma
g = (-8.234450871689381-1.835384888679592*cos(2*pi*a+4.568045799135197)+12.140923790903319*a)/...
    (6.521910293172844*a-4.627546590315015+sin(2*pi*a+2.963335050753793+d*(0.010167241207198103)/(a^2-1.4200382540034495*a+0.5258798683545889))...
    + d*(-22.253725602121175*a^2+47.5730432773046*a-24.334701239930087)*exp(-3.0200712580477367*a) );
f = -g/1.52+1/log(1.52);
gamma = (f*k+g)*log(k)/k;

K = gamma*eye(3);
end
function [pr_0,pv_0] = get_initial_adjoint(a,k,d)
%Adjoint variables approximation
% a - angular distance, revolutions
% k - radii ratio
% d - multiplier of D coefficients, 0 or 1
%MARS
pr_x_M = @(a,d)(0.20074211987250948)/(7.780825627926492*a-0.34159993346898454...
    +sin(2*pi*a+0.1231174475945867+d*0.5853586584957071/(a^2-1.5690074915306857*a+1.9819223815154108))...
    + d*(-41.26470185582318+86.35986725318767*a-51.08371071757919*a^2)*exp(-3.424188046763629*a) ); %OK
pr_y_M = @(a,d)(0.0015463527047873246 - 0.0031678016585249486 * cos(2*pi * a+...
    d*(-2.7056930690628103*a^2+5.176300671261158*a-2.336378597152512)*exp(3.2003639491958467*a-3.2205106566099353*a^2)))/(a^2-0.18328217249773662*a+0.11514749689845788);%OK
pv_x_M = @(a,d)(0.0018661818762522732 - 0.0031537940448213170 * cos(2*pi * a+...
    d*(-7.277086750486235*a^2+13.74021831427378*a-6.262552219285274)*exp(1.381195073394106*a-2.278956114228774*a^2)))/(a^2-0.19152227381595716*a+0.0993237695309777);
pv_y_M = @(a,d)(0.10456877009148134)/(4.050922375438426*a-0.15902042687671913...
    +sin(2*pi*a+0.11429217276196703+d*0.012218372217935064/(a^2-2.8347986183365705*a+2.06664561825753)) ...
    +d*(-540.5981059117084*a^2+958.0032106019531*a-415.15730979365355)*exp(-5.877538308252977*a));
%VENUS
pr_x_V = @(a,d)(-0.2489269721475946)/(7.737402704800786*a-0.32876870077839637...
    +sin(2*pi*a+0.10377485561964216+d*0.5010768929874032/(a^2-1.9514741381243845*a+2.4016982504424953))...
    + d*(-39.740029616989275+73.28781124258543*a-40.57760476908843*a^2)*exp(-3.086139681061194*a) ); %OK
pr_y_V = @(a,d)(-0.006779491667620585 + 0.00432410667098884 * cos(2*pi * a...
    +d*(-5.520942046733101*a^2+9.93390376988929*a-4.541315615665989)*exp(1.3808955605425999*a-2.0000002462936357*a^2)))/(a^2+0.08726450366262754*a-0.20030536675890054);%OK
pv_x_V = @(a,d)(-0.006319471111610397 + 0.004278503661785389 * cos(2*pi * a...
    +d*(-1.786880982902473*a^2+3.1517137043764167*a-1.3891212954291232)*exp(2.9860241498160076*a-2.699241528126778*a^2)))/(a^2+0.13118547856852233*a-0.24254361523409607);
pv_y_V = @(a,d)(-0.12867786064634879)/(4.0002764433454105*a-0.15629991248717665...
    +sin(2*pi*a+0.08483203399561368+d*0.008896394730198013/(a^2-2.9131929806240318*a+2.1648040548553382)) +...
    d*(-1354.631280683438*a^2+2264.6420824634483*a-951.3449842396692)*exp(-6.814266658923464*a));

k1=0.72;
k2=1.52;

%prx
p1 = pr_x_V(a,d);
p2 = pr_x_M(a,d);
f = (p1/log(k1)-(k2/k1)*p2/log(k2))/(1-k2/k1);
g = p2*k2/log(k2)-f*k2;
pr_x = @(k)(log(k)*(f*k+g)/k);
%pry
p1 = pr_y_V(a,d);
p2 = pr_y_M(a,d);
g = (p2*k2^(3/2)/log(k2) - p1*k1*sqrt(k2)/log(k1))/(1-sqrt(k2/k1));
f = p1*k1/log(k1)-g/sqrt(k1);
pr_y = @(k)(log(k)*(f/k+g/(k^(3/2))));
%pvx
p1 = pv_x_V(a,d);
p2 = pv_x_M(a,d);
g = (p2*k2^(3/2)/log(k2) - p1*k1*sqrt(k2)/log(k1))/(1-sqrt(k2/k1));
f = p1*k1/log(k1)-g/sqrt(k1);
pv_x = @(k)(log(k)*(f/k+g/(k^(3/2))));
%pvy
p1 = pv_y_V(a,d);
p2 = pv_y_M(a,d);
f = (p1/log(k1)-(k2/k1)*p2/log(k2))/(1-k2/k1);
g = p2*k2/log(k2)-f*k2;
pv_y = @(k)(log(k)*(f*k+g)/k);

pr_0 = [pr_x(k);pr_y(k);0.];
pv_0 = [pv_x(k);pv_y(k);0.]; 
end