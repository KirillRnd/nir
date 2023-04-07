%скрипт для построения множества критических решений

mug=1;
r0 = [1; 0];
v0 = [0; 1];
T=2*pi;
r = sym('r', [2 1],'real');
assume(norm(r),'positive');
sym_dUdr=-mug*r/norm(r)^3;

sym_ddUdrdr=jacobian(sym_dUdr, r);
pV = sym('pV', [2 1],'real');

ddUdrdrpV=sym_ddUdrdr*pV;
jac_sym_ddUdrdr=jacobian(ddUdrdrpV',r);

dUdr = matlabFunction(sym_dUdr, 'Vars', {r});
ddUdrdr = matlabFunction(sym_ddUdrdr, 'Vars', {r});

z0=zeros([4,1]);
y0 = cat(1, r0,v0,z0);

T_end=T; %время перелёта
options = odeset('AbsTol',1e-10);
options = odeset(options,'RelTol',1e-10);  
tspan=linspace(0, T_end,100);
[~,y_solution] = ode113(@(t,y) internalIntegration2D(t,y),tspan,y0,options);

r_solution = y_solution(:,1:2);
%v_solution = y_solution(:,3:4);
%простая проверка работы интегрирования
% figure(1)
% plot(r_solution(:,1),r_solution(:,2))
% axis equal;


    

%лучше нечётное, чтобы ноль попал в сетку
L_grid=5;

C4D=zeros([L_grid,L_grid,L_grid,L_grid]);
d=1;
PV1=linspace(-d,d,L_grid);
PV2=linspace(-d,d,L_grid);
DPV1=linspace(-d,d,L_grid);
DPV2=linspace(-d,d,L_grid);
%parpool(3)
for i=1:L_grid
    for j=1:L_grid
        for k=1:L_grid
            z00=[PV1(i); PV2(j);DPV1(k); 0];
            parfor v=1:L_grid
                z0=z00;
                z0(4)=DPV2(v);
                C = calculateCondition(y0, z0,T_end);
                C4D(i,j,k,v)=C;
            end
        end
    end
end
%%
zero_coor = round(L_grid/2);
figure(2)
surf(PV1, PV2, C4D(:,:,zero_coor,zero_coor));
%axis equal;
set(gca,'ZScale','log')
xlabel('pv(1)')
ylabel('pv(2)')

figure(3)
surf(DPV1, DPV2, reshape(C4D(zero_coor,zero_coor,:,:), [L_grid,L_grid]));
%axis equal;
set(gca,'ZScale','log')
xlabel('dpvdt(1)')
ylabel('dpvdt(2)')

function C = calculateCondition(y0, z0,T_end)
    y0(5:8) = z0;
    
    options = odeset('AbsTol',1e-10);
    options = odeset(options,'RelTol',1e-10);  
    tspan=linspace(0, T_end,2);
    
    
    %вычисление матрицы чувствительности
    ddeltady0=zeros(4,4);
    step_h=1e-6;
    for i=5:8
        y0_delta=zeros([8,1]);
        y0_delta(i)=step_h;
        [~,y] = ode113(@(t,y) internalIntegration2D(t,y),tspan,y0+y0_delta, options);
        p_plus=y(end,1:4);
        
        [~,y] = ode113(@(t,y) internalIntegration2D(t,y),tspan,y0-y0_delta, options);
        p_minus=y(end,1:4);
        partial=(p_plus-p_minus)/(2*step_h);
        ddeltady0(:,i-4)=partial;
    end
    C=cond(ddeltady0);
end