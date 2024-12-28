load('D:\MATLAB\ipm\mat-files\07-Aug-2024.mat')
display = 0;

for i=7:7%L2 %угловая дальность
    for j=1:3%L4 %величина скорости
        [skip, i_nearest, j_nearest] = checkNear(skipGrid_Vhyp_end, i, j);
        if skip==0
            px_new=reshape(PXevery_Vhyp_end(i_nearest,j_nearest,:),[1,8]);
            if skip > 1 && norm(px_new) == 0
                continue
            end
            phi_new=PHIevery_Vhyp_end(i_nearest,j_nearest);
            omega=OMevery_Vhyp_end(i_nearest,j_nearest); %оптимальная разность фаз
            hi=HIevery_Vhyp_end(i_nearest,j_nearest);    %оптимальное направление
            disp([i,j])
            px=px_new;
            phi=phi_new;
            delta_s=ds(i)*2*pi;
            dV_value_end=dV_range(j);
            %dV_add = dV_value*[cos(dV_rad), sin(dV_rad), 0];
            modifier_p=10^(-4-sqrt(delta_s));
            x0_sec = [px delta_s/(2*pi) phi/(2*pi)];
            integration_acc=1e-12;
            calculate_condition=1;

            checkMethod_params = struct();
            checkMethod_params.t_start = t_start;
            checkMethod_params.delta_s = delta_s;
            checkMethod_params.rad = rad;
            checkMethod_params.UorR = UorR;
            checkMethod_params.decreaseNonPsysical = decreaseNonPsysical;
            checkMethod_params.modifier_p = modifier_p;
            checkMethod_params.modifier_f = modifier_f;
            checkMethod_params.x0_sec = x0_sec;
            checkMethod_params.eta = eta;
            checkMethod_params.case_traj = case_traj;
            checkMethod_params.planet_end = planet_end;
            checkMethod_params.display = display;
            checkMethod_params.terminal_state = terminal_state;
            checkMethod_params.integration_acc = integration_acc;
            checkMethod_params.calculate_condition = calculate_condition;
            checkMethod_params.orbits = orbits;
            checkMethod_params.omega = omega;
            checkMethod_params.hi = hi;
            checkMethod_params.a_rel = 1.52;
            %checkMethod_params.dV_add = dV_add;
            %checkMethod_params.dV_value = dV_value;
            checkMethod_params.dV_value_end = dV_value_end;
            
            minimize_delta_omega = @(delta_omega) find_close_solution(delta_omega, checkMethod_params, omega);
    

            try
                angle_range = pi/32;
                [delta_omega, Jt_end_min] = fminbnd(minimize_delta_omega, -angle_range, angle_range);
                omega_min = delta_omega+omega;
                checkMethod_params.omega = omega_min;
                checkMethod_results = checkMethod_ver2(checkMethod_params);
                %[dr, dV, C, px, s_f, phi, t_end, s, uu, rr, VV, t, Jt, a_ks, evaluation_time] 
            catch
                skipGrid_Vhyp_end(i,j)=2;
                continue
            end
        
            if checkMethod_results.dr < 10000 && length(checkMethod_results.t)==1000
                m=massLP(checkMethod_results.Jt, m0, N);
                Mevery_Vhyp_end(i,j)=m(1)-m(end);
                CONDevery_Vhyp_end(i,j)=checkMethod_results.C;
                PXevery_Vhyp_end(i,j,:)=checkMethod_results.px;
                PHIevery_Vhyp_end(i,j)=checkMethod_results.phi;
                T_end = checkMethod_results.t(end)/(24*60*60); %в днях
                Tevery_Vhyp_end(i,j)=T_end;
                Jevery_Vhyp_end(i,j)=checkMethod_results.Jt(end)*eta;
                ANevery_Vhyp_end(i,j) = calculate_angular_distance(checkMethod_results.rr(:,1:3));
    
                OMevery_Vhyp_end(i,j)=omega_min;
                HIevery_Vhyp_end(i,j)=checkMethod_results.hi_opt;
                skipGrid_Vhyp_end(i,j)=0; %set status only in the end
                savefilename = join(['mat-files/',date]);
                save(savefilename);
                disp('OK')
            else
                skipGrid_Vhyp_end(i,j)=2;
                continue
            end
            
        end
    end
end

function [skip, i_nearest, j_nearest] = checkNear(arr, i, j)
    % checkNear
    s = size(arr);
    if arr(i,j)==0
        %уже есть решение
        skip = 0;
        i_nearest = i;
        j_nearest = j;
    else
        i_j_range = [i,i-1,i,i+1,i-1,i-1,j+1,j+1;
                    j-1,j,j+1,j,j-1,j+1,j+1,j-1];
        idx_rearrange = [2, 4, 1, 3, 5, 6, 7, 8];
        %idx_rearrange = randperm(8);
        i_j_range = i_j_range(:,idx_rearrange);
        for k=1:8
            i_nearest=i_j_range(1,k);
            j_nearest=i_j_range(2,k);
            if i_nearest>=1 && i_nearest <= s(1) && j_nearest>=1 && j_nearest <= s(2)
                if arr(i_nearest, j_nearest) == 0
                    skip=2;
                    return
                end
            end
        end
        %нет соседних решений
        skip = -1;
        i_nearest = 0;
        j_nearest = 0;
    end
end

function Jt_end = find_close_solution(delta_omega, checkMethod_params, omega)
    omega_min = delta_omega+omega;
    checkMethod_params_COPY =checkMethod_params;
    checkMethod_params_COPY.omega = omega_min;
    checkMethod_results = checkMethod_ver2(checkMethod_params_COPY);
    Jt_end = checkMethod_results.Jt(end)*checkMethod_params_COPY.eta;
    if checkMethod_results.dr > 10000
        Jt_end=Jt_end*10;
    end
end
