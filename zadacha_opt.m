function d_best = zadacha_opt(x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%скорость истечения
%global u_dv;


th_mars=x(1);

%максимальный расход
gamma = x(2);

if gamma>4e-03 && gamma<9.1e-03
    [v_opt,t_opt] = zad_exp(gamma);
    d_best = find_zad_opt(th_mars,gamma,v_opt,t_opt);
else
    d_best=1e+15;
end
d_best
end

