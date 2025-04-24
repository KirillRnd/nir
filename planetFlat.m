function [r, V] = planetFlat(t,planet, delta_omega)
%t в днях
ae = 149597870.700; %км
mug = 132712.43994*(10^6); %км3с−2
if strcmp(planet, 'Earth')
    a = 1*ae;
    e = 0;
    Omega=0;
    omega=0;
    i=0;
    M0=0;
    t0=0;
elseif strcmp(planet, 'Mars')
    a = 1.52*ae;
    e = 0;
    Omega=0;
    omega=delta_omega;
    i=0;
    M0=0;
    t0=0;
elseif strcmp(planet, 'Venus')
    a = 0.72*ae;
    e = 0;
    Omega=0;
    omega=delta_omega;
    i=0;
    M0=0;
    t0=0;
elseif strcmp(planet, 'Asteroid2015RE36')
%     Эксцентриситет (e):  0,39;
%     Большая полуось (a): 1,31 а.е.;
%     Наклонение к плоскости эклиптики (i): 0,08°;
%     Долгота восходящего узла: 182°;
%     Аргумент перигелия: 51°;
%     Период обращения: 547 сут (1,5 года);
%     Дата прохождения перигелия: 16.11.2016
%     Абсолютная светимость (H): 21,9 mag;
%     Приблизительные размеры астероида: 130 – 300 м

    a = 1.31*ae;
    e = 0.39;
    Omega=pi*182/180;
    omega=pi*51/180;
    i=pi*(0.08-23.44)/180;
    M0=0;
    t0=juliandate(2016,11,16);
elseif strcmp(planet, 'Polyhymnia')
    % Эксцентриситет (e) 	0,3374924
    % Большая полуось (a) 	428,796 млн км
    % (2,8663230 а. е.)
    % Перигелий (q) 	284,08 млн км
    % (1,8989608 а. е.)
    % Афелий (Q) 	573,511 млн км
    % (3,8336852 а. е.)
    % Период обращения (P) 	1772,499 сут (4,853 г.)
    % Средняя орбитальная скорость 	17,079 км/с
    % Наклонение (i) 	1,87038°
    % Долгота восходящего узла (Ω) 	8,59380°
    % Аргумент перигелия (ω) 	338,14805°
    % Средняя аномалия (M) 	175,24936°
    % https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=33

    a = 2.866*ae;
    e = 0.337;
    Omega=pi*8.59/180;
    omega=pi*338.1/180;
    i=pi*(1.87-23.44)/180;
    M0=0;
    t0=juliandate(2024,06,09);
end
p=a*(1-e^2);
T=2*pi*sqrt(a^3/mug);%секунды
n=2*pi/(T/3600/24);  %рад/день
L_t = length(t);
if L_t == 1
    M = M0+(t-t0)*n;
    E = EAnomaly(e, M);
    nu=2*atan(sqrt((1+e)/(1-e))*tan(E/2));
    
    [r, V] = orb2rv(p,e,i,Omega,omega,nu,mug);
elseif L_t == 2
    t2 = t(1)+t(2);
    M = M0+(t2-t0)*n;
    E = EAnomaly(e, M);
    nu=2*atan(sqrt((1+e)/(1-e))*tan(E/2));
    
    [r, V] = orb2rv(p,e,i,Omega,omega,nu,mug);
elseif L_t > 2
    M = M0+(t-t0)*n;
    E = EAnomaly(e, M);
    nu=2*atan(sqrt((1+e)/(1-e))*tan(E/2));
    [r, V] = orb2rv(repmat(p, 1, L_t),repmat(e, 1, L_t),repmat(i, 1, L_t),repmat(Omega, 1, L_t),repmat(omega, 1, L_t),nu,mug);
end
r=r';
V=V';
end