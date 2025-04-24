function [r, V] = planetModelSimple(t,planet)
%t в днях
ae = 149597870.700; %км
mug = 132712.43994*(10^6); %км3с−2
if strcmp(planet, 'Earth')
    %https://ssd.jpl.nasa.gov/horizons/app.html#/
    a = 1.00*ae;
    e = 1.682781268763978E-02;
    Omega=pi*2.198210990481491E+02/180;
    omega=pi*2.419337242483841E+02/180;
    i=pi*(2.312212390408654E-03)/180;
    M0=4.600206682616619E-02;
    t0=juliandate(2016,01,03);
elseif strcmp(planet, 'Mars')
    a = 1.523783440561007*ae;
    e = 9.336536157600901E-02;
    Omega=pi*4.950841351899370E+01/180;
    omega=pi*2.866264313647392E+02/180;
    i=pi*(1.848382482872302E+00)/180;
    M0=0;
    t0=2457690.998223170638;
elseif strcmp(planet, 'Venus')
    a = 0.723327331758308*ae;
    e = 6.751702316628624E-03;
    Omega=pi*7.663475158642598E+01/180;
    omega=pi*5.508545549669039E+01/180;
    i=pi*(3.394389948608993E+00)/180;
    M0=0;
    t0=2457355.851127521135;
elseif strcmp(planet, 'Asteroid2015RE36')
    % https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=2015RE36

    a = 1.309904433959459*ae;
    e = 0.3903331967297961;
    Omega=pi*182.222933497807/180;
    omega=pi*50.71383675614203/180;
    i=pi*0.08154670244768725/180;
    M0=0;
    t0=2460446.638519069714;
elseif strcmp(planet, 'Polyhymnia')
    % https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=33

    a = 2.873381502364615*ae;
    e = 0.331735924639039;
    Omega=pi*8.42872062375473/180;
    omega=pi*338.6500899677241/180;
    i=pi*(1.851970677506626)/180;
    M0=0;
    t0=2460470.619383189565;
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