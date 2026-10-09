clear; clc; close all;

%% Lunar constants
mu = 4902.8;           % km^3/s^2
R_moon = 1737.4;       % km

%% Initial lunar orbit
h_orbit = 100;         % km
r1 = R_moon + h_orbit; % km

%% Landing site: Site07
site_lat = -88.811008; % deg
site_lon = 123.690068; % deg
h_site = 0.80055;      % km
r2 = R_moon + h_site;  % km

%% Transfer orbit
a = (r1 + r2)/2;       % km

v_circ = sqrt(mu/r1);  % km/s
v_apo = sqrt(mu*(2/r1 - 1/a));
v_peri = sqrt(mu*(2/r2 - 1/a));

%% Delta-V
dv_deorbit = abs(v_circ - v_apo);

% Idealized impulsive braking at periapsis
dv_braking = v_peri;

dv_total = dv_deorbit + dv_braking;

%% Transfer time
t_transfer = pi*sqrt(a^3/mu); % s

%% Results
fprintf('LUNAR DESCENT BASELINE\n');
fprintf('Site: Peak Near Shackleton\n\n');

fprintf('Circular orbit velocity: %.4f km/s\n',v_circ);
fprintf('Transfer apoapsis velocity: %.4f km/s\n',v_apo);
fprintf('Transfer periapsis velocity: %.4f km/s\n',v_peri);

fprintf('\nDeorbit Delta-V: %.2f m/s\n',dv_deorbit*1000);
fprintf('Ideal braking Delta-V: %.2f m/s\n',dv_braking*1000);
fprintf('Total ideal Delta-V: %.2f m/s\n',dv_total*1000);

fprintf('\nTransfer time: %.2f minutes\n',t_transfer/60);
