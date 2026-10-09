clear; clc; close all;

%% Lunar constants
mu = 4902.8;           % km^3/s^2
R = 1737.4;            % km
g0 = 9.80665;          % m/s^2

%% Starship HLS assumptions
m0 = 250000;           % kg
m_dry = 100000;        % kg
Isp = 350;             % s
Tmax = 2e6;            % N

%% Landing site
h_site = 0.80055;      % km
r_site = R + h_site;

%% Transfer orbit
r1 = R + 100;          % km
a = (r1 + r_site)/2;
v_apo = sqrt(mu*(2/r1 - 1/a));

s0 = [r1; 0; 0; v_apo];

%% Coast to 15 km altitude
opts = odeset('RelTol',1e-10,'AbsTol',1e-11, ...
    'Events',@(t,s) altitudeEvent(t,s,R,15));

[t_coast,s_coast,te_coast] = ode45( ...
    @(t,s) coastODE(t,s,mu),[0 4000],s0,opts);

if isempty(te_coast)
    error('Coast did not reach the 15 km altitude.');
end

%% Initial powered descent state
r0 = s_coast(end,1:2)'*1000;  % m
v0 = s_coast(end,3:4)'*1000;  % m/s

% State: x, y, vx, vy, mass, accumulated Delta-V
y0 = [r0; v0; m0; 0];

%% Parameters
p.mu = mu*1e9;         % m^3/s^2
p.R = R*1000;          % m
p.r_site = r_site*1000;
p.Isp = Isp;
p.g0 = g0;
p.Tmax = Tmax;
p.m_dry = m_dry;

%% Powered descent
opts2 = odeset('RelTol',1e-8,'AbsTol',1e-8, ...
    'Events',@(t,y) touchdownEvent(t,y,p));

[t,y,te_touch] = ode45(@(t,y) descentODE(t,y,p), ...
    [0 2000],y0,opts2);

%% Extract results
r = vecnorm(y(:,1:2),2,2);
v = vecnorm(y(:,3:4),2,2);

alt = r-p.r_site;
mass = y(:,5);
dv = y(:,6);

prop_used = m0-mass(end);

%% Calculate thrust history
thrust = zeros(length(t),1);

for i = 1:length(t)
    [~,thrust(i)] = descentODE(t(i),y(i,:)',p);
end

%% Results
fprintf('STARSHIP HLS POWERED DESCENT V2\n\n');

fprintf('Initial altitude: %.2f km\n',alt(1)/1000);
fprintf('Final altitude: %.2f m\n',alt(end));
fprintf('Final speed: %.2f m/s\n',v(end));
fprintf('Descent time: %.2f s\n',t(end));

fprintf('\nPowered Delta-V: %.2f m/s\n',dv(end));
fprintf('Deorbit Delta-V: %.2f m/s\n', ...
    (sqrt(mu/r1)-v_apo)*1000);

fprintf('Total Delta-V: %.2f m/s\n', ...
    dv(end)+(sqrt(mu/r1)-v_apo)*1000);

fprintf('\nPropellant used: %.2f kg\n',prop_used);
fprintf('Remaining mass: %.2f kg\n',mass(end));

if ~isempty(te_touch) && v(end)<=2
    fprintf('\nSOFT LANDING ACHIEVED\n');
elseif ~isempty(te_touch)
    fprintf('\nTOUCHDOWN TOO FAST\n');
else
    fprintf('\nDESCENT DID NOT REACH SURFACE\n');
end

%% Plots
figure;
tiledlayout(3,2);

nexttile;
plot(t,alt/1000,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Altitude (km)');
title('Altitude');
grid on;

nexttile;
plot(t,v,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Speed (m/s)');
title('Velocity');
grid on;

nexttile;
plot(t,mass,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Mass (kg)');
title('Spacecraft Mass');
grid on;

nexttile;
plot(t,dv,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Delta-V (m/s)');
title('Accumulated Powered Delta-V');
grid on;

nexttile;
plot(t,thrust/1000,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Thrust (kN)');
title('Engine Thrust');
grid on;

nexttile;
plot(t,m0-mass,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Propellant Used (kg)');
title('Propellant Consumption');
grid on;

%% Functions
function ds = coastODE(~,s,mu)
    r = s(1:2);
    v = s(3:4);

    ds = [v; -mu*r/norm(r)^3];
end

function [value,isterminal,direction] = ...
    altitudeEvent(~,s,R,h)

    value = norm(s(1:2))-R-h;
    isterminal = 1;
    direction = -1;
end

function [dy,Tmag] = descentODE(~,y,p)

    r = y(1:2);
    v = y(3:4);
    m = y(5);

    rmag = norm(r);
    rhat = r/rmag;

    vr = dot(v,rhat);
    vt = v-vr*rhat;

    h = max(rmag-p.r_site,0);

    %% Desired vertical velocity
    v_des = -min(40,sqrt(2*0.8*h));

    %% Variable controller gains
    if h > 1000
        Kp = 0.15;
    elseif h > 100
        Kp = 0.4;
    else
        Kp = 1.0;
    end

    %% Guidance acceleration
    a_rad = Kp*(v_des-vr);
    a_tan = -0.03*vt;

    g = -p.mu*r/rmag^3;

    a_cmd = a_rad*rhat+a_tan-g;

    %% Thrust limits
    Tvec = m*a_cmd;
    Tmag = min(norm(Tvec),p.Tmax);

    if m <= p.m_dry
        Tmag = 0;
    end

    if norm(Tvec)>0
        u = Tvec/norm(Tvec);
    else
        u = [0;0];
    end

    thrust = Tmag*u;

    %% Equations of motion
    dy = [v;
          g+thrust/m;
          -Tmag/(p.Isp*p.g0);
          Tmag/m];
end

function [value,isterminal,direction] = ...
    touchdownEvent(~,y,p)

    value = norm(y(1:2))-p.r_site;
    isterminal = 1;
    direction = -1;
end


%% Additional validation

r_final = y(end,1:2)';
v_final = y(end,3:4)';

rhat = r_final/norm(r_final);

v_radial = dot(v_final,rhat);
v_tangential = norm(v_final-v_radial*rhat);

dv_powered = dv(end);
dv_deorbit = (sqrt(mu/r1)-v_apo)*1000;

% Rocket equation consistency check
dv_rocket = Isp*g0*log(m0/mass(end));

% Propellant required for initial deorbit burn
m_before_deorbit = m0*exp(dv_deorbit/(Isp*g0));
prop_deorbit = m_before_deorbit-m0;

fprintf('\nLANDING VALIDATION\n');
fprintf('Radial velocity: %.3f m/s\n',v_radial);
fprintf('Tangential velocity: %.3f m/s\n',v_tangential);

fprintf('\nDELTA-V VALIDATION\n');
fprintf('Integrated powered Delta-V: %.2f m/s\n',dv_powered);
fprintf('Rocket equation Delta-V: %.2f m/s\n',dv_rocket);
fprintf('Difference: %.6f m/s\n',abs(dv_powered-dv_rocket));

fprintf('\nPROPELLANT ACCOUNTING\n');
fprintf('Powered descent propellant: %.2f kg\n',prop_used);
fprintf('Estimated deorbit propellant: %.2f kg\n',prop_deorbit);
fprintf('Total propellant: %.2f kg\n', ...
    prop_used+prop_deorbit);
