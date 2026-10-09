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

%% Site07
h_site = 0.80055;      % km
r_site = R + h_site;

%% Transfer orbit
r1 = R + 100;          % km
a = (r1 + r_site)/2;
v_apo = sqrt(mu*(2/r1 - 1/a));

s0 = [r1; 0; 0; v_apo];

%% Propagate to 15 km altitude
opts = odeset('RelTol',1e-10,'AbsTol',1e-11, ...
    'Events',@(t,s) altitudeEvent(t,s,R,15));

[t_coast,s_coast] = ode45( ...
    @(t,s) coastODE(t,s,mu),[0 4000],s0,opts);

%% Initial powered descent state
r0 = s_coast(end,1:2)' * 1000;   % m
v0 = s_coast(end,3:4)' * 1000;   % m/s

y0 = [r0; v0; m0];

%% Powered descent simulation
params.mu = mu*1e9;    % m^3/s^2
params.R = R*1000;     % m
params.r_site = r_site*1000;
params.Isp = Isp;
params.g0 = g0;
params.Tmax = Tmax;
params.m_dry = m_dry;

opts2 = odeset('RelTol',1e-8,'AbsTol',1e-8, ...
    'Events',@(t,y) touchdownEvent(t,y,params));

[t,y] = ode45(@(t,y) descentODE(t,y,params), ...
    [0 2000],y0,opts2);

%% Extract results
r = vecnorm(y(:,1:2),2,2);
v = vecnorm(y(:,3:4),2,2);
alt = r - params.r_site;
mass = y(:,5);

prop_used = m0 - mass(end);

%% Display results
fprintf('STARSHIP HLS POWERED DESCENT\n\n');
fprintf('Initial altitude: %.2f km\n',alt(1)/1000);
fprintf('Final altitude: %.2f m\n',alt(end));
fprintf('Final speed: %.2f m/s\n',v(end));
fprintf('Descent time: %.2f s\n',t(end));
fprintf('Propellant used: %.2f kg\n',prop_used);

%% Plot results
figure;
tiledlayout(3,1);

nexttile;
plot(t,alt/1000,'LineWidth',1.5);
ylabel('Altitude (km)');
grid on;

nexttile;
plot(t,v,'LineWidth',1.5);
ylabel('Speed (m/s)');
grid on;

nexttile;
plot(t,mass,'LineWidth',1.5);
ylabel('Mass (kg)');
xlabel('Time (s)');
grid on;

%% Functions
function ds = coastODE(~,s,mu)
    r = s(1:2);
    v = s(3:4);
    ds = [v; -mu*r/norm(r)^3];
end

function [value,isterminal,direction] = ...
    altitudeEvent(~,s,R,h)
    value = norm(s(1:2)) - R - h;
    isterminal = 1;
    direction = -1;
end

function dy = descentODE(~,y,p)
    r = y(1:2);
    v = y(3:4);
    m = y(5);

    rmag = norm(r);
    rhat = r/rmag;

    vr = dot(v,rhat);
    vt = v - vr*rhat;

    h = max(rmag-p.r_site,0);

    % Desired radial descent speed
    v_des = -min(40,sqrt(2*1.5*h));

    % Guidance accelerations
    a_rad = 0.15*(v_des-vr);
    a_tan = -0.03*vt;

    g = -p.mu*r/rmag^3;
    a_cmd = a_rad*rhat + a_tan - g;

    Tvec = m*a_cmd;
    Tmag = min(norm(Tvec),p.Tmax);

    if m <= p.m_dry
        Tmag = 0;
    end

    if norm(Tvec) > 0
        u = Tvec/norm(Tvec);
    else
        u = [0;0];
    end

    thrust = Tmag*u;

    dy = [v;
          g + thrust/m;
          -Tmag/(p.Isp*p.g0)];
end

function [value,isterminal,direction] = ...
    touchdownEvent(~,y,p)
    value = norm(y(1:2))-p.r_site;
    isterminal = 1;
    direction = -1;
end
