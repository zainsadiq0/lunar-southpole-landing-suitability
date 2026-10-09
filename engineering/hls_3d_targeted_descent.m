clear; clc; close all;

%% Constants (SI units)
mu = 4.9028e12;          % m^3/s^2
R = 1737.4e3;            % m
g0 = 9.80665;            % m/s^2

%% Starship HLS assumptions
m0 = 250000;             % kg, after deorbit
m_dry = 100000;          % kg
Isp = 350;               % s
Tmax = 2e6;              % N

%% Site07: Peak Near Shackleton
site_lat = -88.811008;   % deg
site_lon = 123.690068;   % deg
site_elev = 0.80055e3;   % m

r_site = R + site_elev;

site_hat = [cosd(site_lat)*cosd(site_lon);
            cosd(site_lat)*sind(site_lon);
            sind(site_lat)];

site_pos = r_site*site_hat;

%% Reference parking orbit and transfer
h_orbit = 100e3;         % m
r_apo = R + h_orbit;
r_peri = r_site;

a = (r_apo+r_peri)/2;

v_circ = sqrt(mu/r_apo);
v_apo = sqrt(mu*(2/r_apo-1/a));

dv_deorbit = v_circ-v_apo;

% Orbital plane containing site and lunar Z axis
zaxis = [0;0;1];
q = cross(zaxis,site_hat);
q = q/norm(q);

% Set transfer apoapsis opposite the target
r0 = -r_apo*site_hat;

% Choose transfer velocity direction
v0 = v_apo*q;

% Verify orientation: orbit must approach site
% in the selected direction
hvec = cross(r0,v0);
if dot(hvec,cross(site_hat,q)) < 0
    v0 = -v0;
end

%% Coast to powered-descent altitude
h_start = 14.5e3;       % m above reference sphere
r_start = R+h_start;

s0 = [r0;v0];

opts = odeset('RelTol',1e-10, ...
              'AbsTol',1e-8, ...
              'Events',@(t,s) coastEvent(t,s,r_start));

[t_coast,s_coast,te_coast] = ode45( ...
    @(t,s) coastODE(t,s,mu), ...
    [0 5000],s0,opts);

if isempty(te_coast)
    error('Coast did not reach powered-descent altitude.');
end

%% Coast geometry validation
r_coast = s_coast(end,1:3)';
v_coast = s_coast(end,4:6)';

coast_hat = r_coast/norm(r_coast);

coast_angle = acosd(max(-1,min(1, ...
    dot(coast_hat,site_hat))));

coast_vr = dot(v_coast,coast_hat);
coast_vt = norm(v_coast-coast_vr*coast_hat);

fprintf('\nCOAST GEOMETRY VALIDATION\n');
fprintf('Angular distance to target: %.3f deg\n', ...
    coast_angle);
fprintf('Radial velocity: %.2f m/s\n',coast_vr);
fprintf('Tangential velocity: %.2f m/s\n',coast_vt);

%% Powered descent initial state
y0 = [s_coast(end,1:3)';
      s_coast(end,4:6)';
      m0;
      0];

%% Powered descent initial state
y0 = [s_coast(end,1:3)';
      s_coast(end,4:6)';
      m0;
      0];

p.mu = mu;
p.r_site = r_site;
p.site_hat = site_hat;
p.Isp = Isp;
p.g0 = g0;
p.Tmax = Tmax;
p.m_dry = m_dry;

%% Powered descent integration
opts2 = odeset('RelTol',1e-8, ...
               'AbsTol',1e-7, ...
               'Events',@(t,y) touchdownEvent(t,y,p));

[t,y,te_touch] = ode45( ...
    @(t,y) descentODE(t,y,p), ...
    [0 1500],y0,opts2);

%% Results
rf = y(end,1:3)';
vf = y(end,4:6)';

rhat = rf/norm(rf);

vr = dot(vf,rhat);
vt = norm(vf-vr*rhat);

touchdown_speed = norm(vf);

landing_angle = acos( ...
    max(-1,min(1,dot(rhat,site_hat))));

landing_error = r_site*landing_angle;

dv_powered = y(end,8);
dv_total = dv_deorbit+dv_powered;

prop_powered = m0-y(end,7);

% Account for initial deorbit maneuver
m_pre = m0*exp(dv_deorbit/(Isp*g0));
prop_deorbit = m_pre-m0;
prop_total = prop_powered+prop_deorbit;

%% Output
fprintf('\n3D SITE-TARGETED DESCENT RESULTS\n');
fprintf('Site: Site07 - Peak Near Shackleton\n\n');

fprintf('Coast duration: %.2f min\n',t_coast(end)/60);
fprintf('Powered descent duration: %.2f s\n',t(end));

fprintf('\nTOUCHDOWN\n');
fprintf('Surface reached: %d\n',~isempty(te_touch));
fprintf('Radial velocity: %.3f m/s\n',vr);
fprintf('Tangential velocity: %.3f m/s\n',vt);
fprintf('Touchdown speed: %.3f m/s\n',touchdown_speed);
fprintf('Landing position error: %.3f km\n', ...
    landing_error/1000);

fprintf('\nPROPULSION\n');
fprintf('Deorbit Delta-V: %.2f m/s\n',dv_deorbit);
fprintf('Powered Delta-V: %.2f m/s\n',dv_powered);
fprintf('Total Delta-V: %.2f m/s\n',dv_total);
fprintf('Total propellant: %.2f kg\n',prop_total);

% Preliminary landing criteria
landed = ~isempty(te_touch);
soft = touchdown_speed <= 2;
accurate = landing_error <= 1000;
mass_ok = min(y(:,7)) >= m_dry-1;

if landed && soft && accurate && mass_ok
    fprintf('\nTARGETED LANDING ACHIEVED\n');
else
    fprintf('\nTARGETED LANDING NOT YET ACHIEVED\n');
end

%% 3D trajectory plot
figure;
hold on;

[X,Y,Z] = sphere(100);

surf(R*X/1000,R*Y/1000,R*Z/1000, ...
    'FaceColor',[0.6 0.6 0.6], ...
    'EdgeColor','none', ...
    'FaceAlpha',0.8);

plot3(s_coast(:,1)/1000, ...
      s_coast(:,2)/1000, ...
      s_coast(:,3)/1000, ...
      'b','LineWidth',1.5);

plot3(y(:,1)/1000, ...
      y(:,2)/1000, ...
      y(:,3)/1000, ...
      'r','LineWidth',2);

plot3(site_pos(1)/1000, ...
      site_pos(2)/1000, ...
      site_pos(3)/1000, ...
      'go','MarkerFaceColor','g','MarkerSize',8);

plot3(rf(1)/1000,rf(2)/1000,rf(3)/1000, ...
      'mo','MarkerFaceColor','m','MarkerSize',7);

axis equal;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('3D Starship HLS Lunar Descent');

legend('Moon','Coast','Powered Descent', ...
       'Target Site','Final Position', ...
       'Location','bestoutside');

view(40,-30);
camlight;
lighting gouraud;

%% Descent performance
figure;
tiledlayout(2,2);

alt = (vecnorm(y(:,1:3),2,2)-r_site)/1000;
speed = vecnorm(y(:,4:6),2,2);

nexttile;
plot(t,alt,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Altitude (km)');
title('Altitude');
grid on;

nexttile;
plot(t,speed,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Speed (m/s)');
title('Velocity');
grid on;

nexttile;
plot(t,y(:,7),'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Mass (kg)');
title('Spacecraft Mass');
grid on;

nexttile;
plot(t,y(:,8),'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Delta-V (m/s)');
title('Accumulated Powered Delta-V');
grid on;

%% Local functions
function ds = coastODE(~,s,mu)
    r = s(1:3);
    v = s(4:6);

    ds = [v;-mu*r/norm(r)^3];
end

function [value,isterminal,direction] = ...
    coastEvent(~,s,r_start)

    value = norm(s(1:3))-r_start;
    isterminal = 1;
    direction = -1;
end

function dy = descentODE(~,y,p)

    r = y(1:3);
    v = y(4:6);
    m = y(7);

    rmag = norm(r);
    rhat = r/rmag;

    vr = dot(v,rhat);
    vt = v-vr*rhat;

    h = max(rmag-p.r_site,0);

    %% Vertical velocity guidance
    v_des = -min(40,sqrt(2*0.8*h));

    if h > 1000
        Kp = 0.15;
    elseif h > 100
        Kp = 0.4;
    else
        Kp = 1.0;
    end

    %% Target-position correction
    target_direction = p.site_hat;

    position_error = ...
        p.r_site*target_direction-r;

    lateral_error = position_error- ...
        dot(position_error,rhat)*rhat;

    % Modest lateral correction
    a_position = 0.00005*lateral_error;

    %% Acceleration command
    a_rad = Kp*(v_des-vr);
    a_tan = -0.03*vt;

    gravity = -p.mu*r/rmag^3;

    a_cmd = a_rad*rhat+a_tan+ ...
            a_position-gravity;

    Tvec = m*a_cmd;
    Tmag = min(norm(Tvec),p.Tmax);

    if m <= p.m_dry
        Tmag = 0;
    end

    if norm(Tvec)>0
        u = Tvec/norm(Tvec);
    else
        u = [0;0;0];
    end

    thrust = Tmag*u;

    dy = [v;
          gravity+thrust/m;
          -Tmag/(p.Isp*p.g0);
          Tmag/m];
end

function [value,isterminal,direction] = ...
    touchdownEvent(~,y,p)

    value = norm(y(1:3))-p.r_site;
    isterminal = 1;
    direction = -1;
end
