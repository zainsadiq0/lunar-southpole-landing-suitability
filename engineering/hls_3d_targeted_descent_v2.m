clear; clc; close all;

%% Constants
p.mu = 4.9028e12;       % m^3/s^2
p.R = 1737.4e3;         % m
p.g0 = 9.80665;         % m/s^2

%% Starship HLS assumptions
p.m0 = 250000;          % kg, after deorbit
p.m_dry = 100000;       % kg
p.Isp = 350;            % s
p.Tmax = 2e6;           % N

%% Landing site: Site07
lat = -88.811008;       % deg
lon = 123.690068;       % deg
elev = 800.55;          % m

p.r_site = p.R+elev;

p.site_hat = [cosd(lat)*cosd(lon);
              cosd(lat)*sind(lon);
              sind(lat)];

p.site_pos = p.r_site*p.site_hat;

%% Parking and transfer orbit
r_apo = p.R+100e3;       % m
r_peri = p.r_site;

a = (r_apo+r_peri)/2;

v_circ = sqrt(p.mu/r_apo);
v_apo = sqrt(p.mu*(2/r_apo-1/a));

dv_deorbit = v_circ-v_apo;

% Orbit basis: target direction and tangent
e1 = p.site_hat;
zaxis = [0;0;1];

e2 = cross(zaxis,e1);
e2 = e2/norm(e2);

% At the target, positive e2 is the desired
% direction of orbital motion
% At apoapsis, velocity must be -e2
r_apo_vec = -r_apo*e1;
v_apo_vec = -v_apo*e2;

%% Transfer phase sweep
% Rotates transfer ellipse relative to target
phase_values = -30:5:30;   % deg
n = numel(phase_values);

results = nan(n,5);
success = false(n,1);

best_prop = inf;
best_data = struct();

%% Run cases
for k = 1:n

    phase = phase_values(k);

    % Rotate initial orbit about its normal
    orbit_normal = cross(e1,e2);
    orbit_normal = orbit_normal/norm(orbit_normal);

    Q = axisRotation(orbit_normal,phase);

    r0 = Q*r_apo_vec;
    v0 = Q*v_apo_vec;

    %% Coast to powered descent
    h_start = 14.5e3;    % m above reference sphere
    r_start = p.R+h_start;

    s0 = [r0;v0];

    opts1 = odeset('RelTol',1e-9, ...
        'AbsTol',1e-7, ...
        'Events',@(t,s) coastEvent(t,s,r_start));

    [tc,sc,tec] = ode45( ...
        @(t,s) coastODE(t,s,p.mu), ...
        [0 5000],s0,opts1);

    if isempty(tec)
        continue;
    end

    %% Initial powered state
    y0 = [sc(end,1:3)';
          sc(end,4:6)';
          p.m0;
          0];

    %% Powered descent
    opts2 = odeset('RelTol',1e-8, ...
        'AbsTol',1e-6, ...
        'Events',@(t,y) touchdownEvent(t,y,p));

    [t,y,te] = ode45( ...
        @(t,y) descentODE(t,y,p), ...
        [0 1500],y0,opts2);

    %% Landing results
    rf = y(end,1:3)';
    vf = y(end,4:6)';

    speed = norm(vf);

    angle = acos(max(-1,min(1, ...
        dot(rf/norm(rf),p.site_hat))));

    error_m = p.r_site*angle;

    dv_powered = y(end,8);
    dv_total = dv_deorbit+dv_powered;

    prop_powered = p.m0-y(end,7);

    m_pre = p.m0*exp( ...
        dv_deorbit/(p.Isp*p.g0));

    prop_deorbit = m_pre-p.m0;
    prop_total = prop_powered+prop_deorbit;

    mass_ok = min(y(:,7))>=p.m_dry-1;

    success(k) = ~isempty(te) && ...
        speed<=2 && error_m<=1000 && mass_ok;

    results(k,:) = [phase,speed, ...
        error_m/1000,dv_total,prop_total];

    %% Save best feasible case
    if success(k) && prop_total<best_prop

        best_prop = prop_total;

        best_data.t = t;
        best_data.y = y;
        best_data.tc = tc;
        best_data.sc = sc;
        best_data.phase = phase;
        best_data.speed = speed;
        best_data.error = error_m;
        best_data.dv = dv_total;
        best_data.prop = prop_total;

    end

end

%% Results table
T = array2table(results, ...
    'VariableNames', ...
    {'Phase_deg','Touchdown_mps', ...
     'PositionError_km','DeltaV_mps', ...
     'Propellant_kg'});

T.Success = success;

disp('3D TARGETED DESCENT PHASE SWEEP');
disp(T);

%% Best solution
if isfinite(best_prop)

    fprintf('\nBEST FEASIBLE TRAJECTORY\n');
    fprintf('Phase: %.2f deg\n',best_data.phase);
    fprintf('Touchdown speed: %.3f m/s\n',best_data.speed);
    fprintf('Landing error: %.3f km\n',best_data.error/1000);
    fprintf('Total Delta-V: %.2f m/s\n',best_data.dv);
    fprintf('Propellant: %.2f kg\n',best_data.prop);

else
    fprintf('\nNo feasible trajectory found.\n');
    fprintf('Controller or phase range requires refinement.\n');
end

%% Phase sweep plots
figure;
tiledlayout(2,2);

nexttile;
plot(results(:,1),results(:,2),'-o');
yline(2,'--r','Speed limit');
xlabel('Transfer Phase (deg)');
ylabel('Touchdown Speed (m/s)');
title('Landing Speed');
grid on;

nexttile;
plot(results(:,1),results(:,3),'-o');
yline(1,'--r','Position limit');
xlabel('Transfer Phase (deg)');
ylabel('Position Error (km)');
title('Landing Accuracy');
grid on;

nexttile;
plot(results(:,1),results(:,4),'-o');
xlabel('Transfer Phase (deg)');
ylabel('Delta-V (m/s)');
title('Total Delta-V');
grid on;

nexttile;
plot(results(:,1),results(:,5)/1000,'-o');
xlabel('Transfer Phase (deg)');
ylabel('Propellant (1000 kg)');
title('Propellant Consumption');
grid on;

%% Best trajectory plots
if isfinite(best_prop)

    t = best_data.t;
    y = best_data.y;

    figure;
    tiledlayout(2,2);

    altitude = ...
        (vecnorm(y(:,1:3),2,2)-p.r_site)/1000;

    velocity = vecnorm(y(:,4:6),2,2);

    nexttile;
    plot(t,altitude,'LineWidth',1.5);
    xlabel('Time (s)');
    ylabel('Altitude (km)');
    title('Altitude');
    grid on;

    nexttile;
    plot(t,velocity,'LineWidth',1.5);
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
    title('Powered Delta-V');
    grid on;

    %% 3D plot
    figure;
    hold on;

    [X,Y,Z] = sphere(80);

    surf(p.R*X/1000,p.R*Y/1000,p.R*Z/1000, ...
        'FaceColor',[0.6 0.6 0.6], ...
        'EdgeColor','none','FaceAlpha',0.8);

    sc = best_data.sc;

    plot3(sc(:,1)/1000,sc(:,2)/1000, ...
        sc(:,3)/1000,'b','LineWidth',1.5);

    plot3(y(:,1)/1000,y(:,2)/1000, ...
        y(:,3)/1000,'r','LineWidth',2);

    plot3(p.site_pos(1)/1000, ...
        p.site_pos(2)/1000, ...
        p.site_pos(3)/1000, ...
        'go','MarkerFaceColor','g');

    axis equal;
    grid on;
    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    title('Best 3D Targeted Lunar Descent');
    legend('Moon','Coast','Powered Descent', ...
        'Target','Location','bestoutside');

    view(40,-30);
    camlight;
    lighting gouraud;

end

%% Functions
function Q = axisRotation(axis,angle)

    axis = axis/norm(axis);
    x = axis(1);
    y = axis(2);
    z = axis(3);

    c = cosd(angle);
    s = sind(angle);
    C = 1-c;

    Q = [c+x*x*C, x*y*C-z*s, x*z*C+y*s;
         y*x*C+z*s, c+y*y*C, y*z*C-x*s;
         z*x*C-y*s, z*y*C+x*s, c+z*z*C];
end

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

    %% Vertical guidance
    v_des = -min(40,sqrt(2*0.8*h));

    if h>1000
        Kvr = 0.15;
    elseif h>100
        Kvr = 0.4;
    else
        Kvr = 1.0;
    end

    a_rad = Kvr*(v_des-vr);

    %% Horizontal target guidance
    error_vec = p.site_pos-r;

    error_horizontal = ...
        error_vec-dot(error_vec,rhat)*rhat;

    % Desired horizontal velocity
    Kpos = 0.002;

    v_horizontal_des = Kpos*error_horizontal;

    vmax = 100;          % m/s

    if norm(v_horizontal_des)>vmax
        v_horizontal_des = ...
            vmax*v_horizontal_des/ ...
            norm(v_horizontal_des);
    end

    % Horizontal velocity control
    Kv = 0.08;

    a_horizontal = Kv*(v_horizontal_des-vt);

    %% Gravity compensation
    gravity = -p.mu*r/rmag^3;

    a_cmd = a_rad*rhat+ ...
            a_horizontal-gravity;

    %% Thrust limits
    Tvec = m*a_cmd;

    Tmag = min(norm(Tvec),p.Tmax);

    if m<=p.m_dry
        Tmag = 0;
    end

    if norm(Tvec)>0
        thrust = Tmag*Tvec/norm(Tvec);
    else
        thrust = zeros(3,1);
    end

    %% Dynamics
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
