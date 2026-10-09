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

%% Site07 coordinates
p.lat = -88.811008;     % deg
p.lon = 123.690068;     % deg
p.elev = 800.55;        % m

p.r_site = p.R+p.elev;

p.site_hat = [cosd(p.lat)*cosd(p.lon);
              cosd(p.lat)*sind(p.lon);
              sind(p.lat)];

p.site_pos = p.r_site*p.site_hat;

%% Initial orbit
h_orbit = 100e3;        % m
r_apo = p.R+h_orbit;
r_peri = p.r_site;

a = (r_apo+r_peri)/2;

v_circ = sqrt(p.mu/r_apo);
v_apo = sqrt(p.mu*(2/r_apo-1/a));

dv_deorbit = v_circ-v_apo;

%% Orbit basis
e1 = p.site_hat;
zaxis = [0;0;1];

e2 = cross(zaxis,e1);
e2 = e2/norm(e2);

normal = cross(e1,e2);
normal = normal/norm(normal);

%% Parameter sweep
% Powered descent initiation altitude and transfer phase
altitudes_km = [14.5 20 30 40 50];
phases_deg = -60:10:60;

N = numel(altitudes_km)*numel(phases_deg);

results = nan(N,7);
best = struct();
best_prop = inf;

case_id = 0;

for ia = 1:numel(altitudes_km)
    for ip = 1:numel(phases_deg)

        case_id = case_id+1;

        h_start = altitudes_km(ia)*1000;
        phase = phases_deg(ip);

        % Transfer orbit rotation
        Q = axisRotation(normal,phase);

        r0 = Q*(-r_apo*e1);
        v0 = Q*(-v_apo*e2);

        %% Coast propagation
        r_start = p.R+h_start;

        opts1 = odeset('RelTol',1e-9, ...
            'AbsTol',1e-7, ...
            'Events',@(t,s) coastEvent(t,s,r_start));

        [tc,sc,tec] = ode45( ...
            @(t,s) coastODE(t,s,p.mu), ...
            [0 5000],[r0;v0],opts1);

        if isempty(tec)
            continue;
        end

        %% Powered descent initial state
        y0 = [sc(end,1:3)';
              sc(end,4:6)';
              p.m0;
              0];

        %% Integrate powered descent
        opts2 = odeset('RelTol',1e-8, ...
            'AbsTol',1e-6, ...
            'Events',@(t,y) stopEvent(t,y,p));

        [t,y,te,~,ie] = ode45( ...
            @(t,y) descentODE(t,y,p), ...
            [0 2000],y0,opts2);

        %% Final conditions
        rf = y(end,1:3)';
        vf = y(end,4:6)';

        speed = norm(vf);

        angle = acos(max(-1,min(1, ...
            dot(rf/norm(rf),p.site_hat))));

        error_km = p.r_site*angle/1000;

        dv_powered = y(end,8);
        dv_total = dv_deorbit+dv_powered;

        m_pre = p.m0*exp(dv_deorbit/(p.Isp*p.g0));
        prop_total = m_pre-y(end,7);

        surface_reached = any(ie==1);
        mass_ok = min(y(:,7))>=p.m_dry-1;

        feasible = surface_reached && mass_ok && ...
            speed<=2 && error_km<=1;

        results(case_id,:) = ...
            [h_start/1000,phase,speed,error_km, ...
             dv_total,prop_total,double(feasible)];

        if feasible && prop_total<best_prop

            best_prop = prop_total;

            best.t = t;
            best.y = y;
            best.tc = tc;
            best.sc = sc;
            best.alt = h_start/1000;
            best.phase = phase;
            best.speed = speed;
            best.error = error_km;
            best.dv = dv_total;
            best.prop = prop_total;

        end
    end
end

%% Display results
T = array2table(results, ...
    'VariableNames', ...
    {'StartAlt_km','Phase_deg','Touchdown_mps', ...
     'Error_km','DeltaV_mps','Propellant_kg','Success'});

disp('3D TARGETED DESCENT V3 RESULTS');
disp(T);

if isfinite(best_prop)

    fprintf('\nBEST FEASIBLE LANDING\n');
    fprintf('Start altitude: %.2f km\n',best.alt);
    fprintf('Transfer phase: %.2f deg\n',best.phase);
    fprintf('Touchdown speed: %.3f m/s\n',best.speed);
    fprintf('Landing error: %.3f km\n',best.error);
    fprintf('Total Delta-V: %.2f m/s\n',best.dv);
    fprintf('Total propellant: %.2f kg\n',best.prop);

else
    fprintf('\nNO FEASIBLE LANDING FOUND\n');
    fprintf('The tested guidance and trajectory ');
    fprintf('settings did not satisfy all constraints.\n');
end

%% Parameter sweep visualization
figure;
tiledlayout(2,2);

nexttile;
scatter(results(:,1),results(:,3),45, ...
    results(:,2),'filled');
yline(2,'--r','Touchdown limit');
xlabel('Start Altitude (km)');
ylabel('Touchdown Speed (m/s)');
title('Touchdown Speed');
colorbar; grid on;

nexttile;
scatter(results(:,1),results(:,4),45, ...
    results(:,2),'filled');
yline(1,'--r','Position limit');
xlabel('Start Altitude (km)');
ylabel('Landing Error (km)');
title('Landing Accuracy');
colorbar; grid on;

nexttile;
scatter(results(:,1),results(:,5),45, ...
    results(:,2),'filled');
xlabel('Start Altitude (km)');
ylabel('Total Delta-V (m/s)');
title('Delta-V');
colorbar; grid on;

nexttile;
scatter(results(:,1),results(:,6)/1000,45, ...
    results(:,2),'filled');
xlabel('Start Altitude (km)');
ylabel('Propellant (1000 kg)');
title('Propellant Consumption');
colorbar; grid on;

%% Best trajectory
if isfinite(best_prop)

    t = best.t;
    y = best.y;

    altitude = ...
        (vecnorm(y(:,1:3),2,2)-p.r_site)/1000;

    speed = vecnorm(y(:,4:6),2,2);

    figure;
    tiledlayout(2,2);

    nexttile;
    plot(t,altitude,'LineWidth',1.5);
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
    ylabel('Powered Delta-V (m/s)');
    title('Delta-V');
    grid on;

    %% 3D trajectory
    figure;
    hold on;

    [X,Y,Z] = sphere(80);

    surf(p.R*X/1000,p.R*Y/1000,p.R*Z/1000, ...
        'FaceColor',[0.6 0.6 0.6], ...
        'EdgeColor','none', ...
        'FaceAlpha',0.8);

    sc = best.sc;

    plot3(sc(:,1)/1000,sc(:,2)/1000, ...
        sc(:,3)/1000,'b','LineWidth',1.5);

    plot3(y(:,1)/1000,y(:,2)/1000, ...
        y(:,3)/1000,'r','LineWidth',2);

    plot3(p.site_pos(1)/1000, ...
        p.site_pos(2)/1000, ...
        p.site_pos(3)/1000, ...
        'go','MarkerFaceColor','g','MarkerSize',8);

    plot3(y(end,1)/1000,y(end,2)/1000, ...
        y(end,3)/1000, ...
        'mo','MarkerFaceColor','m','MarkerSize',8);

    axis equal;
    grid on;
    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    title('Best 3D Starship HLS Descent');
    legend('Moon','Coast','Powered Descent', ...
        'Target','Touchdown','Location','bestoutside');

    view(40,-30);
    camlight;
    lighting gouraud;

end

%% Local functions
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

    h = max(rmag-p.r_site,0);

    vr = dot(v,rhat);
    vh = v-vr*rhat;

    %% Horizontal target geometry
    target = p.site_pos;

    err = target-r;

    err_h = err-dot(err,rhat)*rhat;
    distance = norm(err_h);

    %% Time-to-go estimate
    % Based on altitude and descent velocity
    tgo = max(20,min(600, ...
        2*h/max(abs(vr),20)));

    %% Desired horizontal velocity
    % Reduce allowable speed near the target
    a_brake = 3.0;          % m/s^2

    v_allow = sqrt(2*a_brake*distance);

    if distance>1
        v_des_h = min(v_allow,distance/tgo)* ...
            err_h/distance;
    else
        v_des_h = zeros(3,1);
    end

    %% Horizontal acceleration command
    Kh = 0.12;
    a_h = Kh*(v_des_h-vh);

    %% Vertical descent guidance
    % Desired speed approaches zero at touchdown
    a_vertical = 1.0;

    v_des_r = -min(40,sqrt(2*a_vertical*h));

    Kr = 0.5;

    a_r = Kr*(v_des_r-vr);

    %% Gravity compensation
    gravity = -p.mu*r/rmag^3;

    a_cmd = a_h+a_r*rhat-gravity;

    %% Thrust saturation
    Tvec = m*a_cmd;

    if norm(Tvec)>p.Tmax
        Tvec = p.Tmax*Tvec/norm(Tvec);
    end

    if m<=p.m_dry
        Tvec = zeros(3,1);
    end

    Tmag = norm(Tvec);

    %% State derivatives
    dy = [v;
          gravity+Tvec/m;
          -Tmag/(p.Isp*p.g0);
          Tmag/m];
end

function [value,isterminal,direction] = ...
    stopEvent(~,y,p)

    surface = norm(y(1:3))-p.r_site;
    fuel = y(7)-p.m_dry;

    value = [surface;fuel];
    isterminal = [1;1];
    direction = [-1;-1];
end
