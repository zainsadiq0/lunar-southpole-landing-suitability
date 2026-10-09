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

%% Initial transfer orbit
r1 = R + 100;          % km
a = (r1 + r_site)/2;
v_apo = sqrt(mu*(2/r1 - 1/a));

s0 = [r1; 0; 0; v_apo];

dv_deorbit = (sqrt(mu/r1)-v_apo)*1000;

% m0 is the mass immediately after deorbit
m_before_deorbit = m0*exp(dv_deorbit/(Isp*g0));
prop_deorbit = m_before_deorbit-m0;

%% Simulation parameters
p.mu = mu*1e9;         % m^3/s^2
p.R = R*1000;          % m
p.r_site = r_site*1000;
p.Isp = Isp;
p.g0 = g0;
p.Tmax = Tmax;
p.m_dry = m_dry;

%% Powered descent starting altitudes
h_start_values = 10:0.5:16; % km

n = length(h_start_values);

dv_total = nan(n,1);
prop_total = nan(n,1);
touchdown_speed = nan(n,1);
descent_time = nan(n,1);
success = false(n,1);

%% Run simulations
for i = 1:n

    h_start = h_start_values(i);

    %% Coast to starting altitude
    opts = odeset('RelTol',1e-10, ...
        'AbsTol',1e-11, ...
        'Events',@(t,s) altitudeEvent(t,s,R,h_start));

    [~,s_coast,te_coast] = ode45( ...
        @(t,s) coastODE(t,s,mu), ...
        [0 4000],s0,opts);

    if isempty(te_coast)
        fprintf('Coast failed for %.0f km\n',h_start);
        continue;
    end

    %% Initial powered descent state
    r0 = s_coast(end,1:2)'*1000;
    v0 = s_coast(end,3:4)'*1000;

    y0 = [r0; v0; m0; 0];

    %% Powered descent
    opts2 = odeset('RelTol',1e-8, ...
        'AbsTol',1e-8, ...
        'Events',@(t,y) touchdownEvent(t,y,p));

    [t,y,te_touch] = ode45( ...
        @(t,y) descentODE(t,y,p), ...
        [0 2000],y0,opts2);

    %% Results
    r_final = y(end,1:2)';
    v_final = y(end,3:4)';

    rhat = r_final/norm(r_final);

    vr = dot(v_final,rhat);
    vt = norm(v_final-vr*rhat);

    touchdown_speed(i) = norm(v_final);
    descent_time(i) = t(end);

    dv_total(i) = dv_deorbit+y(end,6);

    prop_powered = m0-y(end,5);
    prop_total(i) = prop_deorbit+prop_powered;

    % Check landing requirements
    success(i) = ~isempty(te_touch) && ...
        abs(vr)<=2 && vt<=1 && ...
        touchdown_speed(i)<=2 && ...
        min(y(:,5))>=m_dry-1;

end

%% Display results
results = table( ...
    h_start_values', ...
    dv_total, ...
    prop_total, ...
    touchdown_speed, ...
    descent_time, ...
    success, ...
    'VariableNames', ...
    {'StartAlt_km','DeltaV_mps','Propellant_kg', ...
     'TouchdownSpeed_mps','DescentTime_s','Success'});

disp('POWERED DESCENT OPTIMIZATION RESULTS');
disp(results);

%% Find best feasible altitude
valid = find(success);

if ~isempty(valid)

    [~,idx] = min(prop_total(valid));
    best = valid(idx);

    fprintf('\nBEST FEASIBLE STARTING ALTITUDE\n');
    fprintf('Altitude: %.0f km\n',h_start_values(best));
    fprintf('Total Delta-V: %.2f m/s\n',dv_total(best));
    fprintf('Total propellant: %.2f kg\n',prop_total(best));
    fprintf('Touchdown speed: %.2f m/s\n',touchdown_speed(best));

else
    fprintf('\nNo successful landings found.\n');
end

%% Plot comparisons
figure;
tiledlayout(2,2);

nexttile;
plot(h_start_values(success), ...
    prop_total(success)/1000, ...
    '-o','LineWidth',1.5);
hold on;
plot(h_start_values(~success), ...
    prop_total(~success)/1000, ...
    'rx','MarkerSize',9,'LineWidth',1.5);

xlabel('Starting Altitude (km)');
ylabel('Propellant (1000 kg)');
title('Propellant vs Starting Altitude');
legend('Feasible','Infeasible','Location','best');
grid on;

nexttile;
plot(h_start_values,dv_total, ...
    '-o','LineWidth',1.5);
xlabel('Starting Altitude (km)');
ylabel('Delta-V (m/s)');
title('Delta-V vs Starting Altitude');
grid on;

nexttile;
plot(h_start_values,touchdown_speed, ...
    '-o','LineWidth',1.5);
yline(2,'--r','Touchdown Limit');
xlabel('Starting Altitude (km)');
ylabel('Touchdown Speed (m/s)');
title('Landing Velocity');
grid on;

nexttile;
plot(h_start_values,descent_time, ...
    '-o','LineWidth',1.5);
xlabel('Starting Altitude (km)');
ylabel('Descent Time (s)');
title('Descent Duration');
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

    %% Controller gains
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

    %% Thrust
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
