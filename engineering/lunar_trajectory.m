clear; clc; close all;

%% Lunar constants
mu = 4902.8;           % km^3/s^2
R_moon = 1737.4;       % km

%% Initial orbit and landing site
h_orbit = 100;         % km
h_site = 0.80055;      % km

r1 = R_moon + h_orbit;
r2 = R_moon + h_site;
a = (r1 + r2)/2;

%% Initial conditions at transfer apoapsis
v_apo = sqrt(mu*(2/r1 - 1/a));

x0 = r1;               % km
y0 = 0;                % km
vx0 = 0;               % km/s
vy0 = v_apo;           % km/s

state0 = [x0; y0; vx0; vy0];

%% Transfer time
t_transfer = pi*sqrt(a^3/mu);

%% Numerical integration
options = odeset('RelTol',1e-11,'AbsTol',1e-12);

[t,state] = ode45(@(t,s) orbitalODE(t,s,mu), ...
    [0 t_transfer],state0,options);

%% Extract trajectory
x = state(:,1);
y = state(:,2);
vx = state(:,3);
vy = state(:,4);

r = sqrt(x.^2 + y.^2);
altitude = r - R_moon;
speed = sqrt(vx.^2 + vy.^2);

%% Results
fprintf('NUMERICAL TRAJECTORY RESULTS\n\n');
fprintf('Final altitude: %.4f km\n',altitude(end));
fprintf('Final speed: %.4f km/s\n',speed(end));
fprintf('Transfer time: %.2f minutes\n',t(end)/60);

%% Plot trajectory
figure;
theta = linspace(0,2*pi,500);

plot(R_moon*cos(theta),R_moon*sin(theta), ...
    'k','LineWidth',1.5);
hold on;

plot(x,y,'b','LineWidth',2);
plot(x(1),y(1),'go','MarkerFaceColor','g');
plot(x(end),y(end),'ro','MarkerFaceColor','r');

axis equal;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
title('Lunar Deorbit Transfer Trajectory');
legend('Moon','Trajectory','Start','Periapsis');

%% Plot altitude
figure;
plot(t/60,altitude,'LineWidth',2);
grid on;
xlabel('Time (min)');
ylabel('Altitude (km)');
title('Altitude During Lunar Transfer');

%% Orbital dynamics
function ds = orbitalODE(~,s,mu)
r = s(1:2);
v = s(3:4);

rmag = norm(r);
acceleration = -mu*r/rmag^3;

ds = [v; acceleration];
end
