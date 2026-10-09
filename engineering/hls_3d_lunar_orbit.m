clear; clc; close all;

%% Lunar constants
mu = 4902.8;           % km^3/s^2
R = 1737.4;            % km
h_orbit = 100;         % km
r_orbit = R+h_orbit;   % km

%% Initial lunar orbit
inc = 90;              % deg
RAAN = 0;              % deg

v_circ = sqrt(mu/r_orbit);
T_orbit = 2*pi*sqrt(r_orbit^3/mu);

%% Landing site: Site07
site_lat = -88.811008; % deg
site_lon = 123.690068; % deg
h_site = 0.80055;      % km

r_site = R+h_site;

site_xyz = r_site * [ ...
    cosd(site_lat)*cosd(site_lon);
    cosd(site_lat)*sind(site_lon);
    sind(site_lat)];

site_unit = site_xyz/norm(site_xyz);

%% Generate 3D circular orbit
u = linspace(0,360,500);

x = r_orbit*( ...
    cosd(RAAN).*cosd(u) - ...
    sind(RAAN).*sind(u).*cosd(inc));

y = r_orbit*( ...
    sind(RAAN).*cosd(u) + ...
    cosd(RAAN).*sind(u).*cosd(inc));

z = r_orbit*sind(u).*sind(inc);

%% Orbital plane normal
n = [ ...
    sind(inc)*sind(RAAN);
   -sind(inc)*cosd(RAAN);
    cosd(inc)];

n = n/norm(n);

%% Site-to-plane angular separation
delta = asind(abs(dot(n,site_unit)));

%% Illustrative plane-change Delta-V
dv_plane = 2*v_circ*sind(delta/2)*1000;

%% Closest orbital ground-track point
site_projected = site_unit-dot(site_unit,n)*n;
site_projected = site_projected/norm(site_projected);

closest_xyz = R*site_projected;

%% Results
fprintf('3D LUNAR ORBIT ANALYSIS\n\n');

fprintf('Orbit inclination: %.2f deg\n',inc);
fprintf('RAAN: %.2f deg\n',RAAN);
fprintf('Orbital period: %.2f minutes\n',T_orbit/60);
fprintf('Circular velocity: %.2f m/s\n',v_circ*1000);

fprintf('\nLANDING SITE\n');
fprintf('Latitude: %.6f deg\n',site_lat);
fprintf('Longitude: %.6f deg\n',site_lon);

fprintf('\nGEOMETRIC ANALYSIS\n');
fprintf('Site-to-plane angle: %.4f deg\n',delta);
fprintf('Illustrative plane-change Delta-V: %.2f m/s\n',dv_plane);

%% Plot Moon
figure;
hold on;

[X,Y,Z] = sphere(100);

surf(R*X,R*Y,R*Z, ...
    'FaceColor',[0.65 0.65 0.65], ...
    'EdgeColor','none', ...
    'FaceAlpha',0.85);

%% Plot orbit
plot3(x,y,z,'b','LineWidth',2);

%% Plot landing site
plot3(site_xyz(1),site_xyz(2),site_xyz(3), ...
    'ro','MarkerFaceColor','r','MarkerSize',8);

%% Plot closest ground-track point
plot3(closest_xyz(1),closest_xyz(2),closest_xyz(3), ...
    'go','MarkerFaceColor','g','MarkerSize',7);

%% Connect site to closest point
plot3([site_xyz(1) closest_xyz(1)], ...
      [site_xyz(2) closest_xyz(2)], ...
      [site_xyz(3) closest_xyz(3)], ...
      'k--','LineWidth',1.5);

axis equal;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('3D Lunar Polar Orbit and Artemis Landing Site');

legend('Moon','Polar Orbit','Site07', ...
    'Closest Ground Track','Angular Separation', ...
    'Location','bestoutside');

view(40,25);
camlight;
lighting gouraud;
