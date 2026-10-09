clear; clc; close all;

%% Lunar constants
mu = 4902.8;           % km^3/s^2
R = 1737.4;            % km
h_orbit = 100;         % km
r_orbit = R+h_orbit;

v_circ = sqrt(mu/r_orbit);

%% Initial parking orbit
inc = 90;              % deg
RAAN = 0;              % deg

%% Landing site data
% Site07 is from the supplied dataset.
% Other sites are hypothetical test locations.

site_names = ["Site07";"Test A";"Test B"; ...
              "Test C";"Test D"];

site_lat = [-88.811008;-89.2;-87.5;-85.0;-88.0];
site_lon = [123.690068;45;90;180;270];

site_elev = [0.80055;0;0;0;0]; % km

n_sites = length(site_lat);

%% Orbital plane normal
n = [sind(inc)*sind(RAAN);
    -sind(inc)*cosd(RAAN);
     cosd(inc)];

n = n/norm(n);

%% Calculate site positions and geometry
site_xyz = zeros(3,n_sites);
plane_angle = zeros(n_sites,1);
dv_plane = zeros(n_sites,1);

for i = 1:n_sites

    r_site = R+site_elev(i);

    site_xyz(:,i) = r_site*[ ...
        cosd(site_lat(i))*cosd(site_lon(i));
        cosd(site_lat(i))*sind(site_lon(i));
        sind(site_lat(i))];

    site_unit = site_xyz(:,i)/norm(site_xyz(:,i));

    plane_angle(i) = asind( ...
        abs(dot(n,site_unit)));

    dv_plane(i) = 2*v_circ* ...
        sind(plane_angle(i)/2)*1000;

end

%% Results table
results = table(site_names,site_lat,site_lon, ...
    plane_angle,dv_plane, ...
    'VariableNames', ...
    {'Site','Latitude_deg','Longitude_deg', ...
     'PlaneAngle_deg','PlaneDV_mps'});

disp('LUNAR LANDING SITE GEOMETRY');
disp(results);

%% RAAN sensitivity study
RAAN_values = 0:2:180;

dv_RAAN = zeros(length(RAAN_values),n_sites);

for j = 1:length(RAAN_values)

    O = RAAN_values(j);

    normal = [sind(inc)*sind(O);
             -sind(inc)*cosd(O);
              cosd(inc)];

    for i = 1:n_sites

        site_unit = site_xyz(:,i)/norm(site_xyz(:,i));

        delta = asind(abs(dot(normal,site_unit)));

        dv_RAAN(j,i) = 2*v_circ* ...
            sind(delta/2)*1000;

    end
end

%% Display minimum plane alignment estimate
min_dv = min(dv_RAAN,[],1);

fprintf('\nMINIMUM PLANE-ALIGNMENT ESTIMATES\n');

for i = 1:n_sites
    fprintf('%s: %.2f m/s\n',site_names(i),min_dv(i));
end

%% 3D Moon and orbit
figure;
hold on;

[X,Y,Z] = sphere(100);

surf(R*X,R*Y,R*Z, ...
    'FaceColor',[0.65 0.65 0.65], ...
    'EdgeColor','none', ...
    'FaceAlpha',0.85);

u = linspace(0,360,500);

x = r_orbit*( ...
    cosd(RAAN).*cosd(u)- ...
    sind(RAAN).*sind(u).*cosd(inc));

y = r_orbit*( ...
    sind(RAAN).*cosd(u)+ ...
    cosd(RAAN).*sind(u).*cosd(inc));

z = r_orbit*sind(u).*sind(inc);

plot3(x,y,z,'b','LineWidth',2);

colors = lines(n_sites);

for i = 1:n_sites
    plot3(site_xyz(1,i),site_xyz(2,i), ...
        site_xyz(3,i),'o', ...
        'MarkerFaceColor',colors(i,:), ...
        'MarkerEdgeColor','k', ...
        'MarkerSize',8);
end

axis equal;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('Lunar Landing Sites and Polar Orbit');

legend(['Moon';'Polar Orbit';site_names], ...
    'Location','bestoutside');

view(40,-35);
camlight;
lighting gouraud;

%% Plane-change estimate vs RAAN
figure;
hold on;

for i = 1:n_sites
    plot(RAAN_values,dv_RAAN(:,i), ...
        'LineWidth',1.5, ...
        'DisplayName',site_names(i));
end

xlabel('RAAN (deg)');
ylabel('Illustrative Plane-Change Delta-V (m/s)');
title('Orbital Plane Alignment vs RAAN');
legend('Location','best');
grid on;

%% Fixed orbit comparison
figure;

bar(categorical(site_names),dv_plane);

xlabel('Landing Site');
ylabel('Illustrative Plane-Change Delta-V (m/s)');
title('Site Comparison for Fixed RAAN');
grid on;
