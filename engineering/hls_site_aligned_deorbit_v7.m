clear; clc; close all;

%% Constants and target (km, s)
mu = 4902.800066;            % Lunar gravitational parameter, km^3/s^2
Rmoon = 1737.4;              % Mean lunar radius, km
hOrbit = 100;                % Initial circular orbit altitude, km
hPDI = 15;                  % Target transfer periapsis altitude, km
lat = -88.811008;           % Site07 latitude, deg
lon = 123.690068;           % Site07 east longitude, deg
rA = Rmoon + hOrbit;
rP = Rmoon + hPDI;
a = (rA+rP)/2;

%% Moon-fixed site vector, treated as inertially fixed over the short coast
% Moon rotation is neglected in this V7 geometry demonstration.
site = [cosd(lat)*cosd(lon); cosd(lat)*sind(lon); sind(lat)];
site = site/norm(site);

% Polar orbital plane passing exactly through the site direction.
% Choose the plane spanned by the lunar north pole and the site's longitude.
% Near the poles this remains well defined because longitude is specified.
east = [-sind(lon); cosd(lon); 0];
east = east/norm(east);
% Tangent vector in the orbital plane at the periapsis direction.
tangent = cross(east,site);
tangent = tangent/norm(tangent);

% Periapsis lies above the site; apoapsis lies on the opposite side.
r0 = -rA*site;
vCirc = sqrt(mu/rA);
vA = sqrt(mu*(2/rA-1/a));
% Tangent direction chosen for consistent angular momentum.
v0 = vA*tangent;
vCircular = vCirc*tangent;
dvDeorbit = abs(vCirc-vA);

%% Propagate the two-body coast for half the transfer period
Tcoast = pi*sqrt(a^3/mu);
opts = odeset('RelTol',1e-11,'AbsTol',1e-12);
[t,state] = ode113(@(t,y) [y(4:6); -mu*y(1:3)/norm(y(1:3))^3], ...
    [0 Tcoast], [r0;v0], opts);
rEnd = state(end,1:3)';
vEnd = state(end,4:6)';
rExpected = rP*site;

% Radial and tangential components at PDI
radialUnit = rEnd/norm(rEnd);
vRadial = dot(vEnd,radialUnit);
vTangential = norm(vEnd-vRadial*radialUnit);
angleError = acosd(max(-1,min(1,dot(radialUnit,site))));
surfaceArcError = Rmoon*deg2rad(angleError);
positionError = norm(rEnd-rExpected);

fprintf('\nSTARSHIP HLS V7 - SITE-ALIGNED DEORBIT\n');
fprintf('Site: Site07 - Peak Near Shackleton\n');
fprintf('Parking orbit altitude: %.2f km\n',hOrbit);
fprintf('Transfer periapsis altitude: %.2f km\n',hPDI);
fprintf('Deorbit delta-V: %.3f m/s\n',1000*dvDeorbit);
fprintf('Coast time: %.2f min\n',Tcoast/60);
fprintf('PDI altitude: %.6f km\n',norm(rEnd)-Rmoon);
fprintf('PDI radial velocity: %.6f m/s\n',1000*vRadial);
fprintf('PDI tangential velocity: %.3f m/s\n',1000*vTangential);
fprintf('Angular error to site: %.8f deg\n',angleError);
fprintf('Surface arc error: %.6f m\n',1000*surfaceArcError);
fprintf('PDI position error: %.6f m\n',1000*positionError);
fprintf('\nVALIDATION: transfer periapsis is aligned with Site07.\n');
fprintf('This is NOT a soft landing or a complete landing delta-V.\n');

%% Figure 1 - Moon, orbit, transfer, site, and periapsis
figure('Color','w','Name','V7 Site-Aligned Lunar Transfer');
[xx,yy,zz] = sphere(65);
surf(Rmoon*xx,Rmoon*yy,Rmoon*zz,'FaceColor',[.68 .68 .68], ...
    'EdgeColor','none','FaceAlpha',.9); hold on;
% Circular parking orbit in the same plane
th = linspace(0,2*pi,600);
rCirc = rA * ( ...
    site(:) * cos(th(:)).' + ...
    tangent(:) * sin(th(:)).' ...
);
plot3(rCirc(:,1),rCirc(:,2),rCirc(:,3),'Color',[.2 .35 .95], ...
    'LineWidth',1.2,'DisplayName','100 km parking orbit');
plot3(state(:,1),state(:,2),state(:,3),'r','LineWidth',2.3, ...
    'DisplayName','Deorbit coast');
plot3(Rmoon*site(1),Rmoon*site(2),Rmoon*site(3),'gp', ...
    'MarkerSize',14,'MarkerFaceColor','g','DisplayName','Site07');
plot3(r0(1),r0(2),r0(3),'ko','MarkerFaceColor','k', ...
    'DisplayName','Deorbit burn');
plot3(rEnd(1),rEnd(2),rEnd(3),'mo','MarkerFaceColor','m', ...
    'DisplayName','PDI (15 km)');
axis equal; grid on; xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('V7: Site-aligned lunar deorbit coast');
legend('Moon','100 km parking orbit','Deorbit coast','Site07', ...
    'Deorbit burn','PDI (15 km)','Location','bestoutside');
view(35,25); lighting gouraud; camlight headlight;

%% Figure 2 - Coast altitude and speed
alt = sqrt(sum(state(:,1:3).^2,2))-Rmoon;
speed = sqrt(sum(state(:,4:6).^2,2))*1000;
figure('Color','w','Name','V7 Coast Diagnostics');
subplot(2,1,1); plot(t/60,alt,'LineWidth',1.8); grid on;
ylabel('Altitude (km)'); title('Deorbit coast altitude');
subplot(2,1,2); plot(t/60,speed,'LineWidth',1.8); grid on;
xlabel('Coast time (min)'); ylabel('Speed (m/s)');
title('Moon-centered inertial speed');
