clear; clc; close all;

%% Constants and assumptions
mu_km = 4902.800066;       % Moon GM, km^3/s^2
mu = mu_km*1e9;            % m^3/s^2
R = 1737.4e3;             % Moon radius, m
hOrbit = 100e3;           % Parking orbit altitude, m
hPeri = 15e3;             % Transfer periapsis altitude, m
lat = -88.811008; lon = 123.690068;
m0 = 250000;              % Initial mass at deorbit, kg (assumed)
mDry = 120000;            % Minimum allowable mass, kg (assumed)
Tmax = 2.0e6;             % Maximum thrust, N (assumed)
Isp = 350;                % Specific impulse, s (assumed)
g0 = 9.80665;             % Standard gravity, m/s^2
speedLimit = 2;           % Touchdown speed limit, m/s
posLimit = 100;           % Target position tolerance, m

%% Site-aligned transfer (same geometry as V7)
site = [cosd(lat)*cosd(lon); cosd(lat)*sind(lon); sind(lat)];
site = site/norm(site);
east = [-sind(lon); cosd(lon); 0]; east = east/norm(east);
tangent = cross(east,site); tangent = tangent/norm(tangent);
rA = R+hOrbit; rP = R+hPeri; a = (rA+rP)/2;
vCirc = sqrt(mu/rA);
vA = sqrt(mu*(2/rA-1/a));
r0 = -rA*site; v0 = vA*tangent;
dvDeorbit = vCirc-vA;
Tcoast = pi*sqrt(a^3/mu);
optsCoast = odeset('RelTol',1e-10,'AbsTol',1e-8);
[tCoast,yCoast] = ode45(@(t,y) [y(4:6); -mu*y(1:3)/norm(y(1:3))^3], ...
    linspace(0,Tcoast,1500),[r0;v0],optsCoast);

%% Start braking before periapsis: selectable coast fraction
% 0.65 means powered descent starts 65% of the way to periapsis.
% Adjust this parameter to explore different PDI points.
coastFraction = 0.65;
tPDI = coastFraction*Tcoast;
yPDI = interp1(tCoast,yCoast,tPDI).';
rPDI = yPDI(1:3); vPDI = yPDI(4:6);

%% Powered descent using a bounded PD acceleration controller
% Target is fixed in an inertial frame for this short demonstration.
% Lunar rotation, terrain, attitude, engine restart and throttle minimum
% are not modeled. Controller gains are illustrative, not optimized.
Kp = 1.5e-5;            % Position gain, 1/s^2
Kd = 0.011;             % Velocity gain, 1/s
maxPoweredTime = 2200; % s
initialMass = m0*exp(-dvDeorbit/(Isp*g0));
yInit = [rPDI;vPDI;initialMass;0]; % Last state: accumulated powered dV
optsPower = odeset('RelTol',2e-7,'AbsTol',1e-5, ...
    'Events',@(t,y) stopEvents(t,y,R,mDry));
[tPower,yPower,te,ye,ie] = ode45(@(t,y) poweredODE(y,mu,R,site, ...
    Kp,Kd,Tmax,Isp,g0,mDry),[0 maxPoweredTime],yInit,optsPower);

%% Metrics (do not declare success unless contact and constraints pass)
rFinal = yPower(end,1:3).';
vFinal = yPower(end,4:6).';
mFinal = yPower(end,7);
hFinal = norm(rFinal)-R;
target = R*site;
angularError = acos(max(-1,min(1,dot(rFinal/norm(rFinal),site))));
siteError = R*angularError;
finalSpeed = norm(vFinal);
contact = any(ie==1);
fuelExhausted = any(ie==2);
success = contact && finalSpeed<=speedLimit && siteError<=posLimit ...
    && mFinal>=mDry-1e-3;
dvPowered = yPower(end,8);
propDeorbit = m0-initialMass;
propPowered = initialMass-mFinal;

fprintf('\nSTARSHIP HLS V8 - INTEGRATED DESCENT\n');
fprintf('Site: Site07 - Peak Near Shackleton\n');
fprintf('PDI coast fraction: %.2f\n',coastFraction);
fprintf('Coast duration: %.2f min\n',tPDI/60);
fprintf('PDI altitude: %.2f km\n',(norm(rPDI)-R)/1000);
fprintf('PDI speed: %.2f m/s\n',norm(vPDI));
fprintf('PDI ground-range to target: %.2f km\n', ...
    R*acos(max(-1,min(1,dot(rPDI/norm(rPDI),site))))/1000);
fprintf('Powered descent time: %.2f s\n',tPower(end));
fprintf('Final altitude: %.2f m\n',hFinal);
fprintf('Landing position error: %.2f m\n',siteError);
fprintf('Final speed: %.2f m/s\n',finalSpeed);
fprintf('Deorbit / powered / total dV: %.2f / %.2f / %.2f m/s\n', ...
    dvDeorbit,dvPowered,dvDeorbit+dvPowered);
fprintf('Deorbit / powered / total propellant: %.1f / %.1f / %.1f kg\n', ...
    propDeorbit,propPowered,propDeorbit+propPowered);
fprintf('Remaining mass: %.1f kg\n',mFinal);
fprintf('Surface contact: %d | Fuel limit reached: %d\n',contact,fuelExhausted);
fprintf('SUCCESSFUL LANDING: %d\n',success);
if ~success
    fprintf('NOT FEASIBLE: Do not interpret dV or fuel as a landing requirement.\n');
end

%% Diagnostics
alt = (sqrt(sum(yPower(:,1:3).^2,2))-R)/1000;
spd = sqrt(sum(yPower(:,4:6).^2,2));
range = zeros(size(tPower));
for k=1:numel(tPower)
    u = yPower(k,1:3).'; u=u/norm(u);
    range(k)=R*acos(max(-1,min(1,dot(u,site))))/1000;
end
figure('Color','w','Name','V8 Descent Diagnostics');
subplot(2,2,1); plot(tPower,alt,'LineWidth',1.5); grid on;
xlabel('Powered time (s)'); ylabel('Altitude (km)');
subplot(2,2,2); plot(tPower,spd,'LineWidth',1.5); grid on;
xlabel('Powered time (s)'); ylabel('Speed (m/s)');
subplot(2,2,3); plot(tPower,range,'LineWidth',1.5); grid on;
xlabel('Powered time (s)'); ylabel('Site range (km)');
subplot(2,2,4); plot(tPower,yPower(:,7)/1000,'LineWidth',1.5); grid on;
xlabel('Powered time (s)'); ylabel('Mass (1000 kg)');

figure('Color','w','Name','V8 Orbit and Descent');
[xx,yy,zz]=sphere(45);
surf(R*xx/1000,R*yy/1000,R*zz/1000,'FaceColor',[.65 .65 .65], ...
    'EdgeColor','none','FaceAlpha',0.75); hold on;
plot3(yCoast(:,1)/1000,yCoast(:,2)/1000,yCoast(:,3)/1000, ...
    'b','LineWidth',1.5);
plot3(yPower(:,1)/1000,yPower(:,2)/1000,yPower(:,3)/1000, ...
    'r','LineWidth',2);
plot3(target(1)/1000,target(2)/1000,target(3)/1000, ...
    'gp','MarkerFaceColor','g','MarkerSize',13);
plot3(rFinal(1)/1000,rFinal(2)/1000,rFinal(3)/1000, ...
    'mo','MarkerFaceColor','m','MarkerSize',8);
axis equal; grid on; view(35,25);
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('V8: Site-aligned transfer and powered descent');
legend('Moon','Transfer coast','Powered descent','Site07','Final point', ...
    'Location','bestoutside');

%% Local functions
function dydt = poweredODE(y,mu,R,site,Kp,Kd,Tmax,Isp,g0,mDry)
    r=y(1:3); v=y(4:6); m=y(7);
    gravity=-mu*r/norm(r)^3;
    target=R*site;
    % Acceleration command: gravity compensation + PD position/velocity
    aCmd=-Kp*(r-target)-Kd*v-gravity;
    thrustVec=m*aCmd;
    thrustNorm=norm(thrustVec);
    if thrustNorm>Tmax
        thrustVec=thrustVec*(Tmax/thrustNorm);
        thrustNorm=Tmax;
    end
    if m<=mDry
        thrustVec=zeros(3,1); thrustNorm=0;
    end
    dydt=[v;gravity+thrustVec/m;-thrustNorm/(Isp*g0);thrustNorm/m];
end

function [value,isterminal,direction] = stopEvents(~,y,R,mDry)
    value=[norm(y(1:3))-R; y(7)-mDry];
    isterminal=[1;1];
    direction=[-1;-1];
end
