clear; clc; close all;

%% Parameters (illustrative assumptions)
p.mu = 4.9048695e12;          % lunar gravitational parameter (m^3/s^2)
p.R = 1737.4e3;             % spherical Moon radius (m)
p.g0 = 9.80665;             % standard gravity (m/s^2)
p.Isp = 350;                % assumed effective specific impulse (s)
p.Tmax = 2.0e6;            % assumed maximum total thrust (N)
p.mDry = 100000;           % assumed minimum vehicle mass (kg)
p.m0 = 250000;             % assumed initial vehicle mass (kg)
p.dt = 0.10;               % time step (s)
p.tMax = 1200;             % maximum simulation time (s)
p.siteLat = -88.811008;    % Site07, deg
p.siteLon = 123.690068;    % Site07, deg
p.siteElevation = 0;       % simplified local surface relative to sphere (m)
p.speedLimit = 2;          % landing-speed criterion (m/s)
p.errorLimit = 100;
p.aRef = 0.28;             % reference vertical braking acceleration (m/s^2)
p.vFinalCmd = 0.50;        % minimum target downward speed (m/s)
p.KvUp = 0.55;            % vertical velocity feedback gain (1/s)
p.vUpMax = 65;            % maximum target descent speed (m/s)        % horizontal position tolerance (m)

%% Local frame: east (E), north (N), up (U)
lat = deg2rad(p.siteLat); lon = deg2rad(p.siteLon);
up = [cos(lat)*cos(lon); cos(lat)*sin(lon); sin(lat)];
east = [-sin(lon); cos(lon); 0];
north = cross(up,east);
B = [east north up];       % local-to-Moon-centered inertial rotation
rSite = (p.R+p.siteElevation)*up;

% Assumed start 15 km above site with 3 km east offset.
% These conditions are NOT derived from the prior coast simulation.
r0 = rSite + B*[3000; 0; 15000];
v0 = B*[120; 0; -20];       % initial east/north/up velocity (m/s)
y = [r0; v0; p.m0; 0];     % position, velocity, mass, accumulated DV

N = floor(p.tMax/p.dt)+1;
Y = nan(N,8); time = nan(N,1); thrustHist = nan(N,1);
reason = 'Maximum time reached'; n = 0;

%% Time integration: RK4 with held thrust command each step
for k = 1:N
    t = (k-1)*p.dt;
    r = y(1:3); v = y(4:6); m = y(7);
    localPos = B'*(r-rSite);
    localVel = B'*v;
    altitude = norm(r)-(p.R+p.siteElevation);

    % Contact / failure event checks BEFORE taking another step
    if altitude <= 0
        reason = 'Surface contact';
        break
    end
    if m <= p.mDry+1e-6
        reason = 'Propellant exhausted';
        break
    end

    % Velocity targets depend on remaining distance to the site.
    % Slower targets near the surface reduce landing speed.
    vEastCmd = -max(-60,min(60,localPos(1)/35));
    vNorthCmd = -max(-60,min(60,localPos(2)/35));
    % Braking profile: slow down progressively as altitude decreases.
    % A 0.5 m/s terminal target avoids commanding zero speed before contact.
    vUpCmd = -min(p.vUpMax,max(p.vFinalCmd, ...
                 sqrt(2*p.aRef*max(altitude,0))));
    vCmd = [vEastCmd; vNorthCmd; vUpCmd];

    % Velocity feedback + weak position correction in horizontal axes.
    aLocal = [0.10*(vCmd(1)-localVel(1)); ...
              0.10*(vCmd(2)-localVel(2)); ...
              p.KvUp*(vCmd(3)-localVel(3))];
    aDesired = B*aLocal;
    g = -p.mu*r/norm(r)^3;
    F = m*(aDesired-g);  % gravity compensation
    Fmag = norm(F);
    if Fmag > p.Tmax
        F = F*(p.Tmax/Fmag);
        Fmag = p.Tmax;
    end
    if Fmag < 0, error('Invalid thrust'); end

    n = n+1;
    time(n)=t; Y(n,:)=y'; thrustHist(n)=Fmag;
    if k==N, break; end

    h = p.dt;
    f = @(state) dynamics(state,F,p);
    k1 = f(y);
    k2 = f(y+0.5*h*k1);
    k3 = f(y+0.5*h*k2);
    k4 = f(y+h*k3);
    ynext = y + h*(k1+2*k2+2*k3+k4)/6;

    % Interpolate surface crossing so touchdown metrics are not below ground.
    altNext = norm(ynext(1:3))-(p.R+p.siteElevation);
    if altNext <= 0
        frac = altitude/(altitude-altNext);
        y = y + frac*(ynext-y);
        n=n+1; time(n)=t+frac*h; Y(n,:)=y'; thrustHist(n)=Fmag;
        reason='Surface contact';
        break
    end
    y=ynext;
end

Y=Y(1:n,:); time=time(1:n); thrustHist=thrustHist(1:n);
localEnd = B'*(Y(end,1:3)'-rSite);
velEnd = B'*Y(end,4:6)';
positionError = norm(localEnd(1:2));
touchdownSpeed = norm(velEnd);
usedProp = p.m0-Y(end,7);
DV = Y(end,8);
success = strcmp(reason,'Surface contact') && ...
          touchdownSpeed<=p.speedLimit && positionError<=p.errorLimit;

fprintf('\nSTARSHIP HLS V5 - CONTROLLED DESCENT DEMONSTRATION\n');
fprintf('Site: Site07 - Peak Near Shackleton\n');
fprintf('Termination: %s\n',reason);
fprintf('Duration: %.2f s\n',time(end));
fprintf('Final altitude: %.3f m\n',norm(Y(end,1:3))-p.R-p.siteElevation);
fprintf('Horizontal position error: %.2f m\n',positionError);
fprintf('East/North/Up velocity: %.3f / %.3f / %.3f m/s\n',velEnd);
fprintf('Final speed: %.3f m/s\n',touchdownSpeed);
fprintf('Powered delta-V: %.2f m/s\n',DV);
fprintf('Propellant consumed: %.2f kg\n',usedProp);
fprintf('Remaining mass: %.2f kg\n',Y(end,7));
fprintf('Successful landing: %d\n',success);
if ~success
    fprintf('Landing constraints not met; do not use as a feasible fuel estimate.\n');
end

%% Plots
localXYZ = (B'*(Y(:,1:3)'-rSite))';
localV = (B'*Y(:,4:6)')';
altitude = vecnorm(Y(:,1:3),2,2)-(p.R+p.siteElevation);
figure('Name','HLS V5 descent');
tiledlayout(3,2);
nexttile; plot(time,altitude/1000); grid on; ylabel('Altitude (km)'); xlabel('Time (s)');
nexttile; plot(time,localV); grid on; ylabel('Velocity (m/s)'); xlabel('Time (s)'); legend('East','North','Up');
nexttile; plot(time,Y(:,7)/1000); grid on; ylabel('Mass (1000 kg)'); xlabel('Time (s)');
nexttile; plot(time,Y(:,8)); grid on; ylabel('Powered delta-V (m/s)'); xlabel('Time (s)');
nexttile; plot(time,thrustHist/1000); grid on; ylabel('Thrust (kN)'); xlabel('Time (s)');
nexttile; plot(localXYZ(:,1)/1000,localXYZ(:,2)/1000); hold on;
plot(0,0,'rx','MarkerSize',10,'LineWidth',2); grid on; axis equal;
xlabel('East (km)'); ylabel('North (km)'); title('Ground track');
figure('Name','HLS V5 local 3D trajectory');
plot3(localXYZ(:,1)/1000,localXYZ(:,2)/1000,localXYZ(:,3)/1000,'LineWidth',1.5); hold on;
plot3(0,0,0,'rx','MarkerSize',10,'LineWidth',2);
grid on; axis equal; xlabel('East (km)'); ylabel('North (km)'); zlabel('Up (km)');
title('Site07 powered descent - local frame');

%% Local dynamics function
function dydt = dynamics(y,F,p)
    r=y(1:3); v=y(4:6); m=max(y(7),1);
    g=-p.mu*r/norm(r)^3;
    mdot=-norm(F)/(p.Isp*p.g0);
    dVdot=norm(F)/m;
    dydt=[v; g+F/m; mdot; dVdot];
end
