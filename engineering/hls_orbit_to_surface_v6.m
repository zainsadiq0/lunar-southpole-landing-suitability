clear; clc; close all;

%% User inputs (SI units)
p.mu=4.9048695e12;          % lunar gravitational parameter (m^3/s^2)
p.R=1737.4e3;              % lunar mean radius (m)
p.g0=9.80665;              % standard gravity (m/s^2)
p.Isp=350;                 % illustrative effective Isp (s)
p.Tmax=2e6;               % illustrative maximum thrust (N)
p.mDry=100000;             % minimum spacecraft mass (kg)
p.m0=250000;              % mass before deorbit burn (kg)
p.hOrbit=100e3;            % circular parking orbit altitude (m)
p.hPDI=15e3;              % intended powered-descent start altitude (m)
p.hPeri=0.8e3;            % coast ellipse periapsis altitude (m)
p.dt=0.10;                % powered-descent RK4 step (s)
p.tMax=1600;              % maximum powered-descent duration (s)
p.siteLat=-88.811008;     % Site07 latitude (deg)
p.siteLon=123.690068;     % Site07 longitude (deg)
p.speedLimit=2;           % acceptable touchdown speed (m/s)
p.errorLimit=100;         % acceptable horizontal error (m)
% Controller settings (adjust only after inspecting V6 results)
p.KpHoriz=0.00025;        % horizontal position feedback (1/s^2)
p.KdHoriz=0.055;          % horizontal velocity feedback (1/s)
p.KvUp=0.55;              % vertical velocity feedback (1/s)
p.aRef=0.28;              % vertical speed profile deceleration (m/s^2)
p.vFinalCmd=0.5;          % desired near-ground downward speed (m/s)
p.vUpMax=65;              % maximum commanded descent speed (m/s)

%% Site-centered Moon-fixed frame; treated as inertial (no lunar rotation)
lat=deg2rad(p.siteLat); lon=deg2rad(p.siteLon);
up=[cos(lat)*cos(lon);cos(lat)*sin(lon);sin(lat)];
east=[-sin(lon);cos(lon);0]; north=cross(up,east);
B=[east north up]; rSite=p.R*up;

%% Construct a polar orbit plane through Site07
% The orbital plane is chosen to include the landing site, eliminating
% plane-change cost in this baseline. A fixed-RAAN comparison comes later.
z=[0;0;1]; normal=cross(z,up); normal=normal/norm(normal);
forward=cross(normal,up); forward=forward/norm(forward);
% Orbit passes above the site at true anomaly pi from apoapsis.
% Use a 100 km apoapsis and 0.8 km periapsis coast ellipse.
ra=p.R+p.hOrbit; rp=p.R+p.hPeri;
a=(ra+rp)/2;
vCirc=sqrt(p.mu/ra);
vApo=sqrt(p.mu*(2/ra-1/a));
dVdeorbit=vCirc-vApo;
% At apoapsis opposite site, velocity is in the chosen polar plane.
r0=-ra*up; v0=-vApo*forward;
% Exact impulse propellant via rocket equation
mAfterDeorbit=p.m0*exp(-dVdeorbit/(p.Isp*p.g0));

%% Coast: two-body dynamics until first descending crossing of 15 km
opts=odeset('RelTol',1e-10,'AbsTol',1e-8,'Events',@(t,y) pdiEvent(t,y,p));
[tCoast,yCoast,tEvent,yEvent]=ode45(@(t,y) coastODE(t,y,p), ...
    [0 1.5*pi*sqrt(a^3/p.mu)], [r0;v0], opts);
if isempty(tEvent)
    error('Coast did not reach the PDI altitude on descent. Check orbit settings.');
end
rPDI=yEvent(end,1:3)'; vPDI=yEvent(end,4:6)';
coastTime=tEvent(end);
coastLocal=B'*(rPDI-rSite); coastVel=B'*vPDI;

%% Powered descent: integrate actual coast-end state, without resetting speed
state=[rPDI;vPDI;mAfterDeorbit;0]; % [r(3);v(3);m;powered DV]
N=ceil(p.tMax/p.dt)+2;
Y=nan(N,8); t=nan(N,1); thrust=nan(N,1);
reason='Maximum powered-descent time'; n=0;
for k=1:N
    tk=(k-1)*p.dt;
    r=state(1:3); v=state(4:6); m=state(7);
    h=norm(r)-p.R;
    if h<=0, reason='Surface contact'; break; end
    if m<=p.mDry+1e-6, reason='Propellant exhausted'; break; end
    localPos=B'*(r-rSite); localVel=B'*v;
    % Horizontal PD guidance in the site frame; includes full orbital speed.
    aEast=-p.KpHoriz*localPos(1)-p.KdHoriz*localVel(1);
    aNorth=-p.KpHoriz*localPos(2)-p.KdHoriz*localVel(2);
    vUpCmd=-min(p.vUpMax,max(p.vFinalCmd,sqrt(2*p.aRef*h)));
    aUp=p.KvUp*(vUpCmd-localVel(3));
    g=-p.mu*r/norm(r)^3;
    F=m*(B*[aEast;aNorth;aUp]-g);
    Fmag=norm(F);
    if Fmag>p.Tmax, F=F*(p.Tmax/Fmag); Fmag=p.Tmax; end
    n=n+1; t(n)=tk; Y(n,:)=state'; thrust(n)=Fmag;
    if k==N, break; end
    f=@(s) poweredODE(s,F,p);
    k1=f(state); k2=f(state+p.dt*k1/2);
    k3=f(state+p.dt*k2/2); k4=f(state+p.dt*k3);
    next=state+p.dt*(k1+2*k2+2*k3+k4)/6;
    % Linear event interpolation over the short integration step.
    hNext=norm(next(1:3))-p.R;
    if hNext<=0
        alpha=h/(h-hNext);
        state=state+alpha*(next-state);
        n=n+1; t(n)=tk+alpha*p.dt; Y(n,:)=state'; thrust(n)=Fmag;
        reason='Surface contact'; break;
    end
    if next(7)<=p.mDry
        alpha=(state(7)-p.mDry)/(state(7)-next(7));
        state=state+alpha*(next-state);
        n=n+1; t(n)=tk+alpha*p.dt; Y(n,:)=state'; thrust(n)=Fmag;
        reason='Propellant exhausted'; break;
    end
    state=next;
end
Y=Y(1:n,:);t=t(1:n);thrust=thrust(1:n);
endPos=B'*(Y(end,1:3)'-rSite); endVel=B'*Y(end,4:6)';
positionError=norm(endPos(1:2)); finalSpeed=norm(endVel);
DVpowered=Y(end,8); DVtotal=dVdeorbit+DVpowered;
propDeorbit=p.m0-mAfterDeorbit;
propPowered=mAfterDeorbit-Y(end,7);
propTotal=p.m0-Y(end,7);
success=strcmp(reason,'Surface contact') && ...
    finalSpeed<=p.speedLimit && positionError<=p.errorLimit;

fprintf('\nSTARSHIP HLS V6 - ORBIT TO SURFACE (EDUCATIONAL MODEL)\n');
fprintf('Site: Site07 - Peak Near Shackleton\n');
fprintf('Initial orbit altitude: %.1f km\n',p.hOrbit/1000);
fprintf('Deorbit impulse: %.2f m/s\n',dVdeorbit);
fprintf('Coast duration: %.2f min\n',coastTime/60);
fprintf('PDI altitude: %.3f km\n',(norm(rPDI)-p.R)/1000);
fprintf('PDI east/north/up velocity: %.2f / %.2f / %.2f m/s\n',coastVel);
fprintf('PDI horizontal site offset: %.2f km\n',norm(coastLocal(1:2))/1000);
fprintf('Termination: %s\n',reason);
fprintf('Powered descent duration: %.2f s\n',t(end));
fprintf('Final altitude: %.3f m\n',norm(Y(end,1:3))-p.R);
fprintf('Landing position error: %.2f m\n',positionError);
fprintf('Final east/north/up velocity: %.2f / %.2f / %.2f m/s\n',endVel);
fprintf('Final speed: %.2f m/s\n',finalSpeed);
fprintf('Powered delta-V: %.2f m/s\n',DVpowered);
fprintf('TOTAL delta-V: %.2f m/s\n',DVtotal);
fprintf('Deorbit / powered / total propellant: %.1f / %.1f / %.1f kg\n', ...
    propDeorbit,propPowered,propTotal);
fprintf('Remaining mass: %.1f kg\n',Y(end,7));
fprintf('SUCCESSFUL LANDING: %d\n',success);
if ~success
    fprintf('NOT FEASIBLE: do not interpret delta-V or fuel as landing requirements.\n');
    fprintf('Next task: tune guidance and/or PDI geometry before site comparison.\n');
end

%% Plots
xyz=(B'*(Y(:,1:3)'-rSite))'; vel=(B'*Y(:,4:6)')';
histAlt=vecnorm(Y(:,1:3),2,2)-p.R;
figure('Name','V6 orbital coast and powered descent');
tiledlayout(3,2);
nexttile;plot(t,histAlt/1000);grid on;ylabel('Altitude (km)');xlabel('Powered time (s)');
nexttile;plot(t,vel);grid on;ylabel('Velocity (m/s)');xlabel('Powered time (s)');legend('East','North','Up');
nexttile;plot(t,Y(:,7)/1000);grid on;ylabel('Mass (1000 kg)');xlabel('Powered time (s)');
nexttile;plot(t,dVdeorbit+Y(:,8));grid on;ylabel('Total delta-V (m/s)');xlabel('Powered time (s)');
nexttile;plot(t,thrust/1000);grid on;ylabel('Thrust (kN)');xlabel('Powered time (s)');
nexttile;plot(xyz(:,1)/1000,xyz(:,2)/1000);hold on;
plot(0,0,'rx','MarkerSize',10,'LineWidth',2);axis equal;grid on;
xlabel('East (km)');ylabel('North (km)');title('Ground track');
figure('Name','V6 3D orbit-to-surface');
[xx,yy,zz]=sphere(45);surf(p.R*xx/1000,p.R*yy/1000,p.R*zz/1000, ...
    'FaceColor',[.55 .55 .55],'EdgeColor','none','FaceAlpha',.55);hold on;
plot3(yCoast(:,1)/1000,yCoast(:,2)/1000,yCoast(:,3)/1000,'b','LineWidth',1.4);
plot3(Y(:,1)/1000,Y(:,2)/1000,Y(:,3)/1000,'r','LineWidth',1.5);
plot3(rSite(1)/1000,rSite(2)/1000,rSite(3)/1000,'gx','MarkerSize',10,'LineWidth',2);
axis equal;grid on;view(35,25);xlabel('X (km)');ylabel('Y (km)');zlabel('Z (km)');
legend('Moon','Coast','Powered descent','Site07');title('V6 lunar descent trajectory');

%% Local functions
function dy=coastODE(~,y,p)
 r=y(1:3);dy=[y(4:6);-p.mu*r/norm(r)^3];
end
function [value,isterminal,direction]=pdiEvent(~,y,p)
 value=norm(y(1:3))-(p.R+p.hPDI);isterminal=1;direction=-1;
end
function dy=poweredODE(y,F,p)
 r=y(1:3);m=max(y(7),1);
 dy=[y(4:6);-p.mu*r/norm(r)^3+F/m; ...
     -norm(F)/(p.Isp*p.g0);norm(F)/m];
end
