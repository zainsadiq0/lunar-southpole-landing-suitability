clear; clc; close all;

%% Assumed vehicle and lunar parameters
mu = 4.902800066e12;       % m^3/s^2
R = 1737.4e3;              % m
hOrbit = 100e3; hPeri = 15e3;
lat = -88.811008; lon = 123.690068; % Site07 (deg)
m0 = 250000; mDry = 120000; % kg, illustrative assumptions
Tmax = 2.0e6;              % N, assumed thrust limit
Isp = 350; g0 = 9.80665;   % s, m/s^2
speedLimit = 2; posLimit = 100; % m/s, m

%% Site-aligned deorbit (same two-body setup as V7)
site = [cosd(lat)*cosd(lon);cosd(lat)*sind(lon);sind(lat)];
site = site/norm(site);
east = [-sind(lon);cosd(lon);0]; east=east/norm(east);
tangent = cross(east,site); tangent=tangent/norm(tangent);
rA=R+hOrbit; rP=R+hPeri; a=(rA+rP)/2;
vCirc=sqrt(mu/rA); vA=sqrt(mu*(2/rA-1/a));
dvDeorbit=vCirc-vA; mStart=m0*exp(-dvDeorbit/(Isp*g0));
Tcoast=pi*sqrt(a^3/mu);
y0=[-rA*site;vA*tangent];
optsCoast=odeset('RelTol',1e-10,'AbsTol',1e-7);
[tCoast,yCoast]=ode45(@(t,y) [y(4:6);-mu*y(1:3)/norm(y(1:3))^3], ...
    linspace(0,Tcoast,1600),y0,optsCoast);

%% Choose earlier powered-descent initiation; not an optimized value
coastFraction=0.50;  % Try 0.35 to 0.75 in later sensitivity studies
startTime=coastFraction*Tcoast;
yPDI=interp1(tCoast,yCoast,startTime).';

%% Integrate powered descent with thrust and dry-mass limits
maxPoweredTime=1800;
yInit=[yPDI;mStart;0]; % position(3), velocity(3), mass, integrated dV
opts=odeset('RelTol',1e-7,'AbsTol',1e-5, ...
    'Events',@(t,y) endEvents(y,R,mDry));
[t,y,te,ye,ie]=ode45(@(t,y) dynamics(y,mu,R,site,Tmax,Isp,g0,mDry), ...
    [0 maxPoweredTime],yInit,opts);

%% Evaluate contact and errors (never call an impact a landing)
rF=y(end,1:3).'; vF=y(end,4:6).'; mF=y(end,7);
uF=rF/norm(rF); upSpeed=dot(vF,uF);
vHoriz=vF-upSpeed*uF;
angle=acos(max(-1,min(1,dot(uF,site))));
errorM=R*angle;
contact=any(ie==1); fuelLimit=any(ie==2);
success=contact && norm(vF)<=speedLimit && errorM<=posLimit && ~fuelLimit;
propDeorbit=m0-mStart; propPowered=mStart-mF;

fprintf('\nSTARSHIP HLS V9 - DISTANCE-AWARE DESCENT\n');
fprintf('Site: Site07 - Peak Near Shackleton\n');
fprintf('PDI coast fraction: %.2f\n',coastFraction);
fprintf('Coast duration: %.2f min\n',startTime/60);
fprintf('PDI altitude: %.2f km\n',(norm(yPDI(1:3))-R)/1000);
fprintf('PDI speed: %.2f m/s\n',norm(yPDI(4:6)));
fprintf('PDI site ground range: %.2f km\n',siteRange(yPDI(1:3),site,R)/1000);
fprintf('Powered time: %.2f s\n',t(end));
fprintf('Final altitude: %.2f m\n',norm(rF)-R);
fprintf('Site position error: %.2f m\n',errorM);
fprintf('Horizontal / vertical touchdown speed: %.2f / %.2f m/s\n',norm(vHoriz),upSpeed);
fprintf('Final speed: %.2f m/s\n',norm(vF));
fprintf('Deorbit / powered / total dV: %.2f / %.2f / %.2f m/s\n', ...
    dvDeorbit,y(end,8),dvDeorbit+y(end,8));
fprintf('Deorbit / powered / total propellant: %.1f / %.1f / %.1f kg\n', ...
    propDeorbit,propPowered,propDeorbit+propPowered);
fprintf('Remaining mass: %.1f kg\n',mF);
fprintf('Surface contact: %d | Fuel limit: %d | Success: %d\n',contact,fuelLimit,success);
if ~success
    fprintf('NOT FEASIBLE: reported fuel and dV are NOT landing requirements.\n');
end

%% Diagnostic histories
N=numel(t); alt=zeros(N,1); dist=zeros(N,1);
vH=zeros(N,1); vU=zeros(N,1); thrust=zeros(N,1);
for k=1:N
    r=y(k,1:3).'; v=y(k,4:6).'; u=r/norm(r);
    alt(k)=(norm(r)-R)/1000;
    dist(k)=siteRange(r,site,R)/1000;
    vU(k)=dot(v,u); vH(k)=norm(v-vU(k)*u);
    thrust(k)=norm(controlAccel(y(k,:).',mu,R,site,Tmax,mDry)*y(k,7))/1000;
end

% Dark plots: all figures, axes, titles, grids and legends
bg=[0.06 0.06 0.06]; fg=[0.88 0.88 0.88];
figure('Color',bg,'Name','V9 Descent Diagnostics');
subplot(3,2,1); plot(t,alt,'LineWidth',1.6); title('Altitude'); ylabel('km'); xlabel('Time (s)');
subplot(3,2,2); plot(t,vH,'LineWidth',1.6); hold on; plot(t,vU,'LineWidth',1.6); ...
    title('Velocity components'); ylabel('m/s'); xlabel('Time (s)'); legend('Horizontal','Radial');
subplot(3,2,3); plot(t,dist,'LineWidth',1.6); title('Surface distance to target'); ylabel('km'); xlabel('Time (s)');
subplot(3,2,4); plot(t,y(:,7)/1000,'LineWidth',1.6); title('Vehicle mass'); ylabel('1000 kg'); xlabel('Time (s)');
subplot(3,2,5); plot(t,thrust,'LineWidth',1.6); title('Commanded thrust'); ylabel('kN'); xlabel('Time (s)');
subplot(3,2,6); plot(t,y(:,8),'LineWidth',1.6); title('Accumulated powered delta-V'); ylabel('m/s'); xlabel('Time (s)');
styleDark(gcf,bg,fg);

figure('Color',bg,'Name','V9 Lunar Transfer and Descent');
[xs,ys,zs]=sphere(60);
surf(R*xs/1000,R*ys/1000,R*zs/1000,'FaceColor',[0.45 0.45 0.45], ...
    'EdgeColor','none','FaceAlpha',0.82); hold on;
coastMask=tCoast<=startTime;
plot3(yCoast(coastMask,1)/1000,yCoast(coastMask,2)/1000,yCoast(coastMask,3)/1000, ...
    'Color',[0.15 0.45 1],'LineWidth',1.8);
plot3(y(:,1)/1000,y(:,2)/1000,y(:,3)/1000, ...
    'Color',[1 0.25 0.25],'LineWidth',2);
plot3(R*site(1)/1000,R*site(2)/1000,R*site(3)/1000, ...
    'gp','MarkerFaceColor','g','MarkerSize',14);
plot3(rF(1)/1000,rF(2)/1000,rF(3)/1000, ...
    'mo','MarkerFaceColor','m','MarkerSize',8);
axis equal; grid on; view(35,25); xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('V9: Lunar transfer and distance-aware powered descent');
legend('Moon','Coast','Powered descent','Site07','Final point','Location','bestoutside');
styleDark(gcf,bg,fg);

%% Local functions
function dydt=dynamics(y,mu,R,site,Tmax,Isp,g0,mDry)
    r=y(1:3); v=y(4:6); m=y(7);
    grav=-mu*r/norm(r)^3;
    aThrust=controlAccel(y,mu,R,site,Tmax,mDry);
    T=m*norm(aThrust);
    dydt=[v;grav+aThrust;-T/(Isp*g0);norm(aThrust)];
end

function aThrust=controlAccel(y,mu,R,site,Tmax,mDry)
    r=y(1:3); v=y(4:6); m=y(7);
    if m<=mDry
        aThrust=zeros(3,1); return;
    end
    u=r/norm(r); h=max(norm(r)-R,0);
    % Tangent direction toward the target along the lunar surface
    projected=site-dot(site,u)*u;
    pnorm=norm(projected);
    if pnorm>1e-10
        toward=projected/pnorm;
    else
        toward=zeros(3,1);
    end
    vRad=dot(v,u);
    vTan=v-vRad*u;
    range=siteRange(r,site,R);
    % Horizontal speed command reduces with remaining surface range.
    % Limits prevent the guidance from demanding high approach speed near target.
    desiredTan=min(500,sqrt(2*1.5*max(range-100,0)))*toward;
    aTan=0.018*(desiredTan-vTan);
    % Vertical descent schedule: slow down progressively near the ground.
    desiredRad=-min(65,max(0.6,0.08*sqrt(max(h,0))));
    aRad=0.12*(desiredRad-vRad);
    % Gravity compensation and acceleration request in inertial coordinates
    gravity=-mu*r/norm(r)^3;
    aRequested=aTan+aRad*u-gravity;
    aMax=Tmax/m;
    if norm(aRequested)>aMax
        aThrust=aRequested*aMax/norm(aRequested);
    else
        aThrust=aRequested;
    end
end

function d=siteRange(r,site,R)
    u=r/norm(r);
    d=R*acos(max(-1,min(1,dot(u,site))));
end

function [value,isterminal,direction]=endEvents(y,R,mDry)
    value=[norm(y(1:3))-R;y(7)-mDry];
    isterminal=[1;1]; direction=[-1;-1];
end

function styleDark(fig,bg,fg)
    set(fig,'Color',bg);
    axesList=findall(fig,'Type','axes');
    for j=1:numel(axesList)
        ax=axesList(j);
        set(ax,'Color',[0.085 0.085 0.085],'XColor',fg,'YColor',fg, ...
            'ZColor',fg,'GridColor',[0.48 0.48 0.48], ...
            'GridAlpha',0.35,'FontSize',10);
        grid(ax,'on');
        set(get(ax,'Title'),'Color',fg);
        set(get(ax,'XLabel'),'Color',fg);
        set(get(ax,'YLabel'),'Color',fg);
        set(get(ax,'ZLabel'),'Color',fg);
    end
    legends=findall(fig,'Type','legend');
    for j=1:numel(legends)
        set(legends(j),'Color',[0.10 0.10 0.10],'TextColor',fg, ...
            'EdgeColor',[0.55 0.55 0.55]);
    end
end
