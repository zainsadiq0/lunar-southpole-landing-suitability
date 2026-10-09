clear; clc; close all;

%% Assumptions (SI units)
mu=4.902800066e12; R=1737.4e3;
hOrbit=100e3; hPeri=15e3;
lat=-88.811008; lon=123.690068; % Site07
m0=250000; mDry=120000; Tmax=2e6; Isp=350; g0=9.80665;
site=[cosd(lat)*cosd(lon);cosd(lat)*sind(lon);sind(lat)]; site=site/norm(site);
east=[-sind(lon);cosd(lon);0]; east=east/norm(east);
tangent=cross(east,site); tangent=tangent/norm(tangent);
rA=R+hOrbit; rP=R+hPeri; a=(rA+rP)/2;
vCirc=sqrt(mu/rA); vA=sqrt(mu*(2/rA-1/a));
dvDeorbit=vCirc-vA;
mStart=m0*exp(-dvDeorbit/(Isp*g0));
Tcoast=pi*sqrt(a^3/mu);
y0=[-rA*site;vA*tangent];
opts=odeset('RelTol',1e-9,'AbsTol',1e-7);
[tCoast,yCoast]=ode45(@(t,y) [y(4:6);-mu*y(1:3)/norm(y(1:3))^3],...
    linspace(0,Tcoast,1500),y0,opts);

%% Test braking starts from the second half of transfer
fractions=0.50:0.025:0.975;
N=numel(fractions);
alt0=zeros(N,1); range0=zeros(N,1); speed0=zeros(N,1);
stopDist=zeros(N,1); stopTime=zeros(N,1); propEstimate=zeros(N,1);
finalAlt=zeros(N,1); finalRange=zeros(N,1); finalSpeed=zeros(N,1);
used=zeros(N,1); burnTime=zeros(N,1); termination=strings(N,1);

for k=1:N
    state=interp1(tCoast,yCoast,fractions(k)*Tcoast).';
    r=state(1:3); v=state(4:6);
    alt0(k)=(norm(r)-R)/1000;
    range0(k)=siteRange(r,site,R)/1000;
    speed0(k)=norm(v);
    % Constant-acceleration diagnostic only, not a trajectory solution.
    g=mu/norm(r)^2;
    aConservative=max(Tmax/mStart-g,1e-6);
    stopTime(k)=speed0(k)/aConservative;
    stopDist(k)=speed0(k)^2/(2*aConservative)/1000;
    propEstimate(k)=mStart*(1-exp(-speed0(k)/(Isp*g0)));

    % Numerical 3D propagation: thrust always opposite inertial velocity.
    % Stops at surface, dry mass, or low speed. This is NOT site guidance.
    initial=[state;mStart;0];
    optBurn=odeset('RelTol',1e-7,'AbsTol',1e-5,...
        'Events',@(t,y) stopEvents(y,R,mDry));
    [tb,yb,~,~,ie]=ode45(@(t,y) retroDynamics(y,mu,Tmax,Isp,g0,mDry),...
        [0 1800],initial,optBurn);
    last=yb(end,:).';
    finalAlt(k)=(norm(last(1:3))-R)/1000;
    finalRange(k)=siteRange(last(1:3),site,R)/1000;
    finalSpeed(k)=norm(last(4:6));
    used(k)=mStart-last(7); burnTime(k)=tb(end);
    if any(ie==1)
        termination(k)="Surface impact";
    elseif any(ie==2)
        termination(k)="Propellant limit";
    elseif any(ie==3)
        termination(k)="Near-zero speed";
    else
        termination(k)="Time limit";
    end
end

results=table(fractions.',alt0,range0,speed0,stopDist,stopTime,...
    propEstimate,finalAlt,finalRange,finalSpeed,used,burnTime,termination,...
    'VariableNames',{'CoastFraction','StartAlt_km','SiteRange_km',...
    'StartSpeed_mps','IdealStopDist_km','IdealStopTime_s',...
    'IdealPropellant_kg','EndAlt_km','EndSiteRange_km','EndSpeed_mps',...
    'BurnPropellant_kg','BurnTime_s','Termination'});

fprintf('\nSTARSHIP HLS V10 - BRAKING FEASIBILITY STUDY\n');
fprintf('Site07 - Peak Near Shackleton\n');
fprintf('Parking orbit: %.1f km | Transfer periapsis: %.1f km\n',hOrbit/1000,hPeri/1000);
fprintf('Deorbit delta-V: %.3f m/s | Coast to periapsis: %.2f min\n',dvDeorbit,Tcoast/60);
fprintf('Initial post-deorbit mass: %.1f kg | Dry mass: %.1f kg\n',mStart,mDry);
fprintf('Thrust: %.0f kN | Isp: %.0f s\n\n',Tmax/1000,Isp);
disp(results);
fprintf('\nNOTE: Ideal stopping distance/time ignore orbital curvature, mass change,\n');
fprintf('and thrust direction constraints. Ideal propellant assumes an impulsive\n');
fprintf('rocket-equation velocity change. Neither is a landing fuel estimate.\n');
fprintf('Retrograde trials are diagnostic only: they do NOT steer to Site07.\n');
fprintf('A near-zero speed above the surface is NOT a successful landing.\n');

%% Dark diagnostic plots
bg=[0.06 0.06 0.06]; fg=[0.88 0.88 0.88];
figure('Color',bg,'Name','V10 Braking Feasibility');
subplot(2,2,1);
plot(fractions,range0,'o-','LineWidth',1.5); hold on;
plot(fractions,stopDist,'s-','LineWidth',1.5);
xlabel('Coast fraction'); ylabel('Distance (km)');
title('Site range and ideal braking distance'); legend('Range to Site07','Ideal stop distance');
subplot(2,2,2);
plot(fractions,alt0,'o-','LineWidth',1.5); hold on;
plot(fractions,finalAlt,'s-','LineWidth',1.5);
xlabel('Coast fraction'); ylabel('Altitude (km)');
title('Start and termination altitude'); legend('Start','End');
subplot(2,2,3);
plot(fractions,finalSpeed,'o-','LineWidth',1.5);
xlabel('Coast fraction'); ylabel('Speed (m/s)'); title('Speed at trial termination');
subplot(2,2,4);
plot(fractions,used/1000,'o-','LineWidth',1.5); hold on;
plot(fractions,propEstimate/1000,'s-','LineWidth',1.5);
xlabel('Coast fraction'); ylabel('Propellant (1000 kg)');
title('Retrograde burn vs idealized estimate');
legend('Numerical burn','Ideal rocket equation');
styleDark(gcf,bg,fg);

%% Lunar transfer geometry
figure('Color',bg,'Name','V10 Lunar Transfer Geometry');
[xs,ys,zs]=sphere(60);
surf(R*xs/1000,R*ys/1000,R*zs/1000,'FaceColor',[.45 .45 .45],...
    'EdgeColor','none','FaceAlpha',.85); hold on;
plot3(yCoast(:,1)/1000,yCoast(:,2)/1000,yCoast(:,3)/1000,...
    'Color',[.15 .5 1],'LineWidth',1.8);
plot3(R*site(1)/1000,R*site(2)/1000,R*site(3)/1000,...
    'gp','MarkerFaceColor','g','MarkerSize',13);
for k=1:2:N
    s=interp1(tCoast,yCoast,fractions(k)*Tcoast);
    plot3(s(1)/1000,s(2)/1000,s(3)/1000,'mo','MarkerSize',4);
end
axis equal; grid on; view(35,25);
xlabel('X (km)');ylabel('Y (km)');zlabel('Z (km)');
title('V10: Site-aligned transfer and candidate braking starts');
legend('Moon','Transfer','Site07','Candidate starts','Location','bestoutside');
styleDark(gcf,bg,fg);

%% Local functions
function dydt=retroDynamics(y,mu,Tmax,Isp,g0,mDry)
    r=y(1:3);v=y(4:6);m=y(7);
    grav=-mu*r/norm(r)^3;
    if m<=mDry || norm(v)<0.5
        T=0;aThrust=zeros(3,1);
    else
        T=Tmax;aThrust=-(T/m)*v/norm(v);
    end
    dydt=[v;grav+aThrust;-T/(Isp*g0);norm(aThrust)];
end
function [value,isterminal,direction]=stopEvents(y,R,mDry)
    value=[norm(y(1:3))-R;y(7)-mDry;norm(y(4:6))-0.5];
    isterminal=[1;1;1];direction=[-1;-1;-1];
end
function d=siteRange(r,site,R)
    u=r/norm(r);
    d=R*acos(max(-1,min(1,dot(u,site))));
end
function styleDark(fig,bg,fg)
    set(fig,'Color',bg);
    axesList=findall(fig,'Type','axes');
    for j=1:numel(axesList)
        ax=axesList(j);
        set(ax,'Color',[.085 .085 .085],'XColor',fg,'YColor',fg,...
            'ZColor',fg,'GridColor',[.48 .48 .48],...
            'GridAlpha',.35,'FontSize',10);
        grid(ax,'on');
        set(get(ax,'Title'),'Color',fg);
        set(get(ax,'XLabel'),'Color',fg);
        set(get(ax,'YLabel'),'Color',fg);
        set(get(ax,'ZLabel'),'Color',fg);
    end
    legends=findall(fig,'Type','legend');
    for j=1:numel(legends)
        set(legends(j),'Color',[.10 .10 .10],'TextColor',fg,...
            'EdgeColor',[.55 .55 .55]);
    end
end
