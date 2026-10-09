clear; clc; close all;
% Uses V10's two-body transfer, site coordinates and vehicle assumptions.
% Site is fixed in the Moon-centered inertial frame (lunar rotation omitted).
mu=4.902800066e12; R=1737.4e3;
hOrbit=100e3; hPeri=15e3;
lat=-88.811008; lon=123.690068;
m0=250000; mDry=120000; Tmax=2e6; Isp=350; g0=9.80665;
site=[cosd(lat)*cosd(lon);cosd(lat)*sind(lon);sind(lat)]; site=site/norm(site);
east=[-sind(lon);cosd(lon);0]; east=east/norm(east);
tangent=cross(east,site); tangent=tangent/norm(tangent);
rA=R+hOrbit; rP=R+hPeri; a=(rA+rP)/2;
vCirc=sqrt(mu/rA); vA=sqrt(mu*(2/rA-1/a));
dvDeorbit=vCirc-vA; mStart=m0*exp(-dvDeorbit/(Isp*g0));
Tcoast=pi*sqrt(a^3/mu);
y0=[-rA*site;vA*tangent];
opts=odeset('RelTol',1e-9,'AbsTol',1e-7);
[tCoast,yCoast]=ode45(@(t,y) [y(4:6);-mu*y(1:3)/norm(y(1:3))^3],...
    linspace(0,Tcoast,1500),y0,opts);

%% Sweep ignition timing and guidance time-to-go
fractions=0.85:0.025:0.975;
tgoCandidates=[250 400 600 850];
% Position/velocity feedback with thrust limits; not an optimized solution.
maxBurn=2400; speedLimit=2; errorLimit=100;
rows=[]; trials=struct('score',{},'t',{},'y',{},'fraction',{},'tgo',{},'success',{},'termination',{});
for i=1:numel(fractions)
    state=interp1(tCoast,yCoast,fractions(i)*Tcoast).';
    for j=1:numel(tgoCandidates)
        tgo=tgoCandidates(j);
        initial=[state;mStart;0];
        op=odeset('RelTol',2e-7,'AbsTol',1e-5,...
            'Events',@(t,y) descentEvents(y,R,mDry));
        [tb,yb,~,~,ie]=ode45(@(t,y) guidedDynamics(t,y,mu,R,site,Tmax,Isp,g0,mDry,tgo),...
            [0 maxBurn],initial,op);
        last=yb(end,:).';
        alt=(norm(last(1:3))-R);
        err=siteRange(last(1:3),site,R);
        speed=norm(last(4:6));
        contact=any(ie==1); fuel=any(ie==2);
        success=contact && ~fuel && speed<=speedLimit && err<=errorLimit;
        if contact, reason="Surface contact";
        elseif fuel, reason="Dry mass";
        else, reason="Time limit"; end
        prop=mStart-last(7);
        rows(end+1,:)=[fractions(i),tgo,tb(end),alt,err,speed,prop,last(8),double(success)]; 
        % Score penalizes altitude, site miss, speed, and fuel exhaustion.
        score=err/1000+speed/10+max(alt,0)/1000+1e4*double(fuel);
        if success, score=-1e5+prop/1000; end
        trials(end+1)=struct('score',score,'t',tb,'y',yb,'fraction',fractions(i),...
            'tgo',tgo,'success',success,'termination',reason); 
    end
end
results=array2table(rows,'VariableNames',{'CoastFraction','GuidanceTime_s','BurnTime_s',...
    'EndAltitude_m','SiteError_m','EndSpeed_mps','Propellant_kg','PoweredDeltaV_mps','Success'});
disp('STARSHIP HLS V11 - INTEGRATED DESCENT GUIDANCE STUDY');
fprintf('Deorbit delta-V: %.3f m/s | Coast to periapsis: %.2f min\n',dvDeorbit,Tcoast/60);
fprintf('Touchdown requirements: speed <= %.1f m/s; site error <= %.0f m\n',speedLimit,errorLimit);
disp(results);
[~,bestIdx]=min([trials.score]); best=trials(bestIdx);
fprintf('\nBEST TRIAL: coast fraction %.3f | guidance time %.0f s\n',best.fraction,best.tgo);
fprintf('Termination: %s | Successful landing: %d\n',best.termination,best.success);
fprintf('End altitude %.2f m | site error %.2f km | speed %.2f m/s\n',...
    norm(best.y(end,1:3))-R,siteRange(best.y(end,1:3).',site,R)/1000,norm(best.y(end,4:6)));
fprintf('Powered propellant %.1f kg | powered delta-V %.1f m/s\n',...
    mStart-best.y(end,7),best.y(end,8));
if ~any(results.Success)
    fprintf('NO FEASIBLE LANDING FOUND. These are diagnostic trials, NOT landing requirements.\n');
end

%% Dark plots: best candidate and sweep
bg=[.06 .06 .06]; fg=[.88 .88 .88];
t=best.t; y=best.y; n=size(y,1);
alt=(sqrt(sum(y(:,1:3).^2,2))-R)/1000;
speed=sqrt(sum(y(:,4:6).^2,2));
range=zeros(n,1); thrust=zeros(n,1);
for k=1:n
    range(k)=siteRange(y(k,1:3).',site,R)/1000;
    [~,thrust(k)]=guidedDynamics(t(k),y(k,:).',mu,R,site,Tmax,Isp,g0,mDry,best.tgo);
end
figure('Color',bg,'Name','V11 Best Candidate Diagnostics');
subplot(3,2,1); plot(t,alt,'LineWidth',1.7); xlabel('Powered time (s)');ylabel('Altitude (km)');title('Altitude');
subplot(3,2,2); plot(t,speed,'LineWidth',1.7); xlabel('Powered time (s)');ylabel('Speed (m/s)');title('Inertial speed');
subplot(3,2,3); plot(t,range,'LineWidth',1.7); xlabel('Powered time (s)');ylabel('Site range (km)');title('Surface distance to Site07');
subplot(3,2,4); plot(t,y(:,7)/1000,'LineWidth',1.7); xlabel('Powered time (s)');ylabel('Mass (1000 kg)');title('Vehicle mass');
subplot(3,2,5); plot(t,thrust/1000,'LineWidth',1.7); xlabel('Powered time (s)');ylabel('Thrust (kN)');title('Commanded thrust');
subplot(3,2,6); plot(t,y(:,8),'LineWidth',1.7); xlabel('Powered time (s)');ylabel('Delta-V (m/s)');title('Accumulated powered delta-V');
styleDark(gcf,bg,fg);
figure('Color',bg,'Name','V11 Guidance Sweep');
subplot(1,2,1); scatter(results.CoastFraction,results.SiteError_m/1000,50,results.GuidanceTime_s,'filled');
xlabel('Coast fraction');ylabel('Final site error (km)');title('Landing accuracy');colorbar;
subplot(1,2,2); scatter(results.CoastFraction,results.EndSpeed_mps,50,results.GuidanceTime_s,'filled');
xlabel('Coast fraction');ylabel('Final speed (m/s)');title('Terminal speed');colorbar;
styleDark(gcf,bg,fg);
figure('Color',bg,'Name','V11 Lunar Trajectory');
[xs,ys,zs]=sphere(60);
surf(R*xs/1000,R*ys/1000,R*zs/1000,'FaceColor',[.45 .45 .45],...
    'EdgeColor','none','FaceAlpha',.85); hold on;
plot3(yCoast(:,1)/1000,yCoast(:,2)/1000,yCoast(:,3)/1000,'Color',[.15 .5 1],'LineWidth',1.7);
plot3(y(:,1)/1000,y(:,2)/1000,y(:,3)/1000,'r','LineWidth',2);
plot3(R*site(1)/1000,R*site(2)/1000,R*site(3)/1000,'gp','MarkerFaceColor','g','MarkerSize',13);
plot3(y(end,1)/1000,y(end,2)/1000,y(end,3)/1000,'mo','MarkerFaceColor','m','MarkerSize',7);
axis equal;grid on;view(35,25);xlabel('X (km)');ylabel('Y (km)');zlabel('Z (km)');
title('V11: Orbital transfer and guided powered descent');
legend('Moon','Transfer coast','Powered descent','Site07','Final position','Location','bestoutside');
styleDark(gcf,bg,fg);

%% Local functions
function [dydt,T]=guidedDynamics(t,y,mu,R,site,Tmax,Isp,g0,mDry,tgo)
    r=y(1:3); v=y(4:6); m=y(7);
    grav=-mu*r/norm(r)^3;
    target=R*site;
    % Finite-horizon position/velocity feedback; tgo decreases but is
    % bounded to avoid singular terminal commands.
    tau=max(tgo-t,35);
    aReq=-6*(r-target)/tau^2-4*v/tau-grav;
    % Do not command thrust inward through the surface near touchdown.
    u=r/norm(r);
    if norm(r)-R<2000
        aReq=aReq+max(0,-dot(aReq,u))*u;
    end
    aMag=norm(aReq);
    if m<=mDry || aMag<1e-10
        T=0; aThrust=zeros(3,1);
    else
        T=min(Tmax,m*aMag);
        aThrust=(T/m)*aReq/aMag;
    end
    dydt=[v;grav+aThrust;-T/(Isp*g0);norm(aThrust)];
end
function [value,isterminal,direction]=descentEvents(y,R,mDry)
    value=[norm(y(1:3))-R;y(7)-mDry];
    isterminal=[1;1];direction=[-1;-1];
end
function d=siteRange(r,site,R)
    u=r/norm(r);d=R*acos(max(-1,min(1,dot(u,site))));
end
function styleDark(fig,bg,fg)
    set(fig,'Color',bg);
    axesList=findall(fig,'Type','axes');
    for j=1:numel(axesList)
        ax=axesList(j);
        set(ax,'Color',[.085 .085 .085],'XColor',fg,'YColor',fg,...
            'ZColor',fg,'GridColor',[.48 .48 .48],'GridAlpha',.35,'FontSize',10);
        grid(ax,'on');
        set(get(ax,'Title'),'Color',fg);
        set(get(ax,'XLabel'),'Color',fg);
        set(get(ax,'YLabel'),'Color',fg);
        set(get(ax,'ZLabel'),'Color',fg);
    end
    legends=findall(fig,'Type','legend');
    for j=1:numel(legends)
        set(legends(j),'Color',[.10 .10 .10],'TextColor',fg,'EdgeColor',[.55 .55 .55]);
    end
end
