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

%% V14: constrained engine, refined search and robustness
% Educational lunar descent model; fixed inertial landing site and spherical Moon.
% Engine throttle floor is enforced when ON. Engine OFF is allowed.
% Thrust rate is limited using a continuous actuator state.
% Engine direction changes remain instantaneous (not attitude-constrained).
fractions=0.920:0.001:0.930;
tgoCandidates=680:5:780;
maxBurn=2400; speedLimit=2; targetSpeed=1.5; errorLimit=100;
reserveKg=5000; minFinalMass=mDry+reserveKg;
minThrottle=0.10; thrustRate=2e4; % 10% min ON, 20 kN/s (illustrative)
% Engine shutdown is commanded only when the required thrust is near zero.
% Thrust is continuous through the actuator, including shutdown.
rows=[]; trials=struct('score',{},'t',{},'y',{},'fraction',{},'tgo',{},...
    'success',{},'termination',{});
for i=1:numel(fractions)
    state=interp1(tCoast,yCoast,fractions(i)*Tcoast).';
    for j=1:numel(tgoCandidates)
        tgo=tgoCandidates(j);
        initial=[state;mStart;0;0]; % position, velocity, mass, deltaV, actual thrust
        op=odeset('RelTol',2e-7,'AbsTol',1e-5,'MaxStep',2,...
            'Events',@(t,y) descentEvents14(y,R,mDry));
        [tb,yb,~,~,ie]=ode45(@(t,y) guidedDynamics14(t,y,mu,R,site,Tmax,Isp,g0,mDry,...
            tgo,minThrottle,thrustRate),[0 maxBurn],initial,op);
        last=yb(end,:).';
        alt=norm(last(1:3))-R;
        err=siteRange14(last(1:3),site,R);
        speed=norm(last(4:6));
        contact=any(ie==1); fuel=any(ie==2);
        prop=mStart-last(7);
        success=contact && ~fuel && speed<=speedLimit && ...
            err<=errorLimit && last(7)>=minFinalMass;
        if contact, reason="Surface contact";
        elseif fuel, reason="Dry mass";
        else, reason="Time limit"; end
        radial=dot(last(4:6),last(1:3)/norm(last(1:3)));
        horiz=sqrt(max(0,speed^2-radial^2));
        rows(end+1,:)=[fractions(i),tgo,tb(end),alt,err,speed,prop,...
            last(8),last(7)-mDry,radial,horiz,double(contact),double(success)];
        % Prefer comfortable terminal speed among feasible candidates.
        score=err/1000+speed/10+max(alt,0)/1000+...
            1e4*double(fuel)+1e4*double(last(7)<minFinalMass);
        if success
            score=prop/1000+20*max(0,speed-targetSpeed)^2;
        else
            score=1e5+score;
        end
        trials(end+1)=struct('score',score,'t',tb,'y',yb,...
            'fraction',fractions(i),'tgo',tgo,'success',success,...
            'termination',reason);
    end
end
results=array2table(rows,'VariableNames',{'CoastFraction','GuidanceTime_s',...
    'BurnTime_s','EndAltitude_m','SiteError_m','EndSpeed_mps',...
    'Propellant_kg','PoweredDeltaV_mps','DryMassMargin_kg',...
    'RadialSpeed_mps','HorizontalSpeed_mps','Contact','Success'});
[~,idx]=min([trials.score]); best=trials(idx);
fprintf('\nSTARSHIP HLS V14 - ENGINE CONSTRAINTS AND ROBUSTNESS\n');
fprintf('Engine ON throttle floor %.0f%% | thrust rate %.0f kN/s\n',...
    100*minThrottle,thrustRate/1000);
fprintf('Candidates %d | feasible %d\n',height(results),nnz(results.Success));
fprintf('Best fraction %.4f | guidance time %.0f s | success %d\n',...
    best.fraction,best.tgo,best.success);
fprintf('Termination: %s\n',best.termination);
fprintf('Speed %.3f m/s | error %.3f m | propellant %.1f kg | margin %.1f kg\n',...
    norm(best.y(end,4:6)),siteRange14(best.y(end,1:3).',site,R),...
    mStart-best.y(end,7),best.y(end,7)-mDry);
if ~best.success
    warning('No feasible nominal candidate found with these illustrative engine constraints.');
end
feasible=sortrows(results(results.Success==1,:), 'Propellant_kg');
disp('Lowest propellant feasible candidates:');
disp(feasible(1:min(10,height(feasible)),:));

%% Robustness: deterministic perturbations of ignition state and mass
% This tests local sensitivity, NOT a probabilistic reliability claim.
Ntest=25; rng(14);
robust=zeros(Ntest,7);
baseState=interp1(tCoast,yCoast,best.fraction*Tcoast).';
for k=1:Ntest
    state=baseState;
    state(1:3)=state(1:3)+10*randn(3,1); % 10 m / axis (1 sigma)
    state(4:6)=state(4:6)+0.1*randn(3,1); % 0.1 m/s / axis
    mass0=mStart+100*randn; % 100 kg standard deviation
    initial=[state;mass0;0;0];
    op=odeset('RelTol',2e-7,'AbsTol',1e-5,'MaxStep',2,...
        'Events',@(t,y) descentEvents14(y,R,mDry));
    [~,yb,~,~,ie]=ode45(@(t,y) guidedDynamics14(t,y,mu,R,site,Tmax,Isp,g0,mDry,...
        best.tgo,minThrottle,thrustRate),[0 maxBurn],initial,op);
    last=yb(end,:).';
    err=siteRange14(last(1:3),site,R);
    speed=norm(last(4:6));
    ok=any(ie==1) && ~any(ie==2) && speed<=speedLimit && ...
        err<=errorLimit && last(7)>=minFinalMass;
    robust(k,:)=[k,err,speed,mass0-last(7),last(7)-mDry,...
        norm(last(1:3))-R,double(ok)];
end
robustTable=array2table(robust,'VariableNames',{'Trial','SiteError_m',...
    'EndSpeed_mps','Propellant_kg','DryMassMargin_kg','EndAltitude_m','Success'});
fprintf('Perturbation success: %d / %d\n',nnz(robustTable.Success),Ntest);
fprintf('Worst error %.2f m | worst speed %.3f m/s | minimum margin %.1f kg\n',...
    max(robustTable.SiteError_m),max(robustTable.EndSpeed_mps),...
    min(robustTable.DryMassMargin_kg));

%% Dark plotting style, consistent with V13
bg=[.06 .06 .06]; fg=[.88 .88 .88];
t=best.t; y=best.y;
alt=(sqrt(sum(y(:,1:3).^2,2))-R)/1000;
speed=sqrt(sum(y(:,4:6).^2,2));
range=zeros(numel(t),1);
radial=zeros(numel(t),1); horizontal=zeros(numel(t),1);
for k=1:numel(t)
    range(k)=siteRange14(y(k,1:3).',site,R)/1000;
    u=y(k,1:3).'/norm(y(k,1:3));
    radial(k)=dot(y(k,4:6),u);
    horizontal(k)=norm(y(k,4:6).'-radial(k)*u);
end
figure('Color',bg,'Name','V14 Descent Diagnostics');
subplot(3,2,1);plot(t,alt,'LineWidth',1.7);xlabel('Powered time (s)');ylabel('Altitude (km)');title('Altitude');
subplot(3,2,2);plot(t,speed,'LineWidth',1.7);xlabel('Powered time (s)');ylabel('Speed (m/s)');title('Inertial speed');
subplot(3,2,3);plot(t,range,'LineWidth',1.7);xlabel('Powered time (s)');ylabel('Range (km)');title('Surface distance to Site07');
subplot(3,2,4);plot(t,y(:,7)/1000,'LineWidth',1.7);xlabel('Powered time (s)');ylabel('Mass (1000 kg)');title('Vehicle mass');
subplot(3,2,5);plot(t,y(:,9)/1000,'LineWidth',1.7);xlabel('Powered time (s)');ylabel('Thrust (kN)');title('Actual thrust');
subplot(3,2,6);plot(t,y(:,8),'LineWidth',1.7);xlabel('Powered time (s)');ylabel('Delta-V (m/s)');title('Powered delta-V');
styleDark14(gcf,bg,fg);

figure('Color',bg,'Name','V14 Guidance Sweep');
subplot(1,2,1);scatter(results.CoastFraction,results.SiteError_m,32,results.GuidanceTime_s,'filled');
xlabel('Coast fraction');ylabel('Site error (m)');title('Landing accuracy');colorbar;
subplot(1,2,2);scatter(results.CoastFraction,results.EndSpeed_mps,32,results.GuidanceTime_s,'filled');
hold on;yline(speedLimit,'--');xlabel('Coast fraction');ylabel('Speed (m/s)');title('Terminal speed');colorbar;
styleDark14(gcf,bg,fg);

figure('Color',bg,'Name','V14 Terminal and Robustness');
subplot(2,2,1);plot(t,alt*1000,'LineWidth',1.7);xlim([max(0,t(end)-150),t(end)]);
xlabel('Powered time (s)');ylabel('Altitude (m)');title('Final 150 s');
subplot(2,2,2);plot(t,radial,'LineWidth',1.7);hold on;plot(t,horizontal,'LineWidth',1.7);
xlim([max(0,t(end)-150),t(end)]);xlabel('Powered time (s)');ylabel('Velocity (m/s)');
title('Terminal velocity components');legend('Radial','Horizontal');
subplot(2,2,3);plot(t,(y(:,7)-mDry)/1000,'LineWidth',1.7);hold on;
yline(reserveKg/1000,'--');xlabel('Powered time (s)');ylabel('Mass (1000 kg)');
title('Propellant above dry mass');
subplot(2,2,4);scatter(robustTable.SiteError_m,robustTable.EndSpeed_mps,40,...
    robustTable.Success,'filled');hold on;xline(errorLimit,'--');yline(speedLimit,'--');
xlabel('Site error (m)');ylabel('Touchdown speed (m/s)');title('Perturbed ignition trials');
styleDark14(gcf,bg,fg);

%% Local landing-site visualization: tangent plane, kilometers
figure('Color',bg,'Name','V14 Local Approach');
east=[-sind(lon);cosd(lon);0];east=east/norm(east);
north=cross(site,east);north=north/norm(north);
dr=y(:,1:3)-R*site.';
eastKm=dr*east/1000; northKm=dr*north/1000;
plot(eastKm,northKm,'r','LineWidth',1.8);hold on;
plot(0,0,'gp','MarkerFaceColor','g','MarkerSize',13);
plot(eastKm(end),northKm(end),'mo','MarkerFaceColor','m','MarkerSize',7);
axis equal; xlabel('Local east (km)');ylabel('Local north (km)');
title('V14: Local approach (tangent-plane projection)');
legend('Powered trajectory','Site07','Final point','Location','best');
styleDark14(gcf,bg,fg);

%% Dynamics and utilities
function [dydt,Tcmd]=guidedDynamics14(t,y,mu,R,site,Tmax,Isp,g0,mDry,tgo,minThrottle,thrustRate)
    r=y(1:3);v=y(4:6);m=y(7);Tact=max(0,min(Tmax,y(9)));
    grav=-mu*r/norm(r)^3;target=R*site;
    tau=max(tgo-t,35);
    aReq=-6*(r-target)/tau^2-4*v/tau-grav;
    u=r/norm(r);
    if norm(r)-R<2000
        aReq=aReq+max(0,-dot(aReq,u))*u;
    end
    aMag=norm(aReq);
    if m<=mDry || aMag<1e-10
        Tcmd=0;
    else
        Tcmd=min(Tmax,m*aMag);
        if Tcmd>0, Tcmd=max(minThrottle*Tmax,Tcmd);end
    end
    % Rate-limited continuous actuator; no discontinuous thrust jumps.
    dT=max(-thrustRate,min(thrustRate,20*(Tcmd-Tact)));
    if (Tact<=0 && dT<0) || (Tact>=Tmax && dT>0),dT=0;end
    if aMag>1e-10 && m>mDry
        aThrust=(Tact/m)*aReq/aMag;
    else
        aThrust=zeros(3,1);
    end
    dydt=[v;grav+aThrust;-Tact/(Isp*g0);norm(aThrust);dT];
end
function [value,isterminal,direction]=descentEvents14(y,R,mDry)
    value=[norm(y(1:3))-R;y(7)-mDry];
    isterminal=[1;1];direction=[-1;-1];
end
function d=siteRange14(r,site,R)
    u=r/norm(r);d=R*acos(max(-1,min(1,dot(u,site))));
end
function styleDark14(fig,bg,fg)
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
