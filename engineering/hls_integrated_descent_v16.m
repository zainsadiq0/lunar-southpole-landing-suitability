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

%% V16: thrust-direction slew feasibility
% Standalone; uses V14 lunar transfer and spacecraft constants.
% Direction slew is a proxy for attitude dynamics, NOT a gimbal model.
% The 10% throttle floor is a commanded ON-state approximation.
fractions=0.923:0.001:0.930;
guidanceTimes=660:10:750;
slewDeg=2;                  % illustrative maximum direction rate (deg/s)
minThrottle=0.10; thrustRate=2e4;
maxBurn=2400; speedLimit=2; errorLimit=100;
reserveKg=5000; minFinalMass=mDry+reserveKg;
rows=[]; trials=struct('score',{},'t',{},'y',{},'fraction',{},'tgo',{},'success',{});
for i=1:numel(fractions)
    state=interp1(tCoast,yCoast,fractions(i)*Tcoast).';
    for j=1:numel(guidanceTimes)
        tgo=guidanceTimes(j);
        q0=desiredDirection16(0,[state;mStart;0;0],mu,R,site,tgo);
        initial=[state;mStart;0;0;q0]; % 12 states
        opts16=odeset('RelTol',2e-8,'AbsTol',1e-6,'MaxStep',1,...
            'Events',@(t,y) events16(y,R,mDry));
        [tt,yy,~,~,ie]=ode45(@(t,y) dynamics16(t,y,mu,R,site,Tmax,Isp,g0,...
            mDry,tgo,minThrottle,thrustRate,slewDeg),[0 maxBurn],initial,opts16);
        z=yy(end,:).'; speed=norm(z(4:6));
        err=siteRange16(z(1:3),site,R);
        success=any(ie==1)&&~any(ie==2)&&speed<=speedLimit&&...
            err<=errorLimit&&z(7)>=minFinalMass;
        prop=mStart-z(7);
        score=1e6+speed+err/100;
        if success,score=prop;end
        rows(end+1,:)=[fractions(i),tgo,tt(end),err,speed,prop,...
            z(8),z(7)-mDry,double(any(ie==1)),double(success)];
        trials(end+1)=struct('score',score,'t',tt,'y',yy,...
            'fraction',fractions(i),'tgo',tgo,'success',success);
    end
end
results16=array2table(rows,'VariableNames',{'CoastFraction','GuidanceTime_s',...
    'BurnTime_s','SiteError_m','EndSpeed_mps','Propellant_kg',...
    'PoweredDeltaV_mps','DryMassMargin_kg','Contact','Success'});
[~,idx]=min([trials.score]);best16=trials(idx);
fprintf('\nSTARSHIP HLS V16 - THRUST DIRECTION SLEW STUDY\n');
fprintf('Direction rate limit %.2f deg/s (illustrative)\n',slewDeg);
fprintf('Trials %d | feasible %d\n',height(results16),nnz(results16.Success));
fprintf('Selected fraction %.4f | guidance %.0f s | success %d\n',...
    best16.fraction,best16.tgo,best16.success);
z=best16.y(end,:).';
fprintf('Touchdown/termination speed %.3f m/s | error %.3f m\n',...
    norm(z(4:6)),siteRange16(z(1:3),site,R));
fprintf('Powered propellant %.1f kg | delta-V %.1f m/s | dry margin %.1f kg\n',...
    mStart-z(7),z(8),z(7)-mDry);
fprintf('V14 selected propellant baseline: 111715.1 kg\n');
if ~best16.success
    warning('No feasible landing found with the imposed direction-rate constraint.');
end
disp('Ten lowest-fuel feasible candidates:');
feasible16=sortrows(results16(results16.Success==1,:),'Propellant_kg');
disp(feasible16(1:min(10,height(feasible16)),:));

%% Diagnostics
bg=[.06 .06 .06];fg=[.88 .88 .88];
tt=best16.t;yy=best16.y;n=numel(tt);
alt=(sqrt(sum(yy(:,1:3).^2,2))-R)/1000;
speed=sqrt(sum(yy(:,4:6).^2,2));
rangeKm=zeros(n,1);angleErr=zeros(n,1);dirRate=zeros(n,1);
for k=1:n
    rangeKm(k)=siteRange16(yy(k,1:3).',site,R)/1000;
    q=yy(k,10:12).';q=q/norm(q);
    d=desiredDirection16(tt(k),yy(k,:).',mu,R,site,best16.tgo);
    angleErr(k)=acosd(max(-1,min(1,dot(q,d))));
    dy=dynamics16(tt(k),yy(k,:).',mu,R,site,Tmax,Isp,g0,...
        mDry,best16.tgo,minThrottle,thrustRate,slewDeg);
    dirRate(k)=norm(dy(10:12))*180/pi;
end
figure('Color',bg,'Name','V16 Descent');
subplot(3,2,1);plot(tt,alt,'LineWidth',1.6);xlabel('Time (s)');ylabel('Altitude (km)');title('Altitude');
subplot(3,2,2);plot(tt,speed,'LineWidth',1.6);xlabel('Time (s)');ylabel('Speed (m/s)');title('Speed');
subplot(3,2,3);plot(tt,rangeKm,'LineWidth',1.6);xlabel('Time (s)');ylabel('Range (km)');title('Site range');
subplot(3,2,4);plot(tt,yy(:,7)/1000,'LineWidth',1.6);xlabel('Time (s)');ylabel('Mass (1000 kg)');title('Mass');
subplot(3,2,5);plot(tt,yy(:,9)/1000,'LineWidth',1.6);xlabel('Time (s)');ylabel('Thrust (kN)');title('Actual thrust');
subplot(3,2,6);plot(tt,yy(:,8),'LineWidth',1.6);xlabel('Time (s)');ylabel('Delta-V (m/s)');title('Powered delta-V');
styleDark16(gcf,bg,fg);
figure('Color',bg,'Name','V16 Direction Diagnostics');
subplot(1,2,1);plot(tt,angleErr,'LineWidth',1.6);xlabel('Time (s)');ylabel('Angle (deg)');title('Desired vs actual thrust direction');
subplot(1,2,2);plot(tt,dirRate,'LineWidth',1.6);hold on;yline(slewDeg,'--');
xlabel('Time (s)');ylabel('Direction rate (deg/s)');title('Actual thrust direction rate');
styleDark16(gcf,bg,fg);
figure('Color',bg,'Name','V16 Feasibility');
subplot(1,2,1);scatter(results16.CoastFraction,results16.EndSpeed_mps,30,results16.GuidanceTime_s,'filled');
hold on;yline(speedLimit,'--');xlabel('Coast fraction');ylabel('End speed (m/s)');title('Guidance sweep');colorbar;
subplot(1,2,2);scatter(results16.Propellant_kg/1000,results16.SiteError_m,30,results16.Success,'filled');
hold on;yline(errorLimit,'--');xlabel('Propellant (1000 kg)');ylabel('Site error (m)');title('Feasibility');colorbar;
styleDark16(gcf,bg,fg);
figure('Color',bg,'Name','V16 Terminal Approach');
plot(rangeKm,alt,'LineWidth',1.8);hold on;
plot(rangeKm(end),alt(end),'mo','MarkerFaceColor','m');
xlabel('Surface range to site (km)');ylabel('Altitude (km)');
title('Altitude versus remaining surface range');styleDark16(gcf,bg,fg);

%% Local functions
function d=desiredDirection16(t,y,mu,R,site,tgo)
    r=y(1:3);v=y(4:6);grav=-mu*r/norm(r)^3;
    tau=max(tgo-t,35);
    aReq=-6*(r-R*site)/tau^2-4*v/tau-grav;
    u=r/norm(r);
    if norm(r)-R<2000,aReq=aReq+max(0,-dot(aReq,u))*u;end
    if norm(aReq)<1e-12,d=u;else,d=aReq/norm(aReq);end
end
function dy=dynamics16(t,y,mu,R,site,Tmax,Isp,g0,mDry,tgo,minThrottle,thrustRate,slewDeg)
    r=y(1:3);v=y(4:6);m=y(7);
    Tact=max(0,min(Tmax,y(9)));
    grav=-mu*r/norm(r)^3;
    tau=max(tgo-t,35);
    aReq=-6*(r-R*site)/tau^2-4*v/tau-grav;
    u=r/norm(r);
    if norm(r)-R<2000,aReq=aReq+max(0,-dot(aReq,u))*u;end
    if m<=mDry || norm(aReq)<1e-10
        Tcmd=0;
    else
        Tcmd=min(Tmax,m*norm(aReq));
        Tcmd=max(minThrottle*Tmax,Tcmd);
    end
    dT=max(-thrustRate,min(thrustRate,20*(Tcmd-Tact)));
    if (Tact<=0&&dT<0)||(Tact>=Tmax&&dT>0),dT=0;end
    d=desiredDirection16(t,y,mu,R,site,tgo);
    q=y(10:12);q=q/max(norm(q),eps);
    tangent=d-dot(d,q)*q;
    nt=norm(tangent);
    if nt<1e-12
        qdot=zeros(3,1);
    else
        theta=atan2(nt,dot(d,q));
        omega=min(slewDeg*pi/180,2*theta);
        qdot=omega*tangent/nt;
    end
    aThrust=(Tact/max(m,eps))*q;
    dy=[v;grav+aThrust;-Tact/(Isp*g0);norm(aThrust);dT;qdot];
end
function [value,isterminal,direction]=events16(y,R,mDry)
    value=[norm(y(1:3))-R;y(7)-mDry];
    isterminal=[1;1];direction=[-1;-1];
end
function d=siteRange16(r,site,R)
    u=r/norm(r);d=R*acos(max(-1,min(1,dot(u,site))));
end
function styleDark16(fig,bg,fg)
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
