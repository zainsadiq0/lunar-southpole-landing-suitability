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

%% V17: margin-aware trajectory selection and robustness
% Standalone educational model. This is not a flight-certified HLS simulation.
% All guidance, direction slew, and engine assumptions match V16.
fractions=0.923:0.001:0.930;
guidanceTimes=660:10:750;
slewDeg=2; minThrottle=.10; thrustRate=2e4;
maxBurn=2400; speedLimit=2; errorLimit=100;
reserveKg=5000; minFinalMass=mDry+reserveKg;
% Require a nominal touchdown speed <= 1.85 m/s for margin.
nominalSpeedCap=1.85;
% Monte Carlo tests are exploratory; distributions are illustrative.
nTrials=60; rng(17);
rows=[];
for i=1:numel(fractions)
    state=interp1(tCoast,yCoast,fractions(i)*Tcoast).';
    for j=1:numel(guidanceTimes)
        tgo=guidanceTimes(j);
        [tt,yy,ie]=simulate17(state,mStart,tgo,mu,R,site,Tmax,Isp,g0,...
            mDry,minThrottle,thrustRate,slewDeg,maxBurn,1,2e-8,1e-6);
        z=yy(end,:).'; speed=norm(z(4:6));
        err=siteRange16(z(1:3),site,R);
        success=any(ie==1)&&~any(ie==2)&&speed<=speedLimit&&...
            err<=errorLimit&&z(7)>=minFinalMass;
        rows(end+1,:)=[fractions(i),tgo,tt(end),err,speed,...
            mStart-z(7),z(7)-mDry,double(success),double(success&&speed<=nominalSpeedCap)];
    end
end
T=array2table(rows,'VariableNames',{'CoastFraction','GuidanceTime_s','Duration_s',...
    'SiteError_m','EndSpeed_mps','Propellant_kg','DryMassMargin_kg','Feasible','MarginQualified'});
qualified=sortrows(T(T.MarginQualified==1,:),'Propellant_kg');
if ~isempty(qualified)
    choice=qualified(1,:);selection='Margin-qualified minimum propellant';
else
    feasible=sortrows(T(T.Feasible==1,:),'EndSpeed_mps');
    if isempty(feasible)
        [~,k]=min(T.EndSpeed_mps);choice=T(k,:);selection='DIAGNOSTIC ONLY: no feasible solution';
    else
        choice=feasible(1,:);selection='Fallback: minimum touchdown speed among feasible';
    end
end
fraction=choice.CoastFraction;tgo=choice.GuidanceTime_s;
state=interp1(tCoast,yCoast,fraction*Tcoast).';
[tNom,yNom,ieNom]=simulate17(state,mStart,tgo,mu,R,site,Tmax,Isp,g0,...
    mDry,minThrottle,thrustRate,slewDeg,maxBurn,1,2e-8,1e-6);
nom=metrics17(yNom(end,:).',ieNom,site,R,mStart,mDry,speedLimit,errorLimit,minFinalMass);
fprintf('\nSTARSHIP HLS V17 - MARGIN-AWARE OPTIMIZATION\n');
fprintf('Candidates %d | feasible %d | margin-qualified %d\n',height(T),nnz(T.Feasible),nnz(T.MarginQualified));
fprintf('Selection: %s\n',selection);
fprintf('Fraction %.4f | guidance %.0f s | touchdown speed %.4f m/s | error %.3f m\n',fraction,tgo,nom(2),nom(3));
fprintf('Propellant %.1f kg | dry margin %.1f kg | success %d\n',nom(4),nom(5),nom(1));
fprintf('V16 minimum-propellant reference: 110144.6 kg, 1.972 m/s\n');
fprintf('Top margin-qualified candidates:\n');
disp(qualified(1:min(10,height(qualified)),:));

%% Solver convergence, same selected candidate
steps=[2 1 .5];rel=[1e-8 2e-9 1e-10];
conv=zeros(3,5);
for k=1:3
    [tc,yc,iec]=simulate17(state,mStart,tgo,mu,R,site,Tmax,Isp,g0,...
        mDry,minThrottle,thrustRate,slewDeg,maxBurn,steps(k),rel(k),rel(k)*100);
    m=metrics17(yc(end,:).',iec,site,R,mStart,mDry,speedLimit,errorLimit,minFinalMass);
    conv(k,:)=[steps(k),tc(end),m(2),m(3),m(4)];
end
convergence=array2table(conv,'VariableNames',{'MaxStep_s','Duration_s','EndSpeed_mps','SiteError_m','Propellant_kg'});
fprintf('\nSolver convergence:\n');disp(convergence);

%% Perturbation study: ignition and engine dispersions
% Gaussian 1-sigma variations: 10m position, 0.1m/s velocity,
% 100kg initial mass; 1%% thrust scale; 0.5%% specific impulse.
pert=zeros(nTrials,8);
for k=1:nTrials
    initialState=state;
    initialState(1:3)=initialState(1:3)+10*randn(3,1);
    initialState(4:6)=initialState(4:6)+.1*randn(3,1);
    massPert=mStart+100*randn;
    thrustScale=max(.90,1+.01*randn);
    ispScale=max(.90,1+.005*randn);
    [tp,yp,iep]=simulate17(initialState,massPert,tgo,mu,R,site,...
        Tmax*thrustScale,Isp*ispScale,g0,mDry,minThrottle,...
        thrustRate,slewDeg,maxBurn,1,2e-8,1e-6);
    met=metrics17(yp(end,:).',iep,site,R,massPert,mDry,speedLimit,errorLimit,minFinalMass);
    pert(k,:)=[k,met(1:5),thrustScale,ispScale];
end
perturbations=array2table(pert,'VariableNames',{'Trial','Success','EndSpeed_mps',...
    'SiteError_m','Propellant_kg','DryMassMargin_kg','ThrustScale','IspScale'});
fprintf('\nPerturbation study: %d / %d successful\n',nnz(perturbations.Success),nTrials);
fprintf('Worst speed %.4f m/s | worst error %.3f m | minimum dry margin %.1f kg\n',...
    max(perturbations.EndSpeed_mps),max(perturbations.SiteError_m),min(perturbations.DryMassMargin_kg));
if nnz(perturbations.Success)<nTrials
    fprintf('Unsuccessful perturbations (up to first 10):\n');
    failed=perturbations(perturbations.Success==0,:);
    disp(failed(1:min(10,height(failed)),:));
end

%% Figures
bg=[.06 .06 .06];fg=[.88 .88 .88];
figure('Color',bg,'Name','V17 Propellant versus touchdown speed');
scatter(T.EndSpeed_mps,T.Propellant_kg/1000,40,T.MarginQualified,'filled');hold on;
xline(speedLimit,'--','Speed limit');xline(nominalSpeedCap,':','Nominal target');
plot(nom(2),nom(4)/1000,'mp','MarkerSize',14,'MarkerFaceColor','m');
xlabel('Touchdown speed (m/s)');ylabel('Propellant (1000 kg)');
title('Fuel and touchdown-speed tradeoff');colorbar;styleDark16(gcf,bg,fg);
figure('Color',bg,'Name','V17 Robustness');
subplot(1,2,1);scatter(perturbations.EndSpeed_mps,perturbations.SiteError_m,35,perturbations.Success,'filled');hold on;
xline(speedLimit,'--');yline(errorLimit,'--');
xlabel('Touchdown speed (m/s)');ylabel('Landing error (m)');title('Perturbed landings');
subplot(1,2,2);histogram(perturbations.DryMassMargin_kg/1000,15);
xlabel('Dry mass margin (1000 kg)');ylabel('Cases');title('Remaining mass margin');
styleDark16(gcf,bg,fg);
figure('Color',bg,'Name','V17 Nominal descent');
alt=(sqrt(sum(yNom(:,1:3).^2,2))-R)/1000;
spd=sqrt(sum(yNom(:,4:6).^2,2));
subplot(2,2,1);plot(tNom,alt,'LineWidth',1.6);xlabel('Time (s)');ylabel('Altitude (km)');title('Altitude');
subplot(2,2,2);plot(tNom,spd,'LineWidth',1.6);xlabel('Time (s)');ylabel('Speed (m/s)');title('Speed');
subplot(2,2,3);plot(tNom,yNom(:,9)/1000,'LineWidth',1.6);xlabel('Time (s)');ylabel('Thrust (kN)');title('Thrust');
subplot(2,2,4);plot(tNom,yNom(:,7)/1000,'LineWidth',1.6);xlabel('Time (s)');ylabel('Mass (1000 kg)');title('Mass');
styleDark16(gcf,bg,fg);
fprintf('\nCAVEATS: spherical, nonrotating Moon; no terrain, navigation, true attitude dynamics, or gimbal.\n');
fprintf('Perturbation results apply only to the stated distributions; engine parameters are illustrative.\n');

%% Local functions
function [tt,yy,ie]=simulate17(state,mStart,tgo,mu,R,site,Tmax,Isp,g0,mDry,minThrottle,thrustRate,slewDeg,maxBurn,maxStep,relTol,absTol)
    q0=desiredDirection16(0,[state;mStart;0;0],mu,R,site,tgo);
    y0=[state;mStart;0;0;q0];
    opts=odeset('RelTol',relTol,'AbsTol',absTol,'MaxStep',maxStep,...
        'Events',@(t,y) events16(y,R,mDry));
    [tt,yy,~,~,ie]=ode45(@(t,y) dynamics16(t,y,mu,R,site,Tmax,Isp,g0,...
        mDry,tgo,minThrottle,thrustRate,slewDeg),[0 maxBurn],y0,opts);
end
function m=metrics17(z,ie,site,R,mStart,mDry,speedLimit,errorLimit,minFinalMass)
    speed=norm(z(4:6));err=siteRange16(z(1:3),site,R);
    success=any(ie==1)&&~any(ie==2)&&speed<=speedLimit&&err<=errorLimit&&z(7)>=minFinalMass;
    m=[double(success),speed,err,mStart-z(7),z(7)-mDry];
end

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
