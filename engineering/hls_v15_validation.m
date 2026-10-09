% V15 intentionally does not claim attitude, terrain or flight validation.
% Uses V14 functions as local copies below for independent replay.
if ~exist('best','var') || ~exist('tCoast','var') || ~exist('yCoast','var')
    error('Run hls_integrated_descent_v14.m before this validation script.');
end
fprintf('\nSTARSHIP HLS V15 - VALIDATION STUDY\n');
baseState=interp1(tCoast,yCoast,best.fraction*Tcoast).';
maxSteps=[4 2 1 0.5]; relTols=[2e-7 1e-8 1e-9 1e-10];
conv=zeros(numel(maxSteps),8);
for k=1:numel(maxSteps)
    [tt,yy,ie]=replay15([baseState;mStart;0;0],best.tgo,mu,R,site,Tmax,Isp,g0,mDry,...
        minThrottle,thrustRate,maxBurn,maxSteps(k),relTols(k),max(1e-8,relTols(k)*50));
    z=yy(end,:).';
    conv(k,:)=[maxSteps(k),relTols(k),tt(end),siteRange15(z(1:3),site,R),...
        norm(z(4:6)),mStart-z(7),max(abs(diff(yy(:,9)))./max(diff(tt),eps)),double(any(ie==1))];
end
convergence=array2table(conv,'VariableNames',{'MaxStep_s','RelTol','Duration_s',...
    'SiteError_m','EndSpeed_mps','Propellant_kg','ObservedThrustRate_Nps','Contact'});
disp(convergence);
ref=conv(end,:);
fprintf('Finest-vs-baseline: site error %.4f m, speed %.5f m/s, propellant %.2f kg\n',...
    abs(conv(2,4)-ref(4)),abs(conv(2,5)-ref(5)),abs(conv(2,6)-ref(6)));

%% Deterministic sensitivity tests, with larger disturbances than V14
rng(15); N=100; pert=zeros(N,11);
for k=1:N
    s=baseState;
    s(1:3)=s(1:3)+50*randn(3,1);       % position: 50 m / axis, 1 sigma
    s(4:6)=s(4:6)+0.5*randn(3,1);      % velocity: 0.5 m/s / axis
    mInit=mStart+500*randn;             % initial mass: 500 kg sigma
    thrustScale=1+0.02*randn;          % 2% thrust uncertainty
    ispScale=1+0.01*randn;             % 1% Isp uncertainty
    actualTmax=Tmax*thrustScale;
    actualIsp=Isp*ispScale;
    [~,yy,ie]=replay15([s;mInit;0;0],best.tgo,mu,R,site,actualTmax,actualIsp,g0,mDry,...
        minThrottle,thrustRate,maxBurn,2,2e-7,1e-5);
    z=yy(end,:).'; err=siteRange15(z(1:3),site,R); sp=norm(z(4:6));
    ok=any(ie==1) && ~any(ie==2) && sp<=speedLimit && err<=errorLimit && z(7)>=minFinalMass;
    pert(k,:)=[k,err,sp,mInit-z(7),z(7)-mDry,thrustScale,ispScale,...
        norm(z(1:3))-R,double(any(ie==1)),double(any(ie==2)),double(ok)];
end
robust15=array2table(pert,'VariableNames',{'Trial','SiteError_m','EndSpeed_mps',...
    'Propellant_kg','DryMassMargin_kg','ThrustScale','IspScale','EndAltitude_m',...
    'Contact','DryMassReached','Success'});
fprintf('Expanded perturbations: %d / %d successful\n',nnz(robust15.Success),N);
fprintf('Worst site error %.2f m | worst speed %.3f m/s | minimum dry mass margin %.1f kg\n',...
    max(robust15.SiteError_m),max(robust15.EndSpeed_mps),min(robust15.DryMassMargin_kg));

%% Diagnostics: dark background, no white figure margins
bg=[.06 .06 .06]; fg=[.88 .88 .88];
figure('Color',bg,'Name','V15 Numerical Convergence');
subplot(1,3,1);plot(conv(:,1),conv(:,4),'o-','LineWidth',1.5);xlabel('Max step (s)');ylabel('Site error (m)');title('Position convergence');
subplot(1,3,2);plot(conv(:,1),conv(:,5),'o-','LineWidth',1.5);xlabel('Max step (s)');ylabel('Speed (m/s)');title('Speed convergence');
subplot(1,3,3);plot(conv(:,1),conv(:,6)/1000,'o-','LineWidth',1.5);xlabel('Max step (s)');ylabel('Propellant (1000 kg)');title('Fuel convergence');
styleDark15(gcf,bg,fg);
figure('Color',bg,'Name','V15 Expanded Robustness');
subplot(1,2,1);scatter(pert(:,2),pert(:,3),30,pert(:,11),'filled');hold on;
xline(errorLimit,'--');yline(speedLimit,'--');xlabel('Site error (m)');ylabel('Touchdown speed (m/s)');title('Robustness outcomes');
subplot(1,2,2);histogram(pert(:,5)/1000,20);hold on;xline(reserveKg/1000,'--');
xlabel('Mass above dry (1000 kg)');ylabel('Trials');title('Propellant reserve');
styleDark15(gcf,bg,fg);

% Terminal local trajectory: surface distance and altitude rather than a
% global tangent-plane projection, which is misleading over hundreds of km.
figure('Color',bg,'Name','V15 Terminal Approach');
y=best.y; rnorm=sqrt(sum(y(:,1:3).^2,2));
rangeKm=zeros(size(rnorm));
for k=1:numel(rnorm),rangeKm(k)=siteRange15(y(k,1:3).',site,R)/1000;end
plot(rangeKm,(rnorm-R)/1000,'LineWidth',1.8);hold on;
plot(rangeKm(end),(rnorm(end)-R)/1000,'mo','MarkerFaceColor','m');
xlabel('Surface distance to Site07 (km)');ylabel('Altitude (km)');
title('V15: Descent profile versus target range');styleDark15(gcf,bg,fg);

fprintf('NOTE: This is a simplified fixed-site, spherical-Moon model.\n');
fprintf('No lunar rotation, terrain, attitude dynamics, gimbal, or navigation model.\n');
fprintf('Rate limiting is checked numerically; thrust-direction slew is not constrained.\n');

function [t,y,ie]=replay15(initial,tgo,mu,R,site,Tmax,Isp,g0,mDry,minThrottle,thrustRate,maxBurn,maxStep,rtol,atol)
 op=odeset('RelTol',rtol,'AbsTol',atol,'MaxStep',maxStep,...
     'Events',@(t,y) events15(y,R,mDry));
 [t,y,~,~,ie]=ode45(@(t,y) dynamics15(t,y,mu,R,site,Tmax,Isp,g0,mDry,tgo,minThrottle,thrustRate),...
     [0 maxBurn],initial,op);
end
function dydt=dynamics15(t,y,mu,R,site,Tmax,Isp,g0,mDry,tgo,minThrottle,thrustRate)
 r=y(1:3);v=y(4:6);m=y(7);Tact=max(0,min(Tmax,y(9)));
 grav=-mu*r/norm(r)^3;target=R*site;tau=max(tgo-t,35);
 aReq=-6*(r-target)/tau^2-4*v/tau-grav;
 u=r/norm(r);
 if norm(r)-R<2000,aReq=aReq+max(0,-dot(aReq,u))*u;end
 aMag=norm(aReq);
 if m<=mDry || aMag<1e-10,Tcmd=0;
 else
     Tcmd=min(Tmax,m*aMag);
     if Tcmd>0,Tcmd=max(minThrottle*Tmax,Tcmd);end
 end
 dT=max(-thrustRate,min(thrustRate,20*(Tcmd-Tact)));
 if (Tact<=0 && dT<0)||(Tact>=Tmax && dT>0),dT=0;end
 if aMag>1e-10 && m>mDry,aThrust=(Tact/m)*aReq/aMag;
 else,aThrust=zeros(3,1);end
 dydt=[v;grav+aThrust;-Tact/(Isp*g0);norm(aThrust);dT];
end
function [value,isterminal,direction]=events15(y,R,mDry)
 value=[norm(y(1:3))-R;y(7)-mDry];isterminal=[1;1];direction=[-1;-1];
end
function d=siteRange15(r,site,R)
 u=r/norm(r);d=R*acos(max(-1,min(1,dot(u,site))));
end
function styleDark15(fig,bg,fg)
 set(fig,'Color',bg); ax=findall(fig,'Type','axes');
 for j=1:numel(ax)
  set(ax(j),'Color',[.085 .085 .085],'XColor',fg,'YColor',fg,'ZColor',fg,...
      'GridColor',[.48 .48 .48],'GridAlpha',.35,'FontSize',10);
  grid(ax(j),'on');set(get(ax(j),'Title'),'Color',fg);
  set(get(ax(j),'XLabel'),'Color',fg);set(get(ax(j),'YLabel'),'Color',fg);
 end
end
