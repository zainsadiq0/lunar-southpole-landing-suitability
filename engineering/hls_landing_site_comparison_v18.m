% NOT a powered-descent simulation for each site, mission delta-V, or
% validated Starship HLS performance prediction.
clear; clc; close all;
mu=4.902800066e12; R=1737.4e3; hOrbit=100e3;
m0=250000; Isp=350; g0=9.80665;
% Reference site used in V17, fixed inertial coordinates (Moon nonrotating).
refLat=-88.811008; refLon=123.690068;
% EXAMPLE sites only. Replace these names/coordinates with project data.
siteNames={'V17 reference';'Illustrative A';'Illustrative B';'Illustrative C';'South pole'};
lat=[refLat;-89.0;-88.5;-87.5;-90.0];
lon=[refLon;0;90;180;0];
% Construct ONE fixed orbital plane, aligned with V17 reference geometry.
ref=unit18([cosd(refLat)*cosd(refLon);cosd(refLat)*sind(refLon);sind(refLat)]);
east=unit18([-sind(refLon);cosd(refLon);0]);
tangent=unit18(cross(east,ref));
planeNormal=unit18(cross(ref,tangent));
vCirc=sqrt(mu/(R+hOrbit));
N=numel(lat);crossTrack=zeros(N,1);inclination=zeros(N,1);
planeChangeDV=zeros(N,1);idealProp=zeros(N,1);
for k=1:N
    s=unit18([cosd(lat(k))*cosd(lon(k));cosd(lat(k))*sind(lon(k));sind(lat(k))]);
    % Minimum angular distance between this point and the fixed orbital great circle.
    inclination(k)=asind(min(1,abs(dot(planeNormal,s))));
    crossTrack(k)=R*deg2rad(inclination(k))/1000;
    % Proxy: single instantaneous rotation of a circular-orbit velocity
    % vector by this angle, assuming favorable node geometry.
    % NOT an optimized plane change, transfer delta-V, or mission lower bound.
    planeChangeDV(k)=2*vCirc*sind(inclination(k)/2);
    % Idealized mass-equivalent if this impulse occurred at mass m0.
    % Do NOT add this directly to V17 descent propellant: different mass epoch.
    idealProp(k)=m0*(1-exp(-planeChangeDV(k)/(Isp*g0)));
end
T=table(string(siteNames),lat,lon,inclination,crossTrack,planeChangeDV,idealProp,...
    'VariableNames',{'Site','Latitude_deg','Longitude_deg','PlaneOffset_deg',...
    'CrossTrack_km','CircularPlaneChangeProxy_mps','IdealizedMassEquivalent_kg'});
T=sortrows(T,'PlaneOffset_deg');
fprintf('\nSTARSHIP HLS V18 - COMMON ORBIT SITE ACCESSIBILITY\n');
fprintf('Shared reference orbit: 100 km circular, fixed plane through V17 site\n');
fprintf('Reference site: %.6f deg lat, %.6f deg lon\n',refLat,refLon);
fprintf('Circular speed: %.3f m/s\n',vCirc);
fprintf('V17 reference powered descent: 110850 kg propellant, 1.8028 m/s touchdown\n');
fprintf('Those descent metrics are NOT calculated for other sites.\n\n');
disp(T);
fprintf('\nWARNING: all nonreference site coordinates are illustrative.\n');
fprintf('Plane-change quantities are screening proxies, not transfer solutions.\n');
fprintf('No phasing, orbital timing, terrain, illumination, comms, lunar rotation,\n');
fprintf('full orbital maneuvers, or site-specific descent optimization modeled.\n');
filename='hls_v18_site_screening.csv';writetable(T,filename);
fprintf('Saved table: %s\n',filename);
%% Plots
bg=[.06 .06 .06];fg=[.9 .9 .9];
figure('Color',bg,'Name','V18 Shared Orbit Accessibility');
subplot(1,2,1);bar(categorical(T.Site),T.PlaneOffset_deg);
ylabel('Plane offset (deg)');title('Site offset from common orbit plane');
subplot(1,2,2);bar(categorical(T.Site),T.CircularPlaneChangeProxy_mps);
ylabel('Illustrative plane-change proxy (m/s)');title('Not full mission delta-V');
ax=findall(gcf,'Type','axes');
for j=1:numel(ax)
 set(ax(j),'Color',[.09 .09 .09],'XColor',fg,'YColor',fg,'GridColor',[.5 .5 .5]);
 set(get(ax(j),'Title'),'Color',fg);set(get(ax(j),'YLabel'),'Color',fg);
 grid(ax(j),'on');xtickangle(ax(j),25);
end
figure('Color',bg,'Name','V18 Site Coordinates');
scatter(lon,lat,90,inclination,'filled');hold on;
for k=1:N,text(lon(k)+3,lat(k),siteNames{k},'Color',fg,'FontSize',9);end
xlabel('Longitude (deg)');ylabel('Latitude (deg)');
title('Illustrative south-polar sites (longitude schematic)');colorbar;grid on;
set(gca,'Color',[.09 .09 .09],'XColor',fg,'YColor',fg,'GridColor',[.5 .5 .5]);
set(get(gca,'Title'),'Color',fg);
%% Helpers
function x=unit18(x)
 x=x/norm(x);
end
