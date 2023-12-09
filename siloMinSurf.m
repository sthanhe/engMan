%% Silos with minimal surface area
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the sandTES Engineering Manual
%
%All required files for this script can be found in the software
%repository:
%https://doi.org/10.5281/ZENODO.10207330
% 
%All parameters and results are in SI base units.
%
%
%
%This script calculates the dimensions of silos to achieve minimum possible
%surface area for given loads and sizes. It creates Figures 68-70 and 
%calculates the dimensions shown in Figure 71 of Section 11.2 of the 
%Engineering Manual.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%   - Parallel Computing Toolbox, version 7.8
%Necessary files, classes, functions, and scripts:
%   - @Silo
%   - @implExp


%% Auxiliary structures and functions
beta=struct('name','beta','val',[]);    %Hopper half angle
phi_r=struct('name','phi_r','val',[]);  %Angle or repose
gamma=struct('name','gamma','val',[]);  %Roof angle


%Function handles to generate strings for graphs
degstr=@(ang) ['\',ang.name,'=',num2str(rad2deg(ang.val)),'°, '];   %Degrees
etastr=@(eta) ['\eta_{free}=',num2str(eta)];                        %Free surface effectivity


%% Cylinder diameter over volume for different hopper diameters
%Parameters
V=logspace(5,0,100);    %Silo volume
d_h=(2:-0.5:0.5)';      %Hopper outlet diameter
beta.val=deg2rad(15);   %Hopper half angle
phi_r.val=deg2rad(39);  %Angle of repose
gamma.val=deg2rad(15);  %Roof angle
eta_free=1;             %Free surface effectivity

idx=3;   %Control variable, which hopper diameter should go on the linear scale on the right


%Calculations
[~,d_cyl]=Silo.minSurf(V,d_h,beta.val,phi_r.val,gamma.val,eta_free);


%Create figure
fig=figure(68);
clf(fig);
ax=gca();

lines=loglog(ax,V,d_cyl);
ylabel(ax,'log_{10} (d_{cyl}) (m)');

hold(ax,'on');
yyaxis(ax,'right');
semilogx(ax,V,d_cyl(idx,:));
ylabel(ax,'Linear d_{cyl} (m)');
hold(ax,'off');

grid(ax,'on');

legend(ax,lines,compose('d_{hop}=%.1fm',d_h),'Location','northwest');

xlabel(ax,'Storage volume V (m³)');
title(ax,[degstr(beta),degstr(phi_r),degstr(gamma),etastr(eta_free)]);

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'silo_minSurfDcOverV.tiff']);


%% Slenderness over hopper half angle
%Parameters
V=1e3;                                  %Silo volume
d_h=1;                                  %Hopper outlet diameter
beta.val=deg2rad(linspace(45,5,100));   %Hopper half angle
phi_r.val=deg2rad([45,39,35,30])';      %Angle of repose
gamma.val=deg2rad(15);                  %Roof angle
eta_free=1;                             %Free surface effectivity


%Calculations
[~,d_cyl,h_cyl]=Silo.minSurf(V,d_h,beta.val,phi_r.val,gamma.val,eta_free);
slenderness=h_cyl./d_cyl;


%Create figure
fig=figure(69);
clf(fig);
ax=gca();

plot(ax,rad2deg(beta.val),slenderness);

grid(ax,'on');

legend(ax,compose('\\phi_r=%.0f°',rad2deg(phi_r.val)),'Location','best');

xlabel(ax,'Hopper half angle \beta (°)');
ylabel(ax,'Slenderness h_{cyl} / d_{cyl} (-)');

Vstr=['V=',num2str(V),'m³, '];
d_hStr=['d_{hop}=',num2str(d_h),'m, '];
title(ax,[Vstr,d_hStr,degstr(gamma),etastr(eta_free)]);

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'silo_minSurfSlenderOverBeta.tiff']);


%% Slenderness over roof angle
%Parameters
V=1e3;                                  %Silo volume
d_h=1;                                  %Hopper outlet diameter
beta.val=deg2rad(15);                   %Hopper half angle
phi_r.val=deg2rad(39);                  %Angle of repose
gamma.val=deg2rad(linspace(30,0,100));  %Roof angle
eta_free=(1.2:-0.1:0.8)';               %Free surface effectivity


%Calculations
[~,d_cyl,h_cyl]=Silo.minSurf(V,d_h,beta.val,phi_r.val,gamma.val,eta_free);
slenderness=h_cyl./d_cyl;


%Create figure
fig=figure(70);
clf(fig);
ax=gca();

plot(ax,rad2deg(gamma.val),slenderness);

grid(ax,'on');

legend(ax,compose('\\eta_{free}=%.1f',eta_free),'Location','best');

xlabel(ax,'Roof angle \gamma (°)');
ylabel(ax,'Slenderness h_{cyl} / d_{cyl} (-)');

title(ax,[Vstr,d_hStr,degstr(beta),degstr(phi_r)]);

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'silo_minSurfSlenderOverGamma.tiff']);


%% Silo dimensions for 1000 m³ of sand (Figure 71)
%Parameters
V=1e3;                  %Silo volume
d_h=1;                  %Hopper outlet diameter
beta.val=deg2rad(15);   %Hopper half angle
phi_r.val=deg2rad(39);  %Angle of repose
gamma.val=deg2rad(15);  %Roof angle
eta_free=1;             %Free surface effectivity


%Calculations
[~,d_cyl,h_cyl,h_hop,h_cylEff,h_heap,h_r]=Silo.minSurf(V,d_h,beta.val,phi_r.val,gamma.val,eta_free);


H=h_hop+h_cyl+h_r;                                  %Total height
Vhop=h_hop.*pi./12.*(d_cyl.^2+d_cyl.*d_h+d_h.^2);   %Hopper volume
VcylEff=d_cyl.^2.*pi./4.*h_cylEff;                  %Effective cylinder volume
Vheap=d_cyl.^2.*pi./12.*h_heap;                     %Heap volume




