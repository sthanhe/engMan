%% Aspects of minimum fluidization velocity
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the sandTES Engineering Manual
%
%All required files for this script can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.xxxxxxx
% 
%All parameters and results are in SI base units.
%
%
%
%This script demonstrates the impacts of different factors on the minimum
%fluidization velocity and creates the Figures 8-10 in Section 3.4.1 of the
%Engineering Manual.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%   - Curve Fitting Toolbox, version 3.9
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - @FluBed
%   - @SiO2
%   - @Sinter
%   - @implExp


%% Parameters
d_p=linspace(100,700,4)'.*10^-6;    %Particle diameter
rho_p=2650;                         %Particle density
phi_s=0.8;                          %Particle sphericity

T=linspace(0,900,1000)+273.15;  %Bed temperature
p=1e5;                          %Bed pressure


%% Calculations
wmf=FluBed.wmf(d_p,rho_p,p,T);          %Minimum fluidization velocity
wt=FluBed.w_t(d_p,rho_p,phi_s,p,T);     %Terminal velocity


%% Figure 8: Minimum fluidization over temperature
fig=figure(8);
clf(fig);
ax=gca();

plot(ax,T-273.15,flipud(wmf));

ax.XLim=[min(T),max(T)]-273.15;

legend(ax,compose('d_p=%.0f µm',flipud(d_p).*10^6),'Location','best');

xlabel(ax,'Temperature T (°C)');
ylabel(ax,'Minimum fluidization velocity w_{mf} (m/s)');

title(ax,sprintf('p=%.0f bar, \\rho_p=%.0f kg/m³',p.*10^-5,rho_p));

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'minFluid.tiff']);


%% Figure 9: Terminal velocity over temperature
fig=figure(9);
clf(fig);
ax=gca();

plot(ax,T-273.15,flipud(wt));

ax.XLim=[min(T),max(T)]-273.15;

legend(ax,compose('d_p=%.0f µm',flipud(d_p).*10^6),'Location','best');

xlabel(ax,'Temperature T (°C)');
ylabel(ax,'Terminal velocity w_{t} (m/s)');

title(ax,sprintf('p=%.0f bar, \\rho_p=%.0f kg/m³, \\phi_s=%.2f',...
                p.*10^-5,rho_p,phi_s));

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'terminalVelocity.tiff']);


%% Figure 10: Useful velocity range over temperature
fig=figure(10);
clf(fig);
ax=gca();

plot(ax,T-273.15,wt./wmf);

ax.XLim=[min(T),max(T)]-273.15;

legend(ax,compose('d_p=%.0f µm',d_p.*10^6),'Location','best');

xlabel(ax,'Temperature T (°C)');
ylabel(ax,'Useful velocity range w_{t}/w_{mf} (-)');

title(ax,sprintf('p=%.0f bar, \\rho_p=%.0f kg/m³, \\phi_s=%.2f',...
                p.*10^-5,rho_p,phi_s));

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'usefulVelocityRange.tiff']);




