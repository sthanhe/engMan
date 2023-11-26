%% Particle size distribution and entrainment
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
%This script investigates the impact of particle size distribution on
%entrainment, Figures 33-35 of the Engineering Manual.
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


%% Impact of temperature, Figure 33
%Parameters
d_p=linspace(1e-6,1000e-6,1000);    %Particle diameter
rho_p=2650;                         %Particle density
phi_s=0.75;                         %Particle sphericity
eps_mf=0.45;                        %Porosity at minimum fluidization
p=1e5;                              %Pressue
T=linspace(100,850,4)'+273.15;      %Temperature
FG=4;                               %Degree of fluidization


%Calculate particle size (and below) that gets elutriated
d_pMin=FluBed.elutriation(d_p,rho_p,phi_s,eps_mf,p,T,FG);


%Create figure
fig=figure(33);
clf(fig);
ax=gca();

plot(ax,d_p.*10^6,d_pMin./d_p);

grid(ax,'on');

legend(ax,compose('T=%.0f°C',T-273.15),'Location','best');

xlabel(ax,'Mean particle diameter (µm)');
ylabel(ax,'Fraction of elutriable fines (-)');

title(ax,sprintf(['\\rho_p=%.0f kg/m³, \\phi_s=%.2f, \\epsilon_{mf}=%.2f, ' ...
                    'p=%.0f bar, FG=%.0f'],rho_p,phi_s,eps_mf,p*10^-5,FG));

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'entrainmentTemp.tiff']);


%% Impact of the degree of fluidization, Figure 34
%Parameters
d_p=linspace(1e-6,1000e-6,1000);    %Particle diameter
rho_p=2650;                         %Particle density
phi_s=0.75;                         %Particle sphericity
eps_mf=0.45;                        %Porosity at minimum fluidization
p=1e5;                              %Pressue
T=400+273.15;                       %Temperature
FG=(2:6)';                          %Degree of fluidization


%Calculate particle size (and below) that gets elutriated
d_pMin=FluBed.elutriation(d_p,rho_p,phi_s,eps_mf,p,T,FG);


%Create figure
fig=figure(34);
clf(fig);
ax=gca();

plot(ax,d_p.*10^6,flipud(d_pMin./d_p));

grid(ax,'on');

legend(ax,compose('FG=%.0f',flipud(FG)),'Location','best');

xlabel(ax,'Mean particle diameter (µm)');
ylabel(ax,'Fraction of elutriable fines (-)');

title(ax,sprintf(['\\rho_p=%.0f kg/m³, \\phi_s=%.2f, \\epsilon_{mf}=%.2f, ' ...
                    'p=%.0f bar, T=%.0f°C'],rho_p,phi_s,eps_mf,p*10^-5,T-273.15));

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'entrainmentFG.tiff']);


%% Impact of bed height, Figure 35
%Parameters
d_p=linspace(100e-6,1000e-6,4)';    %Particle diameter
rho_p=2650;                         %Particle density
phi_s=0.75;                         %Particle sphericity
eps=0.5;                            %Porosity at operating conditions
p0=1e5;                             %Ambient pressue
T=350+273.15;                       %Temperature
FG=4;                               %Degree of fluidization
h=linspace(0,4,1000);               %Bed height


%Fluidization velocity at the bed's center
p=p0+FluBed.deltaP(h./2,eps,rho_p);     %Bed pressure
rhoA=DryAir.rho(p,T);                   %Air density
wmf=FluBed.wmf(d_p,rho_p,p,T);          %Minimum fluidization velocity
w=FG.*wmf;                              %Fluidization velocity


%Degree of fluidization at the bed's top (index 0)
rhoA0=DryAir.rho(p0,T);             %Air density
wmf0=FluBed.wmf(d_p,rho_p,p0,T);    %Minimum fluidization velocity
w0=w.*rhoA./rhoA0;                  %Fluidization velocity
FG0=w0./wmf0;                       %Degree of fluidization


%Create figure
fig=figure(35);
clf(fig);
ax=gca();

plot(ax,h,FG0./FG);

grid(ax,'on');

legend(ax,compose('d_p=%.0f µm',d_p.*10^6),'Location','best');

xlabel(ax,'Bed height (m)');
ylabel(ax,'FG_0 / FG_{center} (-)');

title(ax,sprintf(['\\rho_p=%.0f kg/m³, \\epsilon=%.1f, ' ...
                    'p_0=%.0f bar'],rho_p,eps,p0*10^-5));

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'entrainmentBedHeight.tiff']);




