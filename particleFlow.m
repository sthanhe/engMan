%% Bed level gradients and particle flow
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
%This script demonstrates the impacts of different factors on bed level 
%gradients and creates the Figures 62 and 63 in Section 10.1 of the
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
d_p=linspace(100,500,1000).*10^-6;  %Particle diameter
rho_p=2650;                         %Particle density

eps=0.5;    %Porosity at operating conditions
p0=1e5;     %Ambient pressure
FG=4;       %Degree of fluidization

l=1;        %Fluidized bed length
width=1;    %Fluidized bed widht
h0=2;       %Minimal bed level
href=1;     %Reference height

g=9.81;     %Gravitational acceleration


%% Bed level gradients, depending on bed temperature
%Specific parameters
T=linspace(100,850,4)'+273.15;  %Bed temperature
mDot=10;                        %Sand mass flow


%Calculations
mDotS=mDot./(href.*width);  %Specific sand mass flow

p=p0+FluBed.deltaP(h0./2,eps,rho_p);        %Reference bed pressure
wmf=FluBed.wmf(d_p,rho_p,p,T);              %Minimum fluidization velocity
w=FG.*wmf;                                  %Fluidization velocity
Gamma=FluBed.Gamma(w,href,d_p,rho_p,p,T);   %Mass diffusivity

pressGrad=mDotS./Gamma.*g.*href;            %Pressure gradient
deltaH=FluBed.h(pressGrad.*l,eps,rho_p);    %Bed level differences


%Create figure
fig=figure(62);
clf(fig);
ax=gca();

semilogy(ax,d_p.*10^6,flipud(deltaH)./l.*10^3);

grid(ax,'on');

legend(ax,compose('T=%.0f°C',flipud(T)-273.15),'Location','best');

xlabel(ax,'Particle diameter d_p (µm)');
ylabel(ax,'Bed level gradient \Deltah/\Deltax (mm/m)');

title(ax,sprintf('FG=%.0f, h$_0$=%.0f m, $\\dot m_S$=%.0f kg/m$^2$s',...
                    FG,h0,mDotS),'Interpreter','latex');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'levelGradTemp.tiff']);


%% Bed level gradients, depending on specific mass flow
%Specific parameters
T=400+273.15;               %Bed temperature
mDot=linspace(5,20,4)';     %Sand mass flow


%Calculations
mDotS=mDot./(href.*width);  %Specific sand mass flow

p=p0+FluBed.deltaP(h0./2,eps,rho_p);        %Reference bed pressure
wmf=FluBed.wmf(d_p,rho_p,p,T);              %Minimum fluidization velocity
w=FG.*wmf;                                  %Fluidization velocity
Gamma=FluBed.Gamma(w,href,d_p,rho_p,p,T);   %Mass diffusivity

pressGrad=mDotS./Gamma.*g.*href;            %Pressure gradient
deltaH=FluBed.h(pressGrad.*l,eps,rho_p);    %Bed level differences


%Create figure
fig=figure(63);
clf(fig);
ax=gca();

semilogy(ax,d_p.*10^6,flipud(deltaH)./l.*10^3);
% plot(ax,d_p.*10^6,flipud(deltaH).*10^3);

grid(ax,'on');

legend(ax,compose('$\\dot m_S$=%.0f kg/m$^2$s',flipud(mDotS)),...
    'Location','best','Interpreter','latex');

xlabel(ax,'Particle diameter d_p (µm)');
ylabel(ax,'Bed level gradient \Deltah/\Deltax (mm/m)');

title(ax,sprintf('FG=%.0f, h_0=%.0f m, T=%.0f°C',...
                    FG,h0,T-273.15));

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'levelGradMdotS.tiff']);




