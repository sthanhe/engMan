%% Porous plate design
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
%This script shows the design process of the porous plate floor as well as
%the impacts of different factors on the optimal design. It creates Figures
%51-54 in Section 9.1 of the Engineering Manual.
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


%% General arameters
n=75;   %Number of cells for discretization (in sand-flow direction)

rho_p=2650;     %Particle density
eps=0.5;        %Porosity at operating conditions

p0=1e5;         %Ambient pressure

FGdiffRel=0.1;  %Allowable difference in the degree of fluidization within each chamber, relative to FG_ref
deltaT=100;     %Bed temperature difference between the two sides (in sand flow direction)


%% Bed height pressure loss
%Specific parameters
d_p=200e-6;                     %Particle diameter
FG_ref=4;                       %Degree of fluidization in the (controlled) reference point
grade=10;                       %Filter grade
T0=linspace(850,100,4)+273.15;  %Minimum bed temperature
h0=linspace(0,4,n);             %Minimum bed level)
deltaH=100e-3;                  %Bed level difference inside chamber


%Derived parameters
T=[-1/2;1/2].*deltaT+T0;    %Bed temperature left (row 1) and right (row 2)
T=permute(T,[1,3,2]);

h=[deltaH;0]+h0;    %Bed level left (row 1) and right (row 2)

FGdiff=FGdiffRel*FG_ref;        %Allowable difference in the degree of fluidization within each chamber, absolute
FG=[-1/2;1/2].*FGdiff+FG_ref;   %Degree of fluidization left (row 1) and right (row 2)


%Calculations
[~,deltaP]=Sinter.s(d_p,rho_p,eps,p0,grade,T,FG,h);
deltaP=squeeze(deltaP)';


%Create figure
fig=figure(51);
clf(fig);
ax=gca();

plot(ax,h0,flipud(deltaP)*10^-2);

grid(ax,'on');

legend(ax,compose('T_0=%.0f°C',fliplr(T0)-273.15),'Location','best');

xlabel(ax,'Bed height h_0 (m)');
ylabel(ax,'Plate pressure loss \Deltap_{plate0} (mbar)');

title(ax,sprintf(['d_p=%.0f µm, FG_0=%.0f, \\DeltaFG=%.0f%%, ' ...
                    '\\Deltah=%.0f mm, Filter=%.0f'],...
                    d_p*10^6,FG_ref,FGdiffRel*100,deltaH*10^3,grade));

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'porousPressLoss.tiff']);


%% Bed height plate thickness and relative pressure loss
%Specific parameters
d_p=200e-6;             %Particle diameter
FG_ref=4;               %Degree of fluidization in the (controlled) reference point
grade=10;               %Filter grade
T0=400+273.15;          %Minimum bed temperature
h0=linspace(0,4,n);     %Minimum bed level)
deltaH=100e-3;          %Bed level difference inside chamber


%Derived parameters
T=[-1/2;1/2].*deltaT+T0;    %Bed temperature left (row 1) and right (row 2)
T=permute(T,[1,3,2]);

h=[deltaH;0]+h0;    %Bed level left (row 1) and right (row 2)

FGdiff=FGdiffRel*FG_ref;        %Allowable difference in the degree of fluidization within each chamber, absolute
FG=[-1/2;1/2].*FGdiff+FG_ref;   %Degree of fluidization left (row 1) and right (row 2)


%Calculations
[s,deltaP]=Sinter.s(d_p,rho_p,eps,p0,grade,T,FG,h);
deltaPbed=FluBed.deltaP(h0,eps,rho_p);


%Create figure
fig=figure(52);
clf(fig);
ax=gca();

plot(ax,h0,s*10^3);
ylabel(ax,'Plate thickness s (mm)');

yyaxis(ax,'right');
plot(ax,h0,deltaP./deltaPbed);
ylabel(ax,'\Deltap_{plate0} / \Deltap_{bed0} (-)');

grid(ax,'on');

xlabel(ax,'Bed height h_0 (m)');

title(ax,sprintf(['d_p=%.0f µm, FG_0=%.0f, \\DeltaFG=%.0f%%, ' ...
                    'T_0=%.0f°C, \\Deltah=%.0f mm, Filter=%.0f'],...
                    d_p*10^6,FG_ref,FGdiffRel*100,T0-273.15,deltaH*10^3,grade));

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'porousPlateThick.tiff']);


%% Level difference
%Specific parameters
d_p=200e-6;                             %Particle diameter
FG_ref=4;                               %Degree of fluidization in the (controlled) reference point
grade=permute([15,10,7,5],[1,3,2]);     %Filter grade
T0=400+273.15;                          %Minimum bed temperature
h0=2;                                   %Minimum bed level)
deltaH=linspace(1e-3,200e-3,n);         %Bed level difference inside chamber


%Derived parameters
T=[-1/2;1/2].*deltaT+T0;    %Bed temperature left (row 1) and right (row 2)
T=permute(T,[1,3,2]);

h=[deltaH;zeros(1,n)]+h0;    %Bed level left (row 1) and right (row 2)

FGdiff=FGdiffRel*FG_ref;        %Allowable difference in the degree of fluidization within each chamber, absolute
FG=[-1/2;1/2].*FGdiff+FG_ref;   %Degree of fluidization left (row 1) and right (row 2)


%Calculations
s=Sinter.s(d_p,rho_p,eps,p0,grade,T,FG,h);
s=squeeze(s);


%Create figure
fig=figure(53);
clf(fig);
ax=gca();

plot(ax,deltaH*10^3,s*10^3);

grid(ax,'on');

legend(ax,compose('Filter=%.0f',grade),'Location','best');

xlabel(ax,'Level difference \Deltah (mm)');
ylabel(ax,'Plate thickness s (mm)');

title(ax,sprintf(['d_p=%.0f µm, FG_0=%.0f, \\DeltaFG=%.0f%%, ' ...
                    'T_0=%.0f°C, h_0=%.0f m'],...
                    d_p*10^6,FG_ref,FGdiffRel*100,T0-273.15,h0));

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'porousLevelDiff.tiff']);


%% Particle diameter
%Specific parameters
d_p=permute(linspace(100e-6,500e-6,n),[1,3,2]);                             %Particle diameter
FG_ref=2:6;                               %Degree of fluidization in the (controlled) reference point
grade=5;     %Filter grade
T0=400+273.15;                          %Minimum bed temperature
h0=2;                                   %Minimum bed level)
deltaH=100e-3;         %Bed level difference inside chamber


%Derived parameters
T=[-1/2;1/2].*deltaT+T0;    %Bed temperature left (row 1) and right (row 2)
T=permute(T,[1,3,2]);

h=[deltaH;0]+h0;    %Bed level left (row 1) and right (row 2)

FGdiff=FGdiffRel*FG_ref;        %Allowable difference in the degree of fluidization within each chamber, absolute
FG=[-1/2;1/2].*FGdiff+FG_ref;   %Degree of fluidization left (row 1) and right (row 2)


%Calculations
s=Sinter.s(d_p,rho_p,eps,p0,grade,T,FG,h);
s=squeeze(s);


%Create figure
fig=figure(54);
clf(fig);
ax=gca();

plot(ax,squeeze(d_p)*10^6,s*10^3);

grid(ax,'on');

legend(ax,compose('FG_0=%.0f',FG_ref),'Location','best');

xlabel(ax,'Particle diameter d_p (µm)');
ylabel(ax,'Plate thickness s (mm)');

title(ax,sprintf(['\\DeltaFG=%.0f%%, T_0=%.0f°C, h_0=%.0f m, ' ...
                    '\\Deltah=%.0f mm, Filter=%.0f'],...
                    FGdiffRel*100,T0-273.15,h0,deltaH*10^3,grade));

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'porousParticleDiameter.tiff']);




