%% Ideal particle diameter for plain tubes
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
%This script calculates the ideal particle diameter when using plain tubes
%for the  example given in Section 6.2 of the Engineering Manual and 
%creates Figures 31 and 32.
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
%   - @IF97
%   - @SiO2
%   - @Sinter
%   - @implExp
%   - IMfigure.m


%% Parameters
n=1000;         %Number of cells for discretization
DeltaX=0.01;    %Cell length in tube direction

mDotH2O=0.6;    %Water mass flow
Qdot=1e6;       %Transferred heat

Tsand=[480,100]+273.15;     %Sand temperatures at outlet / inlet
TH2O=[500,80]+273.15;       %Water temperatures at inlet / outlet
pH2O=250e5;                 %Water pressure

p=1e5;  %Bed pressure

d_p=linspace(120e-6,1000e-6,1000)';     %Particle diameters
rho_p=SiO2.rho(300);                    %Particle density

FG=4;           %Degree of fluidization
eps_mf=0.45;    %Porosity at minimum fluidization


%Derived parameters
l=n*DeltaX;     %Heat exchanger tube length


%% Temperature difference based on the T-Q graph
QdotVec=linspace(0,Qdot,n);

mDotSand=Qdot/diff(fliplr(SiO2.h(Tsand)));
hSand=SiO2.h(Tsand(1))-QdotVec./mDotSand;
Tsand=SiO2.T_h(hSand);

hH2O=IF97.h(pH2O,TH2O(2))+fliplr(QdotVec)./mDotH2O;
TH2O=IF97.T_ph(pH2O,hH2O);

DeltaT=Tsand-TH2O;


%% Calculate ideal particle diameter
w=FG.*FluBed.wmf(d_p,rho_p,p,Tsand);                        %Fluidization velocity
alpha=FluBed.molerus(w,Tsand,Tsand,p,d_p,rho_p,eps_mf);     %Outside heat transfer coefficient
qDot=alpha.total.*DeltaT;                                   %Specific heat flux

IM=DeltaX*trapz(qDot,2)./(l.*Qdot);     %Integral mean

%Find optimum
[~,idx]=max(IM);
d_pIdeal=d_p(idx);


%% Create figures
%T-Q graph
fig=figure(31);
clf(fig);
ax=gca();
colors=ax.ColorOrder;
hold(ax,'on');


plot(ax,QdotVec./Qdot,Tsand-273.15,'Color',colors(3,:));
plot(ax,QdotVec./Qdot,TH2O-273.15,'Color',colors(1,:));


hold(ax,'off');

legend(ax,{'Sand','H_2O'},'Location','best');

xlabel(ax,'Q/Q_{total} (-)');
ylabel(ax,'Temperature (°C)');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'idealDp_TQ.tiff']);


%Integral mean
fig=IMfigure(d_p,IM,d_pIdeal,32);

exportgraphics(fig,['Figures',filesep,'idealDp.tiff']);




