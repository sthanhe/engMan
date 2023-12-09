%% Generic T-Q graph
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
%This script generates Figure 5 in Section 3.2 of the Engineering Manual, 
%demonstrating basic features of T-Q graphs.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%   - Curve Fitting Toolbox, version 3.9
%Necessary files, classes, functions, and scripts:
%   - @IF97
%   - @SiO2


%% Parameters
n=1000;     %Number of cells for discretization

Qdot=1e6;   %Transferred heat
pinch=10;   %Pinch temperature

TsandIn=50+273.15;      %Sand inlet temperature
TsandIn2=90+273.15;     %Sand inlet temperature for the pinch-point demo
TH2O=[500,100]+273.15;  %Water temperature at inlet / outlet
pH2O=150e5;             %Water pressure


%% Calculations
Tsat=IF97.T_sat(pH2O);          %Water saturation temperature
QdotVec=linspace(0,Qdot,n);     %Heat flux as vector

mDotH2O=Qdot/diff(fliplr(IF97.h(pH2O,TH2O)));           %Water mass flow
hH2O=IF97.h(pH2O,TH2O(2))+fliplr(QdotVec)./mDotH2O;     %Specific enthalpy water
TH2O=IF97.T_ph(pH2O,hH2O);                              %Water temperature

Tsand=[IF97.T_sat(pH2O)-pinch,TsandIn];     %Sand temperature
Tsand2=[IF97.T_sat(pH2O)-pinch,TsandIn2];   %Sand temperature for pinch-point demo

ind=find([false,TH2O>Tsat] & [TH2O==Tsat,false]);   %Index to find pinch-point (approximation)
QdotSand=Qdot-QdotVec(ind);                         %Sand heat flux after pinch-point

mDotSand=QdotSand/diff(fliplr(SiO2.h(Tsand)));      %Sand mass flow
mDotSand2=QdotSand/diff(fliplr(SiO2.h(Tsand2)));    %Sand mass flow for pinch-point demo

hSand=SiO2.h(Tsand(2))+fliplr(QdotVec)./mDotSand;       %Sand specific enthalpy
hSand2=SiO2.h(Tsand2(2))+fliplr(QdotVec)./mDotSand2;    %Sand specific enthalpy for pinch-point demo

TsandVec=SiO2.T_h(hSand);       %Sand temperature
TsandVec2=SiO2.T_h(hSand2);     %Sand temperature for pinch-point demo


%% Create figure
fig=figure(5);
clf(fig);
ax=gca();
colors=ax.ColorOrder;
hold(ax,'on');


plot(ax,QdotVec.*10^-6,TH2O-273.15,'Color',colors(1,:));
plot(ax,QdotVec.*10^-6,TsandVec-273.15,'Color',colors(3,:));
plot(ax,QdotVec.*10^-6,TsandVec2-273.15,'Color',colors(3,:),'LineStyle','--');


hold(ax,'off');

legend(ax,{'H_2O','Sand','Pinch-point demo'},'Location','best');

xlabel(ax,'Q (MW)');
ylabel(ax,'Temperature (°C)');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'TQgeneric.tiff']);




