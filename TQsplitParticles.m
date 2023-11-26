%% T-Q graph for the particle-split system 
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
%This script generates Figure 19 in Section 4.2.2 of the Engineering
%Manual, demonstrating the advantages of a sandTES system with a split
%particle stream.
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

mDotH2Ocharge=1.15;     %Water mass flow during charge
mDotH2Odis=1.22;        %Water mass flow during discharge

QdotSum=3e6;    %Transferred heat
splitDis=0.65;  %Fraction of the transferred heat at which the particles should split off

Tsand=[450,350,100]+273.15;     %Sand temperature at outlet / split / inlet
TH2O=[500,80]+273.15;           %Water temperature at inlet / outlet
pH2O=250e5;                     %Water pressure


%% Calculations
QdotDis=[splitDis,(1-splitDis)]*QdotSum;    %Transferred heat before / after split

split=round(n/2);                           %Index at which the split occurs
QdotVecSandCharge=linspace(0,QdotSum,n);    %Vector of heat fluxes for the sand, charge

%Vector of heat fluxes for the sand, discharge
QdotVecSandDis=linspace(0,QdotDis(1),split+1);
QdotVecSandDis=[QdotVecSandDis(1:end-1),linspace(QdotVecSandDis(end),QdotSum,split)];


%Sand side
mDotSandCharge=QdotSum/diff(SiO2.h(Tsand([3,1])));          %Sand mass flow, charge
hSand=SiO2.h(Tsand(1))-QdotVecSandCharge./mDotSandCharge;   %Specific sand enthalpies, charge
TsandVecCharge=SiO2.T_h(hSand);                             %Vector of sand temperatures, charge

mDotSandDis=QdotDis./-diff(SiO2.h(Tsand));  %Sand mass flow, discharge

%Specific sand enthalpies, discharge
hSand=SiO2.h(Tsand(1))-linspace(0,QdotDis(1),split)./mDotSandDis(1);
hSand=[hSand(1:end-1),hSand(end)-linspace(0,QdotDis(2),split+1)./mDotSandDis(2)];

TsandVecDis=SiO2.T_h(hSand);    %Vector of sand temperatures, discharge


%Water side
QdotVecH2O=linspace(0,QdotSum,n);                           %Heat flux water
hH2Ocharge=IF97.h(pH2O,TH2O(1))-QdotVecH2O./mDotH2Ocharge;  %Specific enthalpy water, charge
TH2Ocharge=IF97.T_ph(pH2O,hH2Ocharge);                      %Water temperature, charge


hH2Odis=IF97.h(pH2O,TH2O(2))+fliplr(QdotVecH2O)./mDotH2Odis;    %Specific enthalpy water, discharge
TH2Odis=IF97.T_ph(pH2O,hH2Odis);                                %Water temperature, discharge


%% Create figure
fig=figure(19);
clf(fig);
ax=gca();
colors=ax.ColorOrder;
hold(ax,'on');


plot(ax,QdotVecSandCharge./QdotSum,TsandVecCharge-273.15,'Color',colors(3,:));
plot(ax,QdotVecSandDis./QdotSum,TsandVecDis-273.15,'Color',colors(3,:),'LineStyle','--');
plot(ax,QdotVecH2O./QdotSum,TH2Ocharge-273.15,'Color',colors(2,:));
plot(ax,QdotVecH2O./QdotSum,TH2Odis-273.15,'Color',colors(1,:),'LineStyle','--');


hold(ax,'off');

legend(ax,{'Sand charge','Sand discharge','H_2O charge','H_2O discharge'},'Location','best');

xlabel(ax,'Q/Q_{total} (-)');
ylabel(ax,'Temperature (°C)');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'TQsplitParticles.tiff']);




