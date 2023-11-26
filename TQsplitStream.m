%% T-Q graph for the fluid-split system 
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
%This script generates Figure 21 in Section 4.3 of the Engineering Manual, 
%demonstrating the advantages of a sandTES system with a split fluid 
%stream.
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

Qdot=1e6;   %Transferred heat, total

QdotSplitCharge=0.3*Qdot;   %Transferred heat after split, charge
mDotSplitCharge=0.3;        %Water mass flow after split, charge

QdotSplitDis=0.7*Qdot;  %Transferred heat after split, discharge
mDotSplitDis=0.2;       %Water mass flow after split, discharge

mDotH2Ocharge=1;    %Water mass flow, charge
mDotH2Odis=0.6;     %Water mass flow, charge

Tsand=[480,100]+273.15;     %Sand temperature at outlet / inlet
TH2O=[500,80]+273.15;       %Water temperature at inlet / outlet
pH2O=250e5;                 %Water pressure


%% Calculations
QdotVec=linspace(0,Qdot,n);     %Heat flux as vector


mDotSand=Qdot/diff(fliplr(SiO2.h(Tsand)));  %Sand mass flow
hSand=SiO2.h(Tsand(1))-QdotVec./mDotSand;   %Sand specific enthalpy
TsandVec=SiO2.T_h(hSand);                   %Sand temperature as vector


hH2Ocharge=IF97.h(pH2O,TH2O(1))-QdotVec./mDotH2Ocharge;     %Water specific enthalpy, charge
TH2Ocharge=IF97.T_ph(pH2O,hH2Ocharge);                      %Water temperature, charge

%Water specific enthalpy, after split, charge
QdotSplitChargeVec=linspace(0,Qdot-QdotSplitCharge,n);
hH2OsplitCharge=IF97.h(pH2O,TH2O(1))-QdotSplitCharge./mDotH2Ocharge-QdotSplitChargeVec./mDotSplitCharge;

TH2OsplitCharge=IF97.T_ph(pH2O,hH2OsplitCharge);    %Water temperature, after split, charge


hH2Odis=IF97.h(pH2O,TH2O(2))+fliplr(QdotVec)./mDotH2Odis;   %Water specific enthalpy, discharge
TH2Odis=IF97.T_ph(pH2O,hH2Odis);                            %Water temperature, discharge

%Water specific enthalpy, after split, discharge
QdotSplitDisVec=linspace(0,Qdot-QdotSplitDis,n);
hH2OsplitDis=IF97.h(pH2O,TH2O(2))+QdotSplitDis./mDotH2Odis+fliplr(QdotSplitDisVec)./mDotSplitDis;

TH2OsplitDis=IF97.T_ph(pH2O,hH2OsplitDis);  %Water temperature, after split, discharge


%% Create figure
fig=figure(21);
clf(fig);
ax=gca();
colors=ax.ColorOrder;
hold(ax,'on');


plot(ax,QdotVec./Qdot,TsandVec-273.15,'Color',colors(3,:));

plot(ax,QdotVec./Qdot,TH2Ocharge-273.15,'Color',colors(2,:));
plot(ax,(QdotSplitCharge+QdotSplitChargeVec)./Qdot,TH2OsplitCharge-273.15,'Color',colors(2,:),'LineStyle','--');

plot(ax,QdotVec./Qdot,TH2Odis-273.15,'Color',colors(1,:));
plot(ax,QdotSplitDisVec./Qdot,TH2OsplitDis-273.15,'Color',colors(1,:),'LineStyle','--');


legItems=repmat(line(),5,1);
legItems(1)=plot(NaN,NaN,'Color',colors(3,:));
legItems(2)=plot(NaN,NaN,'Color',colors(2,:));
legItems(3)=plot(NaN,NaN,'Color',colors(1,:));
legItems(4)=plot(NaN,NaN,'Color','k');
legItems(5)=plot(NaN,NaN,'Color','k','LineStyle','--');


hold(ax,'off');

legend(ax,legItems,{'Sand','H_2O charge','H_2O discharge','Regular','With split'},'Location','best');

xlabel(ax,'Q/Q_{total} (-)');
ylabel(ax,'Temperature (°C)');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'TQsplitstream.tiff']);




