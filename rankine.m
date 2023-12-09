%% Rankine cycle T-s graph
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
%This script creates the T-s graph of a Rankine cycle included in Figure 2
%in Section 3.1 of the Engineering Manual.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%Necessary files, classes, functions, and scripts:
%   - @IF97


%% Parameters
n=100;  %Number of cells for discretization

Tcond=130+273.15;   %Condensation temperature
Tlive=550+273.15;   %Live steam temperature
plive=150e5;        %Live steam pressure
etaPump=0.9;        %Pump isentropic efficiency
etaTurb=0.7;        %Turbine isentropic efficiency


%% Points in the cycle
p1=IF97.p_sat(Tcond);
T1=Tcond;
s1=IF97.s(NaN,T1,0);
h1=IF97.h(NaN,T1,0);

p2=plive;
T2s=IF97.T_ps(p2,s1);
h2s=IF97.h(p2,T2s);
h2=(h2s-h1)./etaPump+h1;
T2=IF97.T_ph(p2,h2);
s2=IF97.s(p2,T2);

p3=plive;
T3=Tlive;
s3=IF97.s(p3,T3);
h3=IF97.h(p3,T3);


p4=p1;
if s3<=IF97.s(p4,NaN,1)
    x4s=1-(IF97.s(p4,NaN,1)-s3)./(IF97.s(p4,NaN,1)-IF97.s(p4,NaN,0));
    T4s=IF97.T_sat(p4);
    h4s=IF97.h(p4,NaN,x4s);
else
    x4s=NaN;
    T4s=IF97.T_ps(p4,s3);
    h4s=IF97.h(p4,T4s);
end

h4=h3-(h3-h4s).*etaTurb;
T4=IF97.T_ph(p4,h4);

if h4<=IF97.h(p4,NaN,1)
    x4=1-(IF97.h(p4,NaN,1)-h4)./(IF97.h(p4,NaN,1)-IF97.h(p4,NaN,0));
    s4=IF97.s(p4,NaN,x4);
else
    x4=NaN;
    s4=IF97.s(p4,T4);
end


%% Lines between points
sPump=linspace(s1,s2,n);
hPump=linspace(h1,h2,n);
Tpump=IF97.T_hs(sPump,hPump);


sSG=linspace(s2,s3,n);
TSG=IF97.T_ps(plive,sSG);


sTurb=linspace(s3,s4,n);
hTurb=linspace(h3,h4,n);
Tturb=IF97.T_hs(hTurb,sTurb);


sCond=linspace(s4,s1,n);
Tcond=IF97.T_ps(p4,sCond);


%% Boiling and evaporation lines
sSatL=linspace(IF97.s(NaN,IF97.T_t,0),IF97.s_c,n);
sSatV=linspace(IF97.s_c,IF97.s(NaN,IF97.T_t,1),n);

TsatL=arrayfun(@(i) ...
            fzero(@(T) ...
                IF97.s(NaN,T,0)-sSatL(i),...
            [IF97.T_t,IF97.T_c]),...
        1:length(sSatL));

TsatV=arrayfun(@(i) ...
            fzero(@(T) ...
                IF97.s(NaN,T,1)-sSatV(i),...
            [IF97.T_t,IF97.T_c]),...
        2:length(sSatV));
TsatV=[IF97.T_c,TsatV];


%% Create figure
fig=figure(2);
clf(fig);
ax=gca();
colors=ax.ColorOrder;
hold(ax,'on');

legItems=repmat(line(ax,'Visible','off'),4,1);

legItems(1)=plot(ax,sPump,Tpump-273.15,'Color',colors(4,:));
legItems(2)=plot(ax,sSG,TSG-273.15,'Color',colors(1,:));
legItems(3)=plot(ax,sTurb,Tturb-273.15,'Color',colors(2,:));
legItems(4)=plot(ax,sCond,Tcond-273.15,'Color',colors(3,:));

plot(ax,[sSatL,sSatV],[TsatL,TsatV]-273.15,'k');


hold(ax,'off');

legend(ax,legItems,{'1\rightarrow2','2\rightarrow3','3\rightarrow4','4\rightarrow1'},'Location','best');


annotation(fig,'arrow','HeadLength',5,'HeadWidth',5,...
    'Units','centimeters','Position',[6,5,1,0],'Color','k');

annotation(fig,'arrow','HeadLength',5,'HeadWidth',5,...
    'Units','centimeters','Position',[6,5,0,1],'Color','k');

text(ax,6.1,4.18,'s','Units','centimeters');
text(ax,4.86,5.4,'T','Units','centimeters');


ax.Visible='off';

fig.Units='centimeters';
fig.Position=[10,5,8,8];

exportgraphics(fig,['Figures',filesep,'rankine.tiff']);




