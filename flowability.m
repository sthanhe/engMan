%% Demonstration of consolidation stress and unconfined yield strength
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
%This script creates the graph that demonstrates the relation between the 
%consolidation stress and the unconfined yield strength of a bulk of
%particles, Figure 15 of the Engineering Manual.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%Necessary files, classes, functions, and scripts:
%   - sigma_1.tiff, in the FigureInserts folder
%   - sigma_c.tiff, in the FigureInserts folder


%% Particle flowability demonstration
%Dummy stresses. x=consolidation stress, y=unconfined yields strength
x=linspace(0,0.8,1000);     %No time consolidation
DeltaX=0.1;                 %x-axis shift
DeltaY=0.2;                 %y-axis shift

x2=linspace(x(1)-DeltaX,x(end)+DeltaX,1000);    %After time consolidation
sigma1=sqrt(abs(x2))/3;                         %Unconfined yield strength


%Create figure
fig=figure(15);
clf(fig);
ax=gca();
colors=ax.ColorOrder;
hold(ax,'on');


plot(x2+x2(1),sigma1);
plot(x2+x2(1),sigma1+DeltaY);

plot(x,x,'k--');
plot(x,x./2,'k:','LineWidth',1);

hold(ax,'off');


ax.XLim=[x(1),x(end)];
ax.YLim=ax.XLim;


ax.XTick=[];
ax.XTickLabel=[];

ax.YTick=[];
ax.YTickLabel=[];


legend(ax,{'\sigma_c @ t_0','\sigma_c @ t_1','ff_c=1','ff_c=2'},'Location','north');


xlabel(ax,'Consolidation stress \sigma_1');
ylabel(ax,'Unconfined yield strength \sigma_c');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];


warning('off','imageio:tifftagsread:badTagValueDivisionByZero');
ax2=axes(fig,'Units','centimeters','Position',[14.5,1,0.9,1.6]);
imshow(['FigureInserts',filesep,'sigma_1.tiff'],'Parent',ax2,'InitialMagnification','fit','Interpolation','bilinear','Reduce',false);

ax3=axes(fig,'Units','centimeters','Position',[2.3,6.2,0.9,1.6]);
imshow(['FigureInserts',filesep,'sigma_c.tiff'],'Parent',ax3,'InitialMagnification','fit','Interpolation','bilinear','Reduce',false);
warning('on','imageio:tifftagsread:badTagValueDivisionByZero');


exportgraphics(fig,['Figures',filesep,'particleConsolidation.tiff']);




