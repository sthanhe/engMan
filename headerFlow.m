%% Vertical header pressures
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
%This script calculates the distribution of pressure inside vertical
%headers depending on the flow direction and creates the graphs included in
%Figure 40 as well as Figure 41 of Section 7.4 of the Engineering Manual.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%Necessary files, classes, functions, and scripts:
%   - @IF97


%% Calculate pressure distributions
[pBotWater,pTopWater]=makeFig(200+273.15,1,'headersWater');
[pBotSteam,pTopSteam,z]=makeFig(500+273.15,2,'headersSteam');

%Setup lines, correct for necessary pressure losses
pBotSteam=pBotSteam-(pBotSteam(end)-pBotWater(end));
pTopSteam=pTopSteam-(pTopSteam(end)-pTopWater(end));
pTopSteam=pTopSteam-max(pTopSteam-pTopWater);


%Create figure
fig=figure(41);
clf(fig);
ax=gca();
colors=ax.ColorOrder;
hold(ax,'on');

plot(ax,pTopWater.*10^-5,z,'Color',colors(1,:));
plot(ax,pTopSteam.*10^-5,z,'Color',colors(1,:));
plot(ax,pBotWater.*10^-5,fliplr(z),'Color',colors(2,:));
plot(ax,pBotSteam.*10^-5,fliplr(z),'Color',colors(2,:));

xline(ax,pTopWater(1).*10^-5);


legItems=repmat(line(ax,'Visible','off'),2,1);
legItems(1)=plot(ax,NaN,NaN,'Color',colors(1,:));
legItems(2)=plot(ax,NaN,NaN,'Color',colors(2,:));

legend(ax,legItems,{'Top','Bottom'},'Location','best');

hold(ax,'off');


fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

annotation(fig,'arrow','Units','centimeters','Position',[12.9,8.5/2,1,-1.5],'Color',colors(1,:));
annotation(fig,'arrow','Units','centimeters','Position',[12.3,4.5,-1,1.5],'Color',colors(1,:));

annotation(fig,'arrow','Units','centimeters','Position',[8.1,4,-1,1.5],'Color',colors(2,:));
annotation(fig,'arrow','Units','centimeters','Position',[4.8,3.5,-1,-1.5],'Color',colors(2,:));


xlabel(ax,'Pressure (bar)');
ylabel(ax,'Position in header (m)');

exportgraphics(fig,['Figures',filesep,'headersComparison.tiff']);


%% Calculation function and addin figures for Figure 40
function [pBot,pTop,z]=makeFig(T,figidx,figname)
    %Parameters
    mDotTotal=25;   %Total mass flow inside header
    pIn=150e5;      %Inlet pressure
    
    A=200e-3.^2.*pi/4;  %Header cross section
    z=[1,0];            %Vertical coordinates of header in- / outlets
    g=9.81;             %Gravitational acceleration
    n=100;              %Number of cells for discretization


    %Derived parameters
    z=linspace(z(1),z(2),n);    %Vertical variable (z-axis)
    zRev=fliplr(z);             %Reversed z-axis

    rho=IF97.v(pIn,T).^-1;          %Fluid density
    v=mDotTotal./(rho.*A);          %Fluid velocity
    mDot=linspace(mDotTotal,0,n);   %Fluid mass flow
    
    
    %Top inlet / outlet
    bern=v.^2./2+g.*z(1)+pIn./rho;                  %Bernoulli equation constant
    pTop=(bern-(mDot./rho./A).^2./2-g.*z).*rho;     %Pressure distribution
    pTopStat=pIn+rho.*(g.*z(1)-g.*z);               %Hydrostatic part
    
    
    % Bottom inlet / outlet
    bern=v.^2./2+g.*zRev(1)+pIn./rho;               %Bernoulli equation constant
    pBot=(bern-(mDot./rho./A).^2./2-g.*zRev).*rho;  %Pressure distribution
    pBotStat=pIn-rho.*g.*z;                         %Hydrostatic part
    
    
    %Create addins for Figure 40
    fig=figure(figidx);
    clf(fig);
    ax=gca();
    colors=ax.ColorOrder;
    hold(ax,'on');
    
    plot(ax,pTop.*10^-5,z);
    plot(ax,pBot.*10^-5,zRev);
    plot(ax,pTopStat.*10^-5,z,'Color',colors(1,:),'LineStyle','--');
    plot(ax,pBotStat.*10^-5,z,'Color',colors(2,:),'LineStyle','--');


    pos=ax.Position;
    x0=pos(1)+pos(3)/2;
    y0=pos(2);

    annotation(fig,'arrow','HeadLength',5,'HeadWidth',5,...
        'Position',[x0,y0,0,pos(4)],'Color','k');

    annotation(fig,'arrow','HeadLength',5,'HeadWidth',5,...
        'Position',[x0,y0,pos(3)/2,0],'Color','k');

    
    hold(ax,'off');
    
    
    ax.XLim=[pIn-0.1e5,pIn+0.1e5].*10^-5;
    
    ax.Visible='off';
    
    fig.Units='centimeters';
    fig.Position=[10,5,6,6];
    
    exportgraphics(fig,['Figures',filesep,figname,'.tiff']);
end







