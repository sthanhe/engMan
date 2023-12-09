%% Fluidization gas distribution between gas boxes
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
%This script calculates the distribution of fluidization gas for the 
%example given in Section 9.2, Gas boxes and orifice plates, of the 
%Engineering Manual and creates Figures 57-60.
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
%   - @GasBox
%   - @Orifice
%   - @SiO2
%   - @Sinter
%   - @implExp


%% Parameters
nDisc=51;   %Number of cells for discretization (in sand-flow direction)

d_p=200e-6;     %Particle diameter
rho_p=2650;     %Particle density
eps=0.5;        %Porosity at operating conditions

p0=1e5;         %Ambient pressure

l=[0.25,1,1,0.25];  %Gas box lengths (in sand-flow direction), from left to right
width=1;            %Fluidized bed width
h0=2;               %Fluidized bed base height (minimum bed level)

AC=[false,true,true,false];     %Vector whether gas box has an air cushion valve, from left to right

nSupply=[1,2,2,1];      %Number of supply tubes per gas box, from left to right
D=88.9e-3-2*3.2e-3;     %Inside diameter of fluidization gas supply tubes

grade=10;   %Filter grade

deltaT=[0,200,200,0];   %Bed temperature difference inside each chamber, from left to right
T0=200+273.15;          %Minimum bed temperature

FG_ref=4;       %Degree of fluidization in the (controlled) reference point
FGdiffRel=0.1;  %Allowable difference in the degree of fluidization within each chamber, relative to FG_ref

deltaH=linspace(25,100,4)'.*10^-3;  %Bed level difference inside largest chamber (linear interpolation according to chamber length for the rest)

deltaPbase=linspace(10,300,100).*10^2;      %Basic orifice plate pressure loss
deltaPbase=sort([deltaPbase,100e2,200e2]);  %Include 100 mbar and 200 mbar for plots


%% Setup, design calculations
%Mean temperatures inside chambers
T0=T0+cumsum(deltaT);
T0(2:end)=mean([T0(1:end-1);T0(2:end)]);

%Start (row 1) and end (row 2) temperatures of each chamber, increasing temperature from left to right
T=[-1/2;1/2].*deltaT+T0;

%Setup cell of gas box configurations and calculate the design cases
gbCell=cell(length(deltaH),length(deltaPbase));
for i=1:length(deltaH)
    for j=1:length(deltaPbase)
        gb(length(l))=GasBox(); %#ok<SAGROW>
        for k=1:length(l)
            gb(k)=GasBox(d_p,rho_p,eps,p0,l(k),width,h0,AC(k),nSupply(k),D,grade);
            gb(k).nDisc=nDisc;
            gb(k).T=T(:,k);
        end
    
    
        gb.design(FG_ref,FGdiffRel,deltaH(i),deltaPbase(j));
    
        gbCell{i,j}=gb;
        clear('gb');
    end
end


%% Pressure distribution in the design case, Figure 57
idx1=find(deltaH==100e-3);
idx2=find(deltaPbase==100e2);
gb=gbCell{idx1,idx2};

[fig,ax]=gb.plotDeltaP(57);

title(ax,sprintf('h_0=%.0f m, \\Deltah=%.0f mm, \\Deltap_{base}=%.0f mbar',...
                    h0,deltaH(idx1)*10^3,deltaPbase(idx2).*10^-2));

exportgraphics(fig,['Figures',filesep,'gasBoxDesignDeltaP.tiff']);


%% FG distribution in the design case, Figure 58
[fig,ax]=gb.plotFG(58);

title(ax,sprintf('h_0=%.0f m, \\Deltah=%.0f mm',...
                    h0,deltaH(idx1)*10^3));

exportgraphics(fig,['Figures',filesep,'gasBoxDesignFG.tiff']);


%% FG distribution in the discharge case, Figure 59
%Switch direction, simulate gas box configurations with chamber 2 as the control chamber
for i=1:numel(gbCell)
    gbCell{i}.direction('discharge');
    gbCell{i}.sim(FG_ref,2);
end


%Setup x and y (=FG) variables for the plot
x=cumsum([0,l]);
x=cell2mat(...
        arrayfun(@(i) ...
            linspace(x(i),x(i+1),nDisc),...
    1:length(l),'UniformOutput',false));

idx=find(any(deltaPbase==[10;100;200].*10^2,1));
y=cell2mat(...
    cellfun(@(x) reshape([x.FG],1,numel([x.FG]))',...
    gbCell(end,idx),'UniformOutput',false));
y=fliplr(y);


%Create figure
fig=figure(59);
clf(fig);
ax=gca();

plot(ax,x,y);

legend(ax,compose('\\Deltap_{base}=%.0f mbar',fliplr(deltaPbase(idx)).*10^-2),'Location','best');

title(ax,sprintf('h_0=%.0f m, \\Deltah=%.0f mm',...
                    h0,deltaH(end)*10^3));

xlabel(ax,'Distance from the left x (m)');
ylabel(ax,'Degree of fluidization FG (-)');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'gasBoxDischargeFG.tiff']);


%% FG difference depending on basic pressure loss, Figure 60
%Find FG differences
FGcell=cellfun(@(x) [x.FG],gbCell,'UniformOutput',false);
deltaFG=flipud(cellfun(@(x) ...
            max(x(:,[2,3]),[],'all')-min(x(:,[2,3]),[],'all'),...
        FGcell));


%Create figure
fig=figure(60);
clf(fig);
ax=gca();

plot(ax,deltaPbase.*10^-2,deltaFG);

legend(ax,compose('\\Deltah=%.0f mm',flipud(deltaH).*10^3),'Location','best');

title(ax,sprintf('h_0=%.0f m',h0));

xlabel(ax,'Basic pressure loss \Deltap_{base} (mbar)');
ylabel(ax,'\DeltaFG (-)');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'gasBoxDischargeDeltaFG.tiff']);




