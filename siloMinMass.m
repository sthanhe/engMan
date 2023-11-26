%% Silos with minimal steel mass
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
%This script calculates the dimensions of silos to achieve minimum possible
%steel mass for given loads and sizes. It creates Figures 78-81 of Section
%11.3.5 of the Engineering Manual.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%   - Parallel Computing Toolbox, version 7.8
%Necessary files, classes, functions, and scripts:
%   - @Silo
%   - @implExp


%% Auxiliary structures and functions
beta=struct('name','beta','val',[]);    %Hopper half angle
phi_r=struct('name','phi_r','val',[]);  %Angle or repose
gamma=struct('name','gamma','val',[]);  %Roof angle


%Function handles to generate strings for graphs
degstr=@(ang) ['\',ang.name,'=',num2str(rad2deg(ang.val)),'°, '];   %Degrees
Tstr=@(T) ['T=',num2str(T-273.15),'°C, '];                          %Temperature
steelstr=@(steel) ['Material=',steel,', '];                         %Wall material
wallstr=@(wall) ['Wall=D',num2str(wall)];                           %Wall type


%% Volume over cylindrical diameter
%Parameters
V=logspace(5,0,100);    %Silo volume
d_h=(2:-0.5:0.5)';      %Hopper outlet diameter
beta.val=deg2rad(15);   %Hopper half angle
phi_r.val=deg2rad(39);  %Angle of repose
gamma.val=deg2rad(15);  %Roof angle
T=400+273.15;           %Storage temperature
steel='S235';           %Wall material
wall=2;                 %Wall type
particles='Sand';       %Storage material

idx=3;  %Control variable, which angle of repose should go on the linear scale on the right


%Calculations
d_cyl=Silo.minMass(V,d_h,beta.val,phi_r.val,gamma.val,T,steel,wall,particles);


%Create figure
fig=figure(78);
clf(fig);
ax=gca();

lines=loglog(ax,V,d_cyl);
ylabel(ax,'log_{10} (d_{cyl}) (m)');

hold(ax,'on');
yyaxis(ax,'right');
semilogx(ax,V,d_cyl(idx,:));
ylabel(ax,'Linear d_{cyl} (m)');
hold(ax,'off');

grid(ax,'on');

legend(ax,lines,compose('d_{hop}=%.1fm',d_h),'Location','northwest');

xlabel(ax,'Storage volume V (m³)');
title(ax,[degstr(beta),degstr(phi_r),degstr(gamma),Tstr(T),steelstr(steel),wallstr(wall)]);

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'silo_minMassDcOverV.tiff']);


%% Slenderness over hopper half angle
%Parameters
V=1e3;                              %Silo volume
d_h=1;                              %Hopper outlet diameter
phi_r.val=deg2rad([45,39,35,30])';  %Angle of repose
gamma.val=deg2rad(15);              %Roof angle
T=400+273.15;                       %Storage temperature
steel='S235';                       %Wall material
wall=2;                             %Wall type
particles='Sand';                   %Storage material

%Hopper half angle: Split discretization to improve code performance
n=[250,50];
split=10;
beta.val=[linspace(45,split,n(1)),linspace(split,5,n(2)+1)];
beta.val(n(1))=[];
beta.val=deg2rad(beta.val);


%Calculations
[d_cyl,h_cyl]=Silo.minMass(V,d_h,beta.val,phi_r.val,gamma.val,T,steel,wall,particles,2.5e3,[0.6,1.7]);
slenderness=h_cyl./d_cyl;


%Create figure
fig=figure(79);
clf(fig);
ax=gca();

plot(ax,rad2deg(beta.val),slenderness);

grid(ax,'on');

legend(ax,compose('\\phi_r=%.0f°',rad2deg(phi_r.val)),'Location','best');

xlabel(ax,'Hopper half angle \beta (°)');
ylabel(ax,'Slenderness h_{cyl}/d_{cyl} (-)');

Vstr=['V=',num2str(V),'m³, '];
d_hStr=['d_{hop}=',num2str(d_h),'m, '];
title(ax,[Vstr,d_hStr,degstr(gamma),Tstr(T),steelstr(steel),wallstr(wall)]);

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'silo_minMassSlenderOverBeta.tiff']);

save('minMassHopperHalfAngle.mat');


%% Slenderness over temperature for different materials
%Parameters
V=1e3;                                  %Silo volume
d_h=1;                                  %Hopper outlet diameter
beta.val=deg2rad(15);                   %Hopper half angle
phi_r.val=deg2rad(39);                  %Angle of repose
gamma.val=deg2rad(15);                  %Roof angle
T=linspace(100,850,300)+273.15;         %Storage temperature
steel={'S235';'S420';'S600';'S700'};    %Wall material
wall=2;                                 %Wall type
particles='Sand';                       %Storage material


%Calculations
[d_cyl,h_cyl]=Silo.minMass(V,d_h,beta.val,phi_r.val,gamma.val,T,steel,wall,particles,2.5e3,[0.2,1.1]);
slenderness=h_cyl./d_cyl;


%Create figure
fig=figure(80);
clf(fig);
ax=gca();

plot(ax,T-273.15,slenderness);

grid(ax,'on');

legend(ax,compose('Material=%s',string(steel)),'Location','best');

xlabel(ax,'Storage temperature T (°C)');
ylabel(ax,'Slenderness h_{cyl}/d_{cyl} (-)');

Vstr=['V=',num2str(V),'m³, '];
d_hStr=['d_{hop}=',num2str(d_h),'m, '];
title(ax,[Vstr,d_hStr,degstr(beta),degstr(phi_r),degstr(gamma),wallstr(wall)]);

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'silo_minMassSlenderOverTempSteel.tiff']);

save('minMassSlenderOverTempSteel.mat');


%% Slenderness over temperature for different wall types
%Parameters
V=1e3;                              %Silo volume
d_h=1;                              %Hopper outlet diameter
beta.val=deg2rad(15);               %Hopper half angle
phi_r.val=deg2rad(39);              %Angle of repose
gamma.val=deg2rad(15);              %Roof angle
T=linspace(100,850,300)+273.15;     %Storage temperature
steel='S235';                       %Wall material
wall=(1:3)';                        %Wall type
particles='Sand';                   %Storage material


%Calculations
[d_cyl,h_cyl]=Silo.minMass(V,d_h,beta.val,phi_r.val,gamma.val,T,steel,wall,particles,2.5e3,[0.6,1.4]);
slenderness=h_cyl./d_cyl;


%Create figure
fig=figure(81);
clf(fig);
ax=gca();

plot(ax,T-273.15,slenderness);

grid(ax,'on');

legend(ax,compose('Wall=D%.0f',wall),'Location','best');

xlabel(ax,'Storage temperature T (°C)');
ylabel(ax,'Slenderness h_{cyl}/d_{cyl} (-)');

Vstr=['V=',num2str(V),'m³, '];
d_hStr=['d_{hop}=',num2str(d_h),'m, '];
titstr=[Vstr,d_hStr,degstr(beta),degstr(phi_r),degstr(gamma),steelstr(steel)];
title(ax,titstr(1:end-2));

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'silo_minMassSlenderOverTempWall.tiff']);

save('minMassSlenderOverTempWall.mat');




