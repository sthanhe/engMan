%% Silos with multiple hoppers
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
%This script calculates the dimensions of silos with multiple hoppers for 
%given loads and sizes. It creates Figures 83-85 and calculates the 
%dimensions shown in Figure 86 of Section 11.4 of the Engineering Manual.
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


%Function handle to generate strings for graphs
degstr=@(ang) ['\',ang.name,'=',num2str(rad2deg(ang.val)),'°, '];   %Degrees


%% Surface over number of hoppers
%Parameters
V=1e3;                  %Silo volume
d_h=1;                  %Hopper outlet diameter
beta.val=deg2rad(15);   %Hopper half angle
phi_r.val=deg2rad(39);  %Angle of repose
gamma.val=deg2rad(15);  %Roof angle
nHop=1:25;              %Number of hoppers


%Calculations
Asingle=Silo.minSurf(V,d_h,beta.val,phi_r.val,gamma.val,1);     %Singlo silo
bat=Silo.minSurfBat(V,d_h,beta.val,phi_r.val,gamma.val,nHop);   %Battery


%Prime number of silos
nHopPrime=[1,primes(max(nHop)*1.5)];
batPrime=Silo.minSurfBat(V,d_h,beta.val,phi_r.val,gamma.val,nHopPrime);


%Square number of silos
nHopSquare=(1:ceil(sqrt(max(nHop)))).^2;
batSquare=Silo.minSurfBat(V,d_h,beta.val,phi_r.val,gamma.val,nHopSquare);


%Make sure that nHop=1 is identical to a single silo
bat.A_tot(1)=Asingle.A_eff; 
batPrime.A_tot(1)=Asingle.A_eff;
batSquare.A_tot(1)=Asingle.A_eff;


%Create figure
fig=figure(83);
clf(fig);
ax=gca();
hold(ax,'on');

scatter(nHop,bat.A_tot./Asingle.A_eff);
plot(nHopPrime,batPrime.A_tot./Asingle.A_eff)
plot(nHopSquare,batSquare.A_tot./Asingle.A_eff)

hold(ax,'off');

grid(ax,'on');


ax.XLim=[min(nHop),max(nHop)];
legend(ax,{'A_{bat} / A_{single} (n_{hop})','n_{hop} \subset P','n_{hop} \subset N^2'},'Location','northwest');

xlabel(ax,'Number of hoppers n_{hop} (-)');
ylabel(ax,'A_{bat} / A_{single} (-)');

Vstr=['V=',num2str(V),'m³, '];
d_hStr=['d_{hop}=',num2str(d_h),'m, '];
titstr=[Vstr,d_hStr,degstr(beta),degstr(phi_r)];
title(ax,titstr(1:end-2));

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'silo_multHopSurfOverNhop.tiff']);


%% Slenderness over hopper half angle
%Parameters
V=1e3;                                  %Silo volume
d_h=1;                                  %Hopper outlet diameter
beta.val=deg2rad(linspace(45,5,1e2));   %Hopper half angle
phi_r.val=deg2rad([45,39,35,30])';      %Angle of repose
gamma.val=deg2rad(15);                  %Roof angle
nHop=4;                                 %Number of hoppers


%Calculations
[~,d_cyl,h_cyl]=Silo.minSurfBat(V,d_h,beta.val,phi_r.val,gamma.val,nHop);
slenderness=h_cyl./d_cyl;


%Create figure
fig=figure(84);
clf(fig);
ax=gca();

plot(ax,rad2deg(beta.val),slenderness);

grid(ax,'on');

legend(ax,compose('\\phi_r=%.0f°',rad2deg(phi_r.val)),'Location','best');

xlabel(ax,'\beta (°)');
ylabel(ax,'Slenderness h_{cyl}/d_{cyl} (-)');

Vstr=['V=',num2str(V),'m³, '];
d_hStr=['d_{hop}=',num2str(d_h),'m, '];
hopstr=['n_{hop}=',num2str(nHop)];
title(ax,[Vstr,d_hStr,hopstr]);

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'silo_multHopSlenderOverBeta.tiff']);


%% Relative total height compared to 1 hopper
%Parameters
V=1e3;                  %Silo volume
d_h=1;                  %Hopper outlet diameter
beta.val=deg2rad(15);   %Hopper half angle
phi_r.val=deg2rad(39);  %Angle of repose
gamma.val=deg2rad(15);  %Roof angle
nHop=1:25;              %Number of hoppers


%Calculations
%Single silo
[~,~,h_cyl,h_hop,~,~,h_r]=Silo.minSurf(V,d_h,beta.val,phi_r.val,gamma.val,1);
Hsingle=h_hop+h_cyl+h_r;


%Battery
[bat,~,~,bat_hop]=Silo.minSurfBat(V,d_h,beta.val,phi_r.val,gamma.val,nHop);
Hbat=bat.h+bat_hop;


%Prime number of silos
nHopPrime=[1,primes(max(nHop)*1.5)];
[bat,~,~,bat_hop]=Silo.minSurfBat(V,d_h,beta.val,phi_r.val,gamma.val,nHopPrime);
HbatPrime=bat.h+bat_hop;


%Square number of silos
nHopSquare=(1:ceil(sqrt(max(nHop)))).^2;
[bat,~,~,bat_hop]=Silo.minSurfBat(V,d_h,beta.val,phi_r.val,gamma.val,nHopSquare);
HbatSquare=bat.h+bat_hop;


%Make sure that nHop=1 is identical to a single silo
Hbat(1)=Hsingle;
HbatPrime(1)=Hsingle;
HbatSquare(1)=Hsingle;


%Create figure
fig=figure(85);
clf(fig);
ax=gca();
hold(ax,'on');

scatter(nHop,Hbat./Hsingle);
plot(nHopPrime,HbatPrime./Hsingle)
plot(nHopSquare,HbatSquare./Hsingle)

hold(ax,'off');

grid(ax,'on');


ax.XLim=[min(nHop),max(nHop)];
legend(ax,{'H_{bat} / H_{single} (n_{hop})','n_{hop} \subset P','n_{hop} \subset N^2'},'Location','best');

xlabel(ax,'Number of hoppers n_{hop} (-)');
ylabel(ax,'H_{bat} / H_{single} (-)');

titstr=[Vstr,d_hStr,degstr(beta),degstr(phi_r)];
title(ax,titstr(1:end-2));

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'silo_multHopHeightOverNhop.tiff']);


%% Battery dimensions for 1000 m³ of sand (Figure 86)
%Parameters
V=1e3;                  %Silo volume
d_h=1;                  %Hopper outlet diameter
beta.val=deg2rad(15);   %Hopper half angle
phi_r.val=deg2rad(39);  %Angle of repose
gamma.val=deg2rad(15);  %Roof angle
nHop=4;                 %Number of hoppers


%Calculations
[bat,d_cyl,h_cyl,h_hop,h_cylEff,h_heap,h_r]=Silo.minSurfBat(V,d_h,beta.val,phi_r.val,gamma.val,nHop);

Hbat=bat.h+h_hop;                                           %Total battery height
Vhop=h_hop.*pi./12.*(d_cyl.^2+d_cyl.*d_h+d_h.^2).*nHop;     %Hopper volume
VcylEff=d_cyl.^2.*pi./4.*h_cylEff.*nHop;                    %Effective cylinder volume
Vheap=d_cyl.^2.*pi./12.*h_heap.*nHop;                       %Heap volume

L=bat.l.*d_cyl;     %Battery length
W=bat.w.*d_cyl;     %Battery width




