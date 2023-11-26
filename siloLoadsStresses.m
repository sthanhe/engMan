%% Loads and stresses in a silo
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
%This script calculates the loads and stresses in a silo to demonstrate
%their impact on silo design. It creates the addin for Figure 72 as well as
%Figures 75 and 77 of Section 11.3 of the Engineering Manual.
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


%% Parameters
h_c=5;  %Cylinder height
d_c=1;  %Characteristic cylinder diameter

T=400+273.15;   %Storage temperature

steel='S235';       %Wall material
wall=2;             %Wall type
particles='Sand';   %Stored particles
critProp=2;         %Load case

t_abr=2e-3;     %Wall thickness added against abrasion

d_h=0.2;            %Hopper outlet diameter
beta=deg2rad(15);   %Hopper half angle

rhoH2O=1000;        %Water density
g=9.81;             %Gravitational acceleration


%% Loads in silo
%Cylinder
[~,p,n,coords]=Silo.cylinder(h_c,d_c,T,steel,wall,particles,critProp);

%Hopper
[~,pStruct]=Silo.cylinder(h_c,d_c,T,steel,wall,particles,3);
p_vfCrit=pStruct.vf;
[~,pHop,x]=Silo.hopper(d_c,d_h,beta,p_vfCrit(end,:),T,steel,particles,wall);

%Combine vectors form cylinder and hopper
p_n=[p.he;flipud(pHop.ne)];     %Normal
p_t=[p.we;flipud(pHop.te)];     %Tangential
p_v=[p.vf;flipud(pHop.vf)];     %Vertical

p_H2O=[rhoH2O.*g.*coords.z;NaN(length(x),1)];   %Water pressure for comparison

p=[p_n,p_t,p_v,p_H2O];                              %Total loads
x=[coords.z;x+coords.z(end)-d_h./2./tan(beta)];     %Vertical coordinates


%Create figure
fig=figure(1);
clf(fig);
ax=gca();
hold(ax,'on');

plot(p,x);

legend(ax,{'p_n','p_t','p_v','p_{H2O}'},'Location','best');

ax.YDir='reverse';
ax.Visible='off';

ax.XLim=[0,max(p,[],'all')];
ax.YLim=[0,max(x)];

fig.Units='centimeters';
fig.Position=[10,5,17,8.28];

exportgraphics(fig,['Figures',filesep,'silo_pressures.tiff']);


%% Stress resultants from loads
n_x=squeeze(n.x);               %Axial stress resultants
n_xTheta=squeeze(n.xTheta);     %Shear stress resultants

z_p=squeeze(coords.z_p);    %Patch load position

%Equivalent stress and where it reaches its maximum
eq=sqrt(n_x.^2+n.theta.^2-n_x.*n.theta+3*n_xTheta.^2);
[~,idx]=max(eq(end,:));
z_pMax=z_p(idx);


%Create figure
fig=figure(75);
clf(fig);
ax=gca();
colors=ax.ColorOrder;


hold(ax,'on');

plot(-n_x(:,[1,end]),coords.z,'Color',colors(1,:));
plot(-n_x(:,idx),coords.z,'Color',colors(2,:));

plot(n.theta,coords.z,'Color',colors(3,:));

plot(n_xTheta(:,[1,end]),coords.z,'Color',colors(1,:),'LineStyle','--');
plot(n_xTheta(:,idx),coords.z,'Color',colors(2,:),'LineStyle','--');


legItems=repmat(line(),2,1);
legItems(1)=plot(NaN,NaN,'Color',colors(1,:));
legItems(2)=plot(NaN,NaN,'Color',colors(2,:));
legItems(3)=plot(NaN,NaN,'Color','k','LineStyle','-');
legItems(4)=plot(NaN,NaN,'Color','k','LineStyle','--');
legItems(5)=plot(NaN,NaN,'Color',colors(3,:));

hold(ax,'off');


ax.YDir='reverse';

legend(ax,legItems,{['z_p=',num2str(round(z_p(1),1)),' m | ',...
                    num2str(round(z_p(end),1)),' m'],...
                    ['z_p=',num2str(round(z_pMax,1)),' m'],...
                    '|n_z|',...
                    'n_{z\theta}',...
                    'n_{\theta}'},'Location','best');

xlabel(ax,'Stress resultant (N/m)');
ylabel(ax,'z (m)');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'silo_resultants.tiff']);


%% Limit states
%Setup
t_minLS1=NaN(1,2);
t_minLS1vec=NaN(100,1,2);
p_hMin=NaN(100,1,2);
silfx=@(critProp) Silo.cylinder(h_c,d_c,T,steel,wall,particles,critProp);


%Load cases
[t,p,~,coords]=silfx(1);
t_minLS1(:,1)=t.min;
t_minLS1vec(:,:,1)=t.vec;
p_hMax=p.he;

[t,p,n]=silfx(2);
t_minLS1(:,2)=t.min;
t_minLS1vec(:,:,2)=t.vec;
p_hMin(:,:,1)=p.he;
n_xMax=n.x;

[~,p]=silfx(3);
p_vf=p.vf;
p_hMin(:,:,2)=p.he;


%Cylinder plastic failure (LS1)
t_minLS1=max(t_minLS1,[],2);
t_minLS1vec=max(t_minLS1vec,[],3);
p_hMin=min(p_hMin,[],3);


%Cylinder buckling (LS3)
n_xMax=max(-n_xMax,[],3);
t_minLS3=fzero(@(t) Silo.buckling(d_c,t,p_hMin,p_hMax,n_xMax,T,steel)-t,...
            [t_minLS1/3,t_minLS1*3]);
[~,t_minLS3vec]=Silo.buckling(d_c,t_minLS3,p_hMin,p_hMax,n_xMax,T,steel);


%Hopper
p_vf=p_vf(end,:)';
[t,~,x]=Silo.hopper(d_c,d_h,beta,p_vf,T,steel,particles,wall);

[t_hop,idx]=max([t.minf,t.mine],[],2);
t_hop=flipud(t_hop)+t_abr;


%% Find critical limit states
%Cylinder: LS1 or LS3
[tCylVec,LSidx]=max([t_minLS1vec,t_minLS3vec],[],2);
LS1=LSidx==1;
LS3=LSidx==2;

tCylVec=tCylVec+t_abr;


%Hopper: filling or extracting
idx=flipud(idx);
filling=idx==1;
extr=idx==2;


%Coordinates
x=x+coords.z(end)-d_h./2./tan(beta);
dz=diff(coords.z(1:2));


%Create figure
fig=figure(77);
clf(fig);
ax=gca();
colors=ax.ColorOrder;
hold(ax,'on');


plot(t_minLS1vec+t_abr,coords.z,'Color',colors(1,:),'LineStyle','--');
plot(t_minLS3vec+t_abr,coords.z,'Color',colors(2,:),'LineStyle','--');

plot(tCylVec(LS1),coords.z(LS1),'Color',colors(1,:));
plot(tCylVec(LS3),coords.z(LS3),'Color',colors(2,:));

plot(flipud(t.minf)+t_abr,x,'Color',colors(3,:),'LineStyle','--');
plot(flipud(t.mine)+t_abr,x,'Color',colors(4,:),'LineStyle','--');

plot(t_hop(filling),x(filling),'Color',colors(3,:));
plot(t_hop(extr),x(extr),'Color',colors(4,:));


yline(coords.z(end),'LineStyle',':');


legItems=repmat(line(ax,'Visible','off'),3,1);
legItems(1)=plot(ax,NaN,NaN,'Color',colors(1,:),'LineStyle','--');
legItems(2)=plot(ax,NaN,NaN,'Color',colors(2,:),'LineStyle','--');
legItems(3)=plot(ax,NaN,NaN,'Color',colors(3,:),'LineStyle','--');
legItems(4)=plot(ax,NaN,NaN,'Color',colors(4,:),'LineStyle','--');
legItems(5)=plot(ax,NaN,NaN,'Color','k');
legItems(6)=plot(ax,NaN,NaN,'Color','k','LineStyle',':');

legend(ax,legItems,{'LS1','LS3','Filling','Discharge','Critical','Transition'},'Location','best');


hold(ax,'off');

ax.YDir='reverse';
ax.YLim=[0,max(x)];

xlabel(ax,'Wall thickness (m)');
ylabel(ax,'z (m)');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'silo_critWallThickness.tiff']);




