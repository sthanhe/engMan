%% Minimum required wall thickness of a silo 
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
%This script calculates the required wall thickness of a silo and its 
%hopper for given silo dimensions. A uniform wall thickness is assumed.
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
h_c=37;     %Cylinder height
d_c=18.5;   %Characteristic cylinder diameter

T=100+273.15;   %Storage temperature

steel='S235';       %Wall material
wall=2;             %Wall type
particles='Sand';   %Stored particles

t_abr=2e-3;     %Wall thickness added against abrasion

d_h=2;              %Hopper outlet diameter
beta=deg2rad(30);   %Hopper half angle



%% Calculations
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

[t_minHop,idx]=max([t.minf,t.mine],[],2);


%% Results in millimeters
t_minLS1=(t_minLS1+t_abr).*10^3;
t_minLS3=(t_minLS3+t_abr).*10^3;

t_minCyl=max([t_minLS1,t_minLS3]);
t_minHop=max(t_minHop+t_abr).*10^3;




