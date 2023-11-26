%% Ideal particle diameter for finned tubes
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
%This script calculates the ideal particle diameter when using finned tubes
%for the  example given in Section 7.3.2 of the Engineering Manual and 
%creates Figure 39.
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
%   - @IF97
%   - @SiO2
%   - @Sinter
%   - @implExp
%   - IMfigure.m
%   - alpha_in.m


%% Parameters
n=1000;     %Number of cells for discretization
    
mDotH2O=0.6;    %Water mass flow
Qdot=1e6;       %Transferred heat

Tsand=[480,100]+273.15;     %Sand temperatures at outlet / inlet
TH2O=[500,80]+273.15;       %Water temperatures at inlet / outlet
pH2O=250e5;                 %Water pressure

p=1e5;  %Bed pressure

d_p=linspace(120e-6,1000e-6,1000)';     %Particle diameters
rho_p=SiO2.rho(300);                    %Particle density

FG=4;           %Degree of fluidization
eps_mf=0.45;    %Porosity at minimum fluidization

d=25e-3;            %Outside tube diameter of tube bundle
h_f=(d-5e-3)/2;     %Fin height
s_f=3e-3;           %Fin thickness
l=1;                %Heat exchanger tube length

lambda=40;          %Thermal conductivity of steel


%Derived parameters
DeltaX=l/n;                         %Cell length in tube direction
D=d+2*h_f;                          %Outside tube diameter including fins
phi=(D./d-1).*(1+0.35*log(D./d));   %Geometry function

pitch=50*d_p+s_f;               %Fin pitch
t=l./pitch;                     %Number of turns
k=pitch./(d*pi);                %Fin slope
x=@(r) 2*r*pi*t.*sqrt(1+k.^2);  %Helix length depending on the reference radius r

A_bottom=s_f.*x(d/2);           %Plain tube area where the fin is attached
A_sides=2*h_f*x(d/2+h_f/2);     %Main fin area on the sides
A_plain=d*pi.*l;                %Plain tube area


%% Temperature difference based on the T-Q graph
QdotVec=linspace(0,Qdot,n);

mDotSand=Qdot/diff(fliplr(SiO2.h(Tsand)));
hSand=SiO2.h(Tsand(1))-QdotVec./mDotSand;
Tsand=SiO2.T_h(hSand);

hH2O=IF97.h(pH2O,TH2O(2))+fliplr(QdotVec)./mDotH2O;
TH2O=IF97.T_ph(pH2O,hH2O);

DeltaT=Tsand-TH2O;


%% Heat transfer areas
%Outside heat transfer coefficient (Molerus)
w=FG.*FluBed.wmf(d_p,rho_p,p,Tsand);    %Fluidization velocity
alpha_o=FluBed.molerus(w,Tsand,Tsand,p,d_p,rho_p,eps_mf);   

%Fin efficiency
X=phi.*d/2.*sqrt(2*alpha_o.total./(lambda*s_f));
eta_f=tanh(X)./X;


A_eff=A_plain-A_bottom+eta_f.*A_sides;  %Effective heat transfer area
A_spec=A_eff/l;                         %Specific heat transfer area


%% Calculate ideal particle diameter
%Get integral mean
qDot=alpha_o.total.*DeltaT.*A_spec;     %Specific heat flux
IM=DeltaX*trapz(qDot,2)./(l.*Qdot);     %Integral mean


%Find optimum
[~,idx]=max(IM);
d_pIdeal=d_p(idx);
alpha_o=alpha_o.total(idx,:);
A_spec=mean(A_spec(idx,:));


%Create figure
fig=IMfigure(d_p,IM,d_pIdeal,39);

exportgraphics(fig,['Figures',filesep,'idealDpFinned.tiff']);




