%% Basic heat exchanger sizing
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
%This script performs the basic heat exchanger sizing for the example given
%in Section 7.7 of the Engineering Manual.
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
%   - idealDpFinned.m


%% Calculate ideal particle diameter, load basic parameters
idealDpFinned;


%Derived parameters
rhoH2O=IF97.v(pH2O,TH2O).^-1;   %Water density


%% Parameters
p_c=pH2O.*1.1;  %Design pressure

R_m=450e6;  %Tensile strength of tube material
R_p=146e6;  %Proof strength of tube material
R_mT=84e6;  %Creep rupture strength at design temperature and 200 000 h of operation of tube material

sRound=0.5e-3;  %Wall thickness to round up to

G=300;  %Mass flux density

F=0.85;     %Flow factor


%% Geometry calculations
f=min([R_m/2.4,R_p/1.5,R_mT/1.25]);     %Yield strength

s_min=p_c.*d./(2*f+p_c);        %Minimal wall thickness
s=ceil(s_min./sRound).*sRound;  %Chosen wall thickness

d_i=d-2.*s;                 %Inner tube diameter
A=d_i.^2.*pi./4;            %Inner tube cross section

n_tubes=round(mDotH2O./(G.*A));     %Number of tubes


%% Required tube length
%Start values for iteration (required since inner HTC depends on tube length)
l_tube=1;   %Tube length
err=1;      %Iteration error
count=0;    %Iteration counter


%Iteration
while err>1e-3 && count<10
    %Inner heat transfer coefficient (Gnielinski)
    alpha_i=alpha_in(rhoH2O,TH2O,mDotH2O,d_i,n_tubes,linspace(0,l_tube,n));
    
    %Thermal transmittance
    k=(1./alpha_o+d./(2.*lambda).*log(d./d_i)+d./(alpha_i.*d_i)).^-1;
    
    A_effReq=Qdot./(F.*mean(k.*DeltaT));    %Required effective heat transfer area
    l_tubeNew=A_effReq./(A_spec.*n_tubes);  %Required tube length

    err=l_tubeNew-l_tube;   %Check iteration error
    l_tube=l_tubeNew;       %Set neq tube length for next iteration

    count=count+1;  %Abort loop when iteration diverges
end
clear('l_tubeNew');




