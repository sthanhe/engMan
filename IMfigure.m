%% Integral Mean (IM) figure creation
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the sandTES Engineering Manual
%
%All required files for this function can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.xxxxxxx
% 
%All parameters and results are in SI base units.
%
%
%
%This function creates the figures that show the ideal particle diameter as
%the point where the integral mean (IM) reaches its maximum. It is used by
%the scripts idealDp and idealDpFinned.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%Necessary files, classes, functions, and scripts:
%   - None


function fig=IMfigure(d_p,IM,d_pIdeal,figidx)
    %d_p        [1,n] Particle diameter
    %IM         [1,n] Integral mean
    %d_pIdeal   [1,1] Ideal particle diameter
    %figidx     [1,1] Figure index


    fig=figure(figidx);
    clf(fig);
    ax=gca();
    hold(ax,'on');
    
    plot(ax,d_p.*1e6,IM);
    
    xline(d_pIdeal.*1e6);
    text(ax,d_pIdeal.*1e6,max(IM)*3/4,[' ',num2str(round(d_pIdeal.*1e6)),' µm']);
    
    hold(ax,'off');
    
    xlabel(ax,'Particle diameter (µm)');
    ylabel(ax,'IM/$\dot Q$ (1/m)','Interpreter','latex');
    
    fig.Units='centimeters';
    fig.Position=[10,5,17,8.5];
end




