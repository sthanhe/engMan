%% Sinter Plate Calculation 
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the sandTES Engineering Manual
%
%All required files for this class can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.xxxxxxx
% 
%All parameters and results are in SI base units.
%
%
%
%This class contains functions to calculate the flow of dry air through a 
%sinter plate
%
%
%Requires all files packaged in the class folder and on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%Data files:
%   - constants.xls
%   - constants.mat
%Additional classes:
%   - DryAir


classdef Sinter
    methods(Static)
        function deltaP=deltaP(w,p,T,s,name)
            %w is the gas velocity at mean pressure inside the filter
            persistent tab
            if isempty(tab)
                tabStruct=coder.load('@Sinter\constants.mat','tab');
                tab=tabStruct.tab;
            end

            if ischar(name)
                name={name};
            end
            
            idx=cellfun(@(x) find(strcmp(tab.Name,x)),name);
            alpha=tab.alpha(idx)';
            beta=tab.beta(idx)';

            eta=DryAir.eta(T);
            rho=DryAir.rho(p,T);

            deltaP=w.*s.*(eta./alpha+rho.*w./beta);
        end


        function [w2,w1,mDotS]=w(p1,p2,T,s,name)
            persistent tab
            if isempty(tab)
                tabStruct=coder.load('@Sinter\constants.mat','tab');
                tab=tabStruct.tab;
            end

            if ischar(name)
                name={name};
            end

            sz=implExp.size(p1,p2,T,s,name);
            [p1,p2,T,s,name]=implExp.normalize(sz,p1,p2,T,s,name);
            
            idx=cellfun(@(x) find(strcmp(tab.Name,x)),name);
            alpha=tab.alpha(idx)';
            beta=tab.beta(idx)';

            p=mean([p1;p2],1);
            eta=DryAir.eta(T);
            rho=DryAir.rho(p,T);

            p=eta.*beta./(rho.*alpha);
            q=-beta.*abs(p1-p2)./(rho.*s);
            w=-p./2+sqrt((p./2).^2-q);
            w=sign(p1-p2).*w;

            mDotS=w.*rho;
            w1=mDotS./DryAir.rho(p1,T);
            w2=mDotS./DryAir.rho(p2,T);

            mDotS=reshape(mDotS,sz);    %in kg/(m²s)
            w1=reshape(w1,sz);
            w2=reshape(w2,sz);
        end


        function [s,deltaP,p1,p2]=s(d_p,rho_p,eps,p0,grade,Tbed,FG,h)
            %T=temperature distribution as a column vector
            %at least: T=[Tleft;Tcenter;Tright]
            %if T has only 2 rows, it's T=[Tleft;Tright] and Tcenter is
            %calculated as the mean between the two
            %same for FG and h
            %multiple columns indicate multiple distributions
            
            reg={d_p,rho_p,eps,p0,grade};
            disc={Tbed,FG,h};

            disc=cellfun(@(x) getMean(x),disc,'UniformOutput',false);

            sz=implExp.sizeDisc(reg,disc,1);
            [reg,disc]=implExp.normDisc(sz,reg,disc,1);

            reg=cell2struct(reg,{'d_p','rho_p','eps','p0','grade'},2);
            disc=cell2struct(disc,{'Tbed','FG','h'},2);


            center=round(mean([1,size(disc.Tbed,1)]));

            FG_ref=disc.FG(center,:);
            FGdiff=disc.FG(end,:)-disc.FG(1,:);            
            
            name=compose('SIKA-R %.0f AX',reg.grade);

            
            s=arrayfun(@(i) ...
                    fzero(@(s) ...
                        FluBed.FGdist(reg.d_p(i),reg.rho_p(i),reg.eps(i),...
                        FG_ref(i),reg.p0(i),name(i),s,disc.Tbed(:,i),disc.h(:,i))-FGdiff(i),...
                    [1e-6,10]),...
                1:length(reg.d_p));

            [~,~,p1,p2]=FluBed.FGdist(reg.d_p,reg.rho_p,reg.eps,FG_ref,reg.p0,name,s,disc.Tbed,disc.h);
            deltaP=p1-p2(center,:);

            s=reshape(s,sz);
            deltaP=reshape(deltaP,sz);


            function x=getMean(x)
                if size(x,1)<3
                    col=repmat({':'},1,ndims(x)-1);
                    x=[x(1,col{:});mean([x(1,col{:});x(end,col{:})],1);x(end,col{:})];
                end
            end

        end


        function loadConstants()
            tab=readtable('@Sinter/constants.xls','Range','A:F');
            tab.alpha=tab.alpha*10^-12;
            tab.beta=tab.beta*10^-7;
            tab.x=tab.x*10^-6;
            tab.DeltaRho=tab.DeltaRho*10^2;
            tab.tau=tab.tau*10^-6;

            save('@Sinter/constants.mat','tab');
        end
    end
end




