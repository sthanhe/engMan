%% Fluidization gas distribution box
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the sandTES Engineering Manual
%
%All required files for this class can be found in the software
%repository:
%https://doi.org/10.5281/ZENODO.10207330
% 
%All parameters and results are in SI base units.
%
%
%
%This class contains the functions to design a combination of gas boxes for
%a sandTES heat exchanger. 
%
%
%Requires all files packaged in the class folder and on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%   - Curve Fitting Toolbox, version 3.9
%Additional classes:
%   - DryAir
%   - FluBed
%   - GasBox
%   - Orifice
%   - SiO2
%   - Sinter
%   - implExp


classdef GasBox < handle
    properties
        nDisc=51;       %[1,1] double, number of cells for discretization

        d_p             %[1,1] double, mean particle diameter
        rho_p=2650;     %[1,1] double, raw particle density
        eps=0.5;        %[1,1] double, mean bed porosity

        p0=1e5;         %[1,1] double, ambient pressure

        l               %[1,n] double, length of each chamber
        width           %[1,1] double, chamber width
        h0              %[1,1] double, basic bed height

        AC              %[1,n] boolean, does chamber have an air cushion

        nSupply         %[1,n] double, number of air supply tubes into the gas box, for each chamber
        D               %[1,1] double, inner diameter of air supply tube
        d               %[1,n] double, inner diameter of orifice plate, for each chamber

        grade           %[1,1] double, filter grade of porous plate
        s               %[1,n] double, thickness of porous filter, for each chamber

        pSupply         %[1,1] double, common supply pressure into all gas boxes
        p1              %[1,n] double, pressure inside each gas box
        

        %Automatic update
        p2      %depends on h
        mDot    %depends on mDotS
        FG      %depends on mDotS
    end


    properties(Dependent)
        center          %[1,1] double, center index of discretization
        name            %[1,1] char, Porous plate name
        Tsupply         %[1,1] double, common supply temperature into all gas boxes
        qm              %[1,n] double, mass flow through each orifice plate


        %With private storage / automatic update
        T       %creates a discretized T-vector automatically (depends on nDisc)
        h       %updates p2
        mDotS   %updates mDot and FG
    end


    properties(Access=private)
        strg_T
        strg_h
        strg_mDotS
    end


    methods
        %Class constructor
        function obj=GasBox(d_p,rho_p,eps,p0,l,width,h0,AC,nSupply,D,...
                            grade)

            if nargin==0
                return
            end
            
            obj.d_p=d_p;
            obj.rho_p=rho_p;
            obj.eps=eps;

            obj.p0=p0;

            obj.l=l;
            obj.width=width;
            obj.h0=h0;

            obj.AC=AC;

            obj.nSupply=nSupply;
            obj.D=D;
            
            obj.grade=grade;            
        end


        %Design the gas boxes based on a certain use case
        function design(obj,FG_ref,FGdiffRel,deltaH,deltaPbase)
            %FG_ref         reference degree of fluidization
            %FGdiffRel      allowable difference in / relative to FG_ref
            %deltaH         Bed level difference
            %deltaPbase     Basic (minimum) pressure loss


            extrap=[obj.l]./max([obj.l]);
            if isscalar(deltaH)
                deltaH=deltaH.*extrap;
            end

            if isscalar(FGdiffRel)
                FGdiffRel=FGdiffRel.*extrap;
            end
            FGdiff=FGdiffRel.*FG_ref;

            
            hvec=cumsum([deltaH(2:end),obj(end).h0],'reverse');
            hvec=[deltaH;zeros(1,length(deltaH))]+hvec;

            FGvec=[-1/2;1/2].*FGdiff+FG_ref;

            for i=1:length(obj)
                obj(i).h=GasBox.expand(hvec(:,i),obj(i).nDisc);
                obj(i).FG=GasBox.expand(FGvec(:,i),obj(i).nDisc);
            end
            

            obj.plateThick;

            
            pSupplyScal=max([obj.p1])+deltaPbase;
            deltaP=pSupplyScal-[obj.p1];
            
            for i=1:length(obj)
                obj(i).pSupply=pSupplyScal;
                obj(i).d=fzero(@(d) ...
                            Orifice.deltaOmega(obj(i).qm,obj(i).pSupply,obj(i).Tsupply,d,obj(i).D)-deltaP(i),...
                        [obj(i).D*0.1,obj(i).D*0.75]);
            end
        end


        %Simulate the distribution of fluidization gas for a different use case
        function sim(obj,FG_ref,contr)
            %FG_ref         reference degree of fluidization
            %contr          %controlled chamber, index from left


            obj(contr).FGdist(FG_ref);
            
            pSupplyScal=fzero(@(p) ...
                    Orifice.deltaOmega(obj(contr).qm,p,obj(contr).Tsupply,obj(contr).d,...
                        obj(contr).D)+obj(contr).p1-p,...
                    [obj(contr).p1,2*obj(contr).p1]);
            
            
            % p2_ref=arrayfun(@(i) obj(i).p2(obj(i).center),1:length(obj));
            % for i=1:length(obj)
            %     obj(i).p1=fzero(@(p1) deltaMdot(obj(i),p1),...
            %         [p2_ref(i)*2,p2_ref(i)]);
            %     [~,~,obj(i).mDotS]=Sinter.w(obj(i).p1,obj(i).p2,obj(i).T,obj(i).s,obj(i).name);
            % end
            % 
            % 
            % function deltaMdot=deltaMdot(obj,p1)
            %     deltaPorif=pSupplyScal-p1;
            %     mDotOrif=fzero(@(mDot)...
            %                     Orifice.deltaOmega(mDot./obj.nSupply,pSupplyScal,...
            %                     obj.Tsupply,obj.d,obj.D)-deltaPorif,...
            %                 [obj.mDot/100,obj.mDot*100]);
            % 
            %     [~,~,mDotSplate]=Sinter.w(p1,obj.p2,obj.T,obj.s,obj.name);
            %     mDotPlate=obj.width.*obj.l./obj.nDisc.*trapz(mDotSplate,1);
            % 
            %     deltaMdot=mDotOrif-mDotPlate;
            % end


            p2_ref=arrayfun(@(i) obj(i).p2(obj(i).center),1:length(obj));
            T2=arrayfun(@(i) obj(i).T(obj(i).center),1:length(obj));
            A=[obj.l].*[obj.width];
            deltaP=pSupplyScal-p2_ref;
            for i=1:length(obj)
                [~]=fzero(@(mDot) deltaPOrifPorous(obj(i),mDot,i)-deltaP(i),...
                    [obj(i).mDot/100,obj(i).mDot*100]);
                [~,~,obj(i).mDotS]=Sinter.w(obj(i).p1,obj(i).p2,obj(i).T,obj(i).s,obj(i).name);
            end


            function deltaP=deltaPOrifPorous(obj,mDot,i)
                deltaPorif=Orifice.deltaOmega(mDot./obj.nSupply,pSupplyScal,...
                    obj.Tsupply,obj.d,obj.D);
                obj.p1=pSupplyScal-deltaPorif;

                pPorous=mean([obj.p1;p2_ref(i)],1);
                w=mDot./(DryAir.rho(pPorous,T2(i)).*A(i));
                deltaPporous=Sinter.deltaP(w,pPorous,T2(i),obj.s,obj.name);

                deltaP=deltaPorif+deltaPporous;
            end
        end


        %Set the required plate thickness based on a certain distribution
        %of fluidization gas
        function plateThick(obj)
            FGmat=[obj.FG];
            FG_ref=arrayfun(@(i) FGmat(obj(i).center,i),1:length(obj));
            FGdiff=FGmat(end,:)-FGmat(1,:);

            for i=1:length(obj)
                obj(i).s=fzero(@(s) ...
                            obj(i).FGdist(FG_ref(i),s)-FGdiff(i),...
                        [1e-6,10]);
            end
        end


        %Calculate the distribution of FG based on the current operating
        %conditions
        function FGdiff=FGdist(obj,FG_ref,s)
            if nargin>2
                obj.s=s;
            end

            p=([obj.p0]+[obj.p2])./2;

            
            %Reference point (ref): bed center
            p_ref=arrayfun(@(i) p(obj(i).center,i),1:length(obj));
            T_ref=arrayfun(@(i) obj(i).T(obj(i).center),1:length(obj));
            wmf_ref=FluBed.wmf([obj.d_p],[obj.rho_p],p_ref,T_ref);
            w_ref=FG_ref.*wmf_ref;
            rho_ref=DryAir.rho(p_ref,T_ref);
            mDotS_ref=w_ref.*rho_ref;

            
            p2_ref=arrayfun(@(i) obj(i).p2(obj(i).center),1:length(obj));
            rho2_ref=DryAir.rho(p2_ref,T_ref);
            w2_ref=mDotS_ref./rho2_ref;

            for i=1:length(obj)
                obj(i).p1=fzero(@(p1) ...
                            Sinter.w(p1,p2_ref(i),T_ref(i),obj(i).s,obj(i).name)-w2_ref(i),...
                        [p2_ref(i),p2_ref(i).*100]);

                [~,~,obj(i).mDotS]=Sinter.w(obj(i).p1,obj(i).p2,obj(i).T,obj(i).s,obj(i).name);
            end
            
            FGmat=[obj.FG];
            FGdiff=FGmat(end,:)-FGmat(1,:);
        end


        %Change the orientation of some vectors based on the operating
        %direction
        function direction(obj,dir)
            %dir        Operating direction, 'charge' or 'discharge'

            hmat=[obj.h];
            hvec=reshape(hmat,1,numel(hmat));
            switch dir
                case 'charge'
                    hvec=sort(hvec,2,'descend');
                case 'discharge'
                    hvec=sort(hvec,2,'ascend');
            end
            hmat=reshape(hvec,size(hmat));

            for i=1:length(obj)
                obj(i).h=hmat(:,i);
            end
        end


        %Create a bar plot of the distribution of pressure losses
        function [fig,ax,deltaP]=plotDeltaP(obj,figidx)
            % deltaPorif=Orifice.deltaOmega([obj.qm],[obj.pSupply],[obj.Tsupply],[obj.d],[obj.D]);
            deltaPorif=[obj.pSupply]-[obj.p1];

            p2_ref=arrayfun(@(i) obj(i).p2(obj(i).center),1:length(obj));
            deltaPplate=[obj.p1]-p2_ref;
            
            hmat=[obj.h];
            hmean=arrayfun(@(i) hmat(obj(i).center,i),1:length(obj));
            deltaH=abs(hmat(end,:)-hmat(1,:));
            deltaHAC=hmean-deltaH./2-[obj.h0];
            deltaHAC(~[obj.AC])=0;
            
            deltaPbed=FluBed.deltaP(hmean-deltaHAC,[obj.eps],[obj.rho_p]);
            deltaPAC=FluBed.deltaP(deltaHAC,[obj.eps],[obj.rho_p]);
            
            deltaP=[deltaPorif',deltaPplate',deltaPbed',deltaPAC'];
            
            
            fig=figure(figidx);
            clf(fig);
            ax=gca();
            
            b=bar(ax,1:length(obj),deltaP.*10^-2,'stacked');
            
            legend(ax,fliplr(b),{'AC','Bed','Porous','Orifice'},'Location','bestoutside');
            % legend(ax,fliplr(b),{'\Deltap_{AC}','\Deltap_{bed}','\Deltap_{plate}','\Deltap_{orifice}'},'Location','bestoutside');
            
            xlabel(ax,'Chamber');
            ylabel(ax,'Pressure loss \Deltap (mbar)');
            
            fig.Units='centimeters';
            fig.Position=[10,5,17,8.5];
        end


        %Plot the distribution of fluidization gas
        function [fig,ax]=plotFG(obj,figidx)
            x=cumsum([0,obj.l]);
            x=cell2mat(...
                    arrayfun(@(i) ...
                        linspace(x(i),x(i+1),obj(i).nDisc),...
                1:length(obj),'UniformOutput',false));

            y=reshape([obj.FG],1,numel([obj.FG]));


            fig=figure(figidx);
            clf(fig);
            ax=gca();

            plot(ax,x,y);

            xlabel(ax,'Distance from the left x (m)');
            ylabel(ax,'Degree of fluidization FG (-)');
            
            fig.Units='centimeters';
            fig.Position=[10,5,17,8.5];
        end
    end


    %Dependent variables (no private storage)
    methods
        function val=get.center(obj)
            val=round(mean([1,obj.nDisc]));
        end

        function val=get.name(obj)
            val=compose('SIKA-R %.0f AX',obj.grade);
        end

        function val=get.Tsupply(obj)
            val=mean(obj.T,'all');
        end

        function val=get.qm(obj)
            val=obj.mDot./obj.nSupply;
        end
    end


    %Automatic updates (with private storage)
    methods
        function set.T(obj,val)
            if size(val,1)==obj.nDisc
                obj.strg_T=val;
            else
                obj.strg_T=GasBox.expand(val,obj.nDisc);
            end
        end

        function val=get.T(obj)
            val=obj.strg_T;
        end


        function set.mDotS(obj,val)
            obj.strg_mDotS=val;
            obj.mDot=obj.width.*obj.l./obj.nDisc.*trapz(val,1);

            p=(obj.p0+obj.p2)./2;
            w=val./DryAir.rho(p,obj.T);
            wmf=FluBed.wmf(obj.d_p,obj.rho_p,p,obj.T);
            obj.FG=w./wmf;
        end

        function val=get.mDotS(obj)
            val=obj.strg_mDotS;
        end


        function set.h(obj,val)
            obj.strg_h=val;
            obj.p2=obj.p0+FluBed.deltaP(val,obj.eps,obj.rho_p);
        end

        function val=get.h(obj)
            val=obj.strg_h;
        end
    end


    %Auxiliary methods
    methods(Static,Access=protected)
        function x=expand(x,nDisc)
            x=cell2mat(...
                arrayfun(@(i)...
                    linspace(x(1,i),x(end,i),nDisc)',...
                1:size(x,2),'UniformOutput',false));
        end
    end
end




