%% International steam table IF97
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
%This class contains property functions for water. 
%
%
%Requires all files packaged in the class folder and on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14


classdef IF97
    properties(Constant)
        R=0.461526*1000;
        
        T_c=647.096;
        p_c=22.064*10^6;
        rho_c=322;
        s_c=4.41202148223476*10^3;
        
        T_t=273.16;
        p_t=611.657;
    end


    properties(Constant)
        cf = 0.72;    %properties of the fluid                                                    
        q0 = 150000;  %substance-specific values
        d0 = 10^-2;
        Ra0 = 10^-6;
        alpha_0 = 25580;
    end
    
    
    methods(Static)
        createConstants()
        loadPersistent(n)
    end
    
    
    %% Saturation Functions
    methods(Static)
        function p=p_sat(T)
            val=273.15<=T & T<=IF97.T_c;
            
            p=NaN(size(T));
            p(val)=IF97.p_satIntern(T(val));
        end
        
        
        function T=T_sat(p)
            val=6.112126774443449*10^2<=p & p<=IF97.p_c;
            
            T=NaN(size(p));
            T(val)=IF97.T_satIntern(p(val));
        end
        
        
        function rho=rho_satL(T)
            val=273.15<=T & T<=IF97.T_c;
            
            rho=NaN(size(T));
            rho(val)=IF97.rho_satLintern(T(val));
        end
        
        
        function rho=rho_satV(T)
            val=273.15<=T & T<=IF97.T_c;
            
            rho=NaN(size(T));
            rho(val)=IF97.rho_satVintern(T(val));
        end
        
        
        function x=x(rho,T)
            sz=IF97.implExpSize(rho,T);

            val=273.15<=T & T<=IF97.T_c;
            
            rho0=NaN(size(T));
            rho0(val)=IF97.rho_satLintern(T(val));
            
            rho1=NaN(size(T));
            rho1(val)=IF97.rho_satVintern(T(val));
            
            val=val & rho1<=rho & rho<=rho0;

            [rho,rho0,rho1]=IF97.normalize(sz,rho,rho0,rho1);
            
            x=NaN(sz);
            x(val)=IF97.x_intern(rho(val),rho0(val),rho1(val));
        end
        
        
        function is2phase=is2phase(rho,T)
            val=273.15<=T & T<=IF97.T_c;
            sz=size(T);
            
            rho0=NaN(sz);
            rho0(val)=IF97.rho_satLintern(T(val));
            
            rho1=NaN(sz);
            rho1(val)=IF97.rho_satVintern(T(val));
            
            is2phase=val & rho1<rho & rho<rho0;
        end
    end
    
    
    %% Melting and Sublimation Lines
    methods(Static)
        function [p_low,p_upp]=p_melt(T,region)
            persistent a b T_t p_t
            if isempty(a)
                folder=IF97.getFolder();
                vars=load([folder,filesep,'MeltSublConstants.mat'],'a','b','T_t','p_t');
                a=vars.a;
                b=vars.b;
                T_t=vars.T_t;
                p_t=vars.p_t;
            end
            
            
            rangeIh=T_t(1)<=T & T<=IF97.T_t;
            rangeIII=T_t(1)<=T & T<=T_t(2);
            rangeV=T_t(2)<T & T<=T_t(3);
            rangeVI=T_t(3)<T & T<=T_t(4);
            rangeVII=T_t(4)<T & T<=715;
            
            
            if nargin<2
                sz=size(T);
                T=reshape(T,1,numel(T));
                rangeUnamb=(rangeV & ~rangeIh) | rangeVI | rangeVII;
                
                p_low=NaN(sz);
                p_upp=NaN(sz);
                
                if nargout>1
                    p_upp(rangeIII)=IF97.meltfx1(T(rangeIII),a.III,b.III,p_t(1),T_t(1));
                end
                p_upp(rangeV)=IF97.meltfx1(T(rangeV),a.V,b.V,p_t(2),T_t(2));
                p_upp(rangeVI)=IF97.meltfx1(T(rangeVI),a.VI,b.VI,p_t(3),T_t(3));
                p_upp(rangeVII)=IF97.meltfx2(T(rangeVII),a.VII,b.VII,p_t(4),T_t(4));
                
                p_low(rangeIh)=IF97.meltfx1(T(rangeIh),a.Ih,b.Ih,IF97.p_t,IF97.T_t);
                p_low(rangeUnamb)=p_upp(rangeUnamb);
                
            else
                if ischar(region)
                    region={region};
                end
                sz=IF97.implExpSize(T,region);
                T=IF97.normalize(sz,T);
                
                
                regionmx=false([sz,5]);
                page=IF97.pagefx(sz);
                
                regionmx(page(1))=rangeIh & getreg('Ih');
                regionmx(page(2))=rangeIII & getreg('III');
                regionmx(page(3))=rangeV & getreg('V');
                regionmx(page(4))=rangeVI & getreg('VI');
                regionmx(page(5))=rangeVII & getreg('VII');
                
                reg=@(x) regionmx(page(x));
                
                
                p_low=NaN(sz);
                p_low(reg(1))=IF97.meltfx1(T(reg(1)),a.Ih,b.Ih,IF97.p_t,IF97.T_t);
                p_low(reg(2))=IF97.meltfx1(T(reg(2)),a.III,b.III,p_t(1),T_t(1));
                p_low(reg(3))=IF97.meltfx1(T(reg(3)),a.V,b.V,p_t(2),T_t(2));
                p_low(reg(4))=IF97.meltfx1(T(reg(4)),a.VI,b.VI,p_t(3),T_t(3));
                p_low(reg(5))=IF97.meltfx2(T(reg(5)),a.VII,b.VII,p_t(4),T_t(4));
            end
            
            
            function reg=getreg(str)
                reg=reshape(cellfun(@(x) strcmp(x,str),region),size(region));
            end
        end
        
        
        function T=T_melt(p)
            persistent a b T_t p_t pmax
            if isempty(a)
                folder=IF97.getFolder();
                vars=load([folder,filesep,'MeltSublConstants.mat'],'a','b','T_t','p_t');
                a=vars.a;
                b=vars.b;
                T_t=vars.T_t;
                p_t=vars.p_t;
                
                pmax=IF97.p_melt(715);
            end
            
            T=NaN(size(p));
            p=reshape(p,1,numel(p));
            
            
            rangeIh=IF97.p_t<=p & p<=p_t(1);
            rangeIII=p_t(1)<p & p<=p_t(2);
            rangeV=p_t(2)<p & p<=p_t(3);
            rangeVI=p_t(3)<p & p<=p_t(4);
            rangeVII=p_t(4)<p & p<=pmax;
            
            
            T(rangeIII)=inv(p(rangeIII),a.III,b.III,p_t(1),T_t(1));
            T(rangeV)=inv(p(rangeV),a.V,b.V,p_t(2),T_t(2));
            T(rangeVI)=inv(p(rangeVI),a.VI,b.VI,p_t(3),T_t(3));
            
            for i=1:find(rangeIh)
                T(i)=fzero(@(T) IF97.meltfx1(T,a.Ih,b.Ih,IF97.p_t,IF97.T_t)-p(i),[T_t(1),IF97.T_t]);
            end
            
            for i=1:find(rangeVII)
                T(i)=fzero(@(T) IF97.meltfx2(T,a.VII,b.VII,p_t(4),T_t(4))-p(i),[T_t(4),715]);
            end
            
            
            function T=inv(p,a,b,p_t,T_t)
                T=T_t.*(-((p./p_t-1)./a-1)).^(b.^-1);
            end
        end
        
        
        function p=p_subl(T)
            persistent a b
            if isempty(a)
                folder=IF97.getFolder();
                vars=load([folder,filesep,'MeltSublConstants.mat'],'aSub','bSub');
                a=vars.aSub;
                b=vars.bSub;
            end
            
            p=NaN(size(T));
            T=reshape(T,1,numel(T));
            
            val=50<=T & T<=IF97.T_t;
            theta=T(val)./IF97.T_t;
            
            p(val)=IF97.p_t.*exp(theta.^-1.*sum(a.*theta.^b,1));
        end
    end
    
    
    %% rho / p Transformation Functions
    methods(Static)
        function v=v(p,T,x)
            if nargin<3 || isempty(x)
                x=NaN;
            end
            
            [p,T,sz]=IF97.pTcheck(p,T,x);
            [regions,x0,x1]=IF97.pTregions(p,T,x,sz,false);
            [p,T,x,x0,x1]=IF97.normalize(sz,p,T,x,x0,x1);
            [p,T]=IF97.pT2norm(p,T);
            
            gammaregion=regions(1) | regions(2) | regions(5);
            
            [pi,tau]=IF97.pT2pitau(p,T,regions);
            
            v=NaN(sz);
            
            gamma=@(dpi,dtau) IF97.gammaSwitch(pi,tau,dpi,dtau,regions);
            v(gammaregion)=IF97.v_pTx(p(gammaregion),T(gammaregion),pi(gammaregion),gamma);
            
            v(regions(3))=IF97.v3(p(regions(3)),T(regions(3)));
            
            
            v(x0)=IF97.rho_satLintern(T(x0)).^-1;
            v(x1)=IF97.rho_satVintern(T(x1)).^-1;

            sat=NaN(2,nnz(regions(4)));
            sat(1,:)=IF97.rho_satLintern(T(regions(4))).^-1;
            sat(2,:)=IF97.rho_satVintern(T(regions(4))).^-1;

            v(regions(4))=x(regions(4)).*sat(2,:)+(1-x(regions(4))).*sat(1,:);
        end
        
        
        function p=p(rho,T)
            sz=IF97.implExpSize(rho,T);
            [regions,region2]=IF97.rhoTregions(rho,T,sz,false);
            [rho,T]=IF97.normalize(sz,rho,T);
            
            gammaregion=regions(1) | regions(2) | regions(5);
            
            
            p=NaN(sz);
            p(gammaregion)=IF97.p_gamma(rho,T,regions,region2,gammaregion);
            
            delta=IF97.delta3(rho(regions(3)));
            tau=IF97.tau3(T(regions(3)));
            p(regions(3))=rho(regions(3)).*IF97.R.*T(regions(3)).*delta.*IF97.phi3(delta,tau,1,0);
            
            p(regions(4))=IF97.p_satIntern(T(regions(4)));
        end
    end
    
    
    %% prop=prop(p,T,x) Functions
    %  v=v(p,T,x) is in the Transformation Functions section
    methods(Static)
        function u=u(p,T,x)
            if nargin<3 || isempty(x)
                x=NaN;
            end
            
            [p,T,sz]=IF97.pTcheck(p,T,x);
            [regions,x0,x1]=IF97.pTregions(p,T,x,sz,true);
            [p,T,x]=IF97.normalize(sz,p,T,x);
            [p,T]=IF97.pT2norm(p,T);
            
            gammaregion=regions(1) | regions(2) | regions(5);
            [rho,rho0,rho1]=IF97.pTx2rhoT(p,T,regions,x0,x1);
            
            u=NaN(sz);
            u(gammaregion | regions(4))=IF97.u_pTx(p,T,x,regions,gammaregion);
            u(regions(3) | regions(6))=IF97.u_rhoTintern(rho,T,regions,rho0,rho1,x);
        end
        
        
        function s=s(p,T,x)
            if nargin<3 || isempty(x)
                x=NaN;
            end
            
            [p,T,sz]=IF97.pTcheck(p,T,x);
            [regions,x0,x1]=IF97.pTregions(p,T,x,sz,true);
            [p,T,x]=IF97.normalize(sz,p,T,x);
            [p,T]=IF97.pT2norm(p,T);
            
            gammaregion=regions(1) | regions(2) | regions(5);
            [rho,rho0,rho1]=IF97.pTx2rhoT(p,T,regions,x0,x1);
            
            s=NaN(sz);
            s(gammaregion | regions(4))=IF97.s_pTx(p,T,x,regions,gammaregion);
            s(regions(3) | regions(6))=IF97.s_rhoTintern(rho,T,regions,rho0,rho1,x);
        end
        
        
        function h=h(p,T,x)
            if nargin<3 || isempty(x)
                x=NaN;
            end
            
            [p,T,sz]=IF97.pTcheck(p,T,x);
            [regions,x0,x1]=IF97.pTregions(p,T,x,sz,true);
            [p,T,x]=IF97.normalize(sz,p,T,x);
            [p,T]=IF97.pT2norm(p,T);
            
            gammaregion=regions(1) | regions(2) | regions(5);
            [rho,rho0,rho1]=IF97.pTx2rhoT(p,T,regions,x0,x1);
            
            h=NaN(sz);
            h(gammaregion | regions(4))=IF97.h_pTx(p,T,x,regions,gammaregion);
            h(regions(3) | regions(6))=IF97.h_rhoTintern(rho,T,regions,rho0,rho1,x);
        end
        
        
        function c_p=c_p(p,T,x)
            if nargin<3 || isempty(x)
                x=NaN;
            end
            
            [p,T,sz]=IF97.pTcheck(p,T,x);
            [regions,x0,x1]=IF97.pTregions(p,T,x,sz,true);
            [p,T]=IF97.normalize(sz,p,T);
            [p,T]=IF97.pT2norm(p,T);
            
            gammaregion=regions(1) | regions(2) | regions(5);
            rho=IF97.pTx2rhoT(p,T,regions,x0,x1);
            
            c_p=NaN(sz);
            c_p(gammaregion)=IF97.c_p_pTx(p,T,regions,gammaregion);
            c_p(regions(3))=IF97.c_p_rhoTintern(rho,T,regions);
        end
        
        
        function c_v=c_v(p,T,x)
            if nargin<3 || isempty(x)
                x=NaN;
            end
            
            [p,T,sz]=IF97.pTcheck(p,T,x);
            [regions,x0,x1]=IF97.pTregions(p,T,x,sz,true);
            [p,T]=IF97.normalize(sz,p,T);
            [p,T]=IF97.pT2norm(p,T);
            
            gammaregion=regions(1) | regions(2) | regions(5);
            rho=IF97.pTx2rhoT(p,T,regions,x0,x1);
            
            c_v=NaN(sz);
            c_v(gammaregion)=IF97.c_v_pTx(p,T,regions,gammaregion);
            c_v(regions(3))=IF97.c_v_rhoTintern(rho,T,regions);
        end
        
        
        function w=w(p,T,x)
            if nargin<3 || isempty(x)
                x=NaN;
            end
            
            [p,T,sz]=IF97.pTcheck(p,T,x);
            [regions,x0,x1]=IF97.pTregions(p,T,x,sz,true);
            [p,T]=IF97.normalize(sz,p,T);
            [p,T]=IF97.pT2norm(p,T);
            
            gammaregion=regions(1) | regions(2) | regions(5);
            rho=IF97.pTx2rhoT(p,T,regions,x0,x1);
            
            w=NaN(sz);
            w(gammaregion)=IF97.w_pTx(p,T,regions,gammaregion);
            w(regions(3))=IF97.w_rhoTintern(rho,T,regions);
        end
    end
    
    
    %% prop=prop(rho,T) Functions
    %  p=p(rho,T) is in the Transformation Functions section
    methods(Static)
        function u=u_rhoT(rho,T)
            sz=IF97.implExpSize(rho,T);
            [regions,region2,rho0,rho1]=IF97.rhoTregions(rho,T,sz,true);
            [rho,T,rho0,rho1]=IF97.normalize(sz,rho,T,rho0,rho1);
            
            gammaregion=regions(1) | regions(2) | regions(5);
            [p,x]=IF97.rhoT2pTx(rho,T,regions,region2,gammaregion,rho0,rho1);
            
            u=NaN(sz);
            u(gammaregion | regions(4))=IF97.u_pTx(p,T,x,regions,gammaregion);
            u(regions(3) | regions(6))=IF97.u_rhoTintern(rho,T,regions,rho0,rho1,x);
        end
        
        
        function s=s_rhoT(rho,T)
            sz=IF97.implExpSize(rho,T);
            [regions,region2,rho0,rho1]=IF97.rhoTregions(rho,T,sz,true);
            [rho,T,rho0,rho1]=IF97.normalize(sz,rho,T,rho0,rho1);
            
            gammaregion=regions(1) | regions(2) | regions(5);
            [p,x]=IF97.rhoT2pTx(rho,T,regions,region2,gammaregion,rho0,rho1);
            
            s=NaN(sz);
            s(gammaregion | regions(4))=IF97.s_pTx(p,T,x,regions,gammaregion);
            s(regions(3) | regions(6))=IF97.s_rhoTintern(rho,T,regions,rho0,rho1,x);
        end
        
        
        function h=h_rhoT(rho,T)
            sz=IF97.implExpSize(rho,T);
            [regions,region2,rho0,rho1]=IF97.rhoTregions(rho,T,sz,true);
            [rho,T,rho0,rho1]=IF97.normalize(sz,rho,T,rho0,rho1);
            
            gammaregion=regions(1) | regions(2) | regions(5);
            [p,x]=IF97.rhoT2pTx(rho,T,regions,region2,gammaregion,rho0,rho1);
            
            h=NaN(sz);
            h(gammaregion | regions(4))=IF97.h_pTx(p,T,x,regions,gammaregion);
            h(regions(3) | regions(6))=IF97.h_rhoTintern(rho,T,regions,rho0,rho1,x);
        end
        
        
        function c_p=c_p_rhoT(rho,T)
            sz=IF97.implExpSize(rho,T);
            [regions,region2]=IF97.rhoTregions(rho,T,sz,true);
            [rho,T]=IF97.normalize(sz,rho,T);
            
            gammaregion=regions(1) | regions(2) | regions(5);
            p=IF97.rhoT2pTx(rho,T,regions,region2,gammaregion);
            
            c_p=NaN(sz);
            c_p(gammaregion)=IF97.c_p_pTx(p,T,regions,gammaregion);
            c_p(regions(3))=IF97.c_p_rhoTintern(rho,T,regions);
        end
        
        
        function c_v=c_v_rhoT(rho,T)
            sz=IF97.implExpSize(rho,T);
            [regions,region2]=IF97.rhoTregions(rho,T,sz,true);
            [rho,T]=IF97.normalize(sz,rho,T);
            
            gammaregion=regions(1) | regions(2) | regions(5);
            p=IF97.rhoT2pTx(rho,T,regions,region2,gammaregion);
            
            c_v=NaN(sz);
            c_v(gammaregion)=IF97.c_v_pTx(p,T,regions,gammaregion);
            c_v(regions(3))=IF97.c_v_rhoTintern(rho,T,regions);
        end
        
        
        function w=w_rhoT(rho,T)
            sz=IF97.implExpSize(rho,T);
            [regions,region2]=IF97.rhoTregions(rho,T,sz,true);
            [rho,T]=IF97.normalize(sz,rho,T);
            
            gammaregion=regions(1) | regions(2) | regions(5);
            p=IF97.rhoT2pTx(rho,T,regions,region2,gammaregion);
            
            w=NaN(sz);
            w(gammaregion)=IF97.w_pTx(p,T,regions,gammaregion);
            w(regions(3))=IF97.w_rhoTintern(rho,T,regions);
        end
    end
    
    
    %% Backward Equations
    methods(Static)
        T=T_ph(p,h)
        
        
        function T=T_ps(p,s)
            sz=IF97.implExpSize(p,s);
            [regions,region2,region3]=IF97.psregions(p,s,sz);
            [p,s]=IF97.normalize(sz,p,s);

            T=IF97.T_psIntern(p,s,regions,region2,region3);
            T=reshape(T,sz);
        end
        
        
        function p=p_hs(h,s)
            sz=IF97.implExpSize(h,s);
            [regions,region2,region3]=IF97.hsregions(h,s,sz);
            [h,s]=IF97.normalize(sz,h,s);

            p=IF97.p_hsIntern(h,s,regions,region2,region3);
            p=reshape(p,sz);
        end
        
        
        function T=T_hs(h,s)
            sz=IF97.implExpSize(h,s);
            [regions,region2,region3]=IF97.hsregions(h,s,sz);
            [h,s]=IF97.normalize(sz,h,s);

            p=IF97.p_hsIntern(h,s,regions,region2,region3);
            
            T=IF97.T_psIntern(p,s,regions,region2,region3);
            T=reshape(T,sz);
        end
        
        
        function x=x_hs(h,s)
            sz=IF97.implExpSize(h,s);
            regions=IF97.hsregions(h,s,sz);
            [h,s]=IF97.normalize(sz,h,s);

            T=IF97.T_sathsIntern(h,s,regions);
            
            
            gammaregion=regions(4) & T<=623.15;
            reg3=regions(4) & 623.15<T;
            
            
            h0=NaN(sz);
            h1=NaN(sz);
            
            p=IF97.p_satIntern(T(gammaregion));
            
            n=nnz(gammaregion);
            exp=@(x) repmat(x,1,n);
            h0(gammaregion)=IF97.h_pTx(p,T(gammaregion),exp(0),@(x) exp(x==1),true(1,n));
            h1(gammaregion)=IF97.h_pTx(p,T(gammaregion),exp(1),@(x) exp(x==2),true(1,n));
            
            
            rho0=IF97.rho_satLintern(T(reg3));
            rho1=IF97.rho_satVintern(T(reg3));
            
            exp=@(x) repmat(x,1,nnz(reg3));
            h0(reg3)=IF97.h_rhoTintern(rho0,T(reg3),@(x) exp(x==3),rho0,rho1,exp(0));
            h1(reg3)=IF97.h_rhoTintern(rho1,T(reg3),@(x) exp(x==3),rho0,rho1,exp(1));
            

            x=(reshape(h,sz)-h0)./(h1-h0);
        end
        
        
        function T=T_saths(h,s)
            sz=IF97.implExpSize(h,s);
            regions=IF97.hsregions(h,s,sz);
            [h,s]=IF97.normalize(sz,h,s);

            T=IF97.T_sathsIntern(h,s,regions);
            T=reshape(T,sz);
        end
        
        
        function p=p_saths(h,s)
            p=IF97.p_satIntern(IF97.T_saths(h,s));
        end
    end
    
    
    %% Additional Properties
    methods(Static)
        function my=my(rho,T)   %Dynamic Viscosity
            sz=IF97.implExpSize(rho,T);
            val=IF97.checkAddFx(rho,T,sz);
            [rho,T]=IF97.normalize(sz,rho,T);
            
            my=IF97.myintern(rho,T,sz,val);
        end
        
        
        function ny=ny(rho,T)   %Kinematic Viscosity
            ny=IF97.my(rho,T)./rho;
        end
        
        
        function lambda=lambda(rho,T)   %Thermal Conductivity
            sz=IF97.implExpSize(rho,T);
            [val,reg5]=IF97.checkAddFx(rho,T,sz);
            [rho,T]=IF97.normalize(sz,rho,T);
            
            lambda=IF97.lambdaintern(rho,T,sz,val,reg5);
        end
        
        
        function a=a(rho,T) %Thermal Diffusivity
            sz=IF97.implExpSize(rho,T);
            [val,reg5]=IF97.checkAddFx(rho,T,sz);
            [rho,T]=IF97.normalize(sz,rho,T);
            
            [lambda,c_p]=IF97.lambdaintern(rho,T,sz,val,reg5);
            
            a=lambda./reshape(rho,sz)./c_p;
        end
        
        
        function Pr=Pr(rho,T)   %Prandtl Number
            sz=IF97.implExpSize(rho,T);
            [val,reg5]=IF97.checkAddFx(rho,T,sz);
            [rho,T]=IF97.normalize(sz,rho,T);
            
            [lambda,c_p,my]=IF97.lambdaintern(rho,T,sz,val,reg5);
            
            Pr=c_p.*my./lambda;
        end
        
        
        function [kappa,c_p,c_v]=kappa(rho,T)   %Heat Capacity Ratio
            sz=IF97.implExpSize(rho,T);
            [regions,region2]=IF97.rhoTregions(rho,T,sz,true);
            [rho,T]=IF97.normalize(sz,rho,T);
            
            gammaregion=regions(1) | regions(2) | regions(5);
            p=IF97.rhoT2pTx(rho,T,regions,region2,gammaregion);
            
            c_p=NaN(sz);
            c_p(gammaregion)=IF97.c_p_pTx(p,T,regions,gammaregion);
            c_p(regions(3))=IF97.c_p_rhoTintern(rho,T,regions);
            
            c_v=NaN(sz);
            c_v(gammaregion)=IF97.c_v_pTx(p,T,regions,gammaregion);
            c_v(regions(3))=IF97.c_v_rhoTintern(rho,T,regions);
            
            kappa=c_p./c_v;
        end
        
        
        function chi=chi(rho,T) %Isothermal Compressibility
            sz=IF97.implExpSize(rho,T);
            val=~IF97.is2phase(rho,T);
            [rho,T]=IF97.normalize(sz,rho,T);
            
            chi=IF97.chiintern(rho,T,sz,val);
        end
        
        
        function beta=beta(rho,T)   %Isobaric Expansion
            persistent dT
            if isempty(dT)
                dT=0.001;
            end
            
            sz=IF97.implExpSize(rho,T);
            val=~IF97.is2phase(rho,T);
            beta=NaN(sz);
            
            if any(val)
                [rho,T]=IF97.normalize(sz,rho,T);
                
                p=IF97.p(rho(val),T(val));
                drho=diff(IF97.v(p,[T(val);T(val)+dT]).^-1,1,1);
                beta(val)=-drho./dT./rho(val);
            end
        end
        
        
        function sigma=sigma(T) %Surface Tension
            persistent B b my
            if isempty(B)
                B=235.8e-3;
                b=-0.625;
                my=1.256;
            end
            
            val=IF97.T_t<=T & T<=IF97.T_c;
            
            tau=1-T(val)./IF97.T_c;
            sigma=NaN(size(T));
            sigma(val)=B.*tau.^my.*(1+b.*tau);
        end
    end
    
    
    %% Internal: Boundary Functions
    methods(Static, Access=private)
        function boundary=B23(x,asPressure) 
            persistent n p_star
            if isempty(p_star)
                folder=IF97.getFolder();
                nstruct=load([folder,filesep,'BoundaryConstants.mat'],'n_B23');
                n=nstruct.n_B23;
                p_star=10^6;
            end
            
            boundary=NaN(size(x));
            if asPressure
                limit=623.15<=x & x<=863.15;
                x=x(limit);
                boundary(limit)=p_star*(n(1)+n(2)*x+n(3)*x.^2);
            else
                limit=16.5292*10^6<=x & x<=100*10^6;
                x=x(limit)/p_star;
                boundary(limit)=n(4)+sqrt((x-n(5))/n(3));
            end
        end
    end
    
    
    %% Internal: Input Analyzing Functions (p,T,x)
    methods(Static, Access=private)
        function [p,T,sz]=pTcheck(p,T,x)
            if isempty(p)
                p=NaN;
            end
            if isempty(T)
                T=NaN;
            end
            
            sz=IF97.implExpSize(p,T,x);
        end
        
        
        function [regions,x0,x1]=pTregions(p,T,x,sz,sat2regions)
            persistent p_sat623
            if isempty(p_sat623)
                p_sat623=IF97.p_satIntern(623.15);
            end
            
            x0=x==0;
            x1=x==1;
            xexp=true(size(x));
            
            
            page=IF97.pagefx(sz);
            
            if sat2regions
                regionmx=false([sz,6]);
                
                regionmx(page(1))=(273.15<=T & T<=623.15 & IF97.p_sat(T)<=p & p<=100*10^6) | ...
                                (x0 & 273.15<=T & T<=623.15) | ...
                                (x0 & 611.213<=p & p<=p_sat623);
                regionmx(page(2))=(273.15<=T & T<=623.15 & 0<p & p<=IF97.p_sat(T)) | ...
                                (623.15<T & T<=863.15 & 0<p & p<=IF97.B23(T,true)) | ...
                                (863.15<T & T<=1073.15 & 0<p & p<=100*10^6) | ...
                                (x1 & 273.15<=T & T<=623.15) | ...
                                (x1 & 611.213<=p & p<=p_sat623);
                regionmx(page(3))=(623.15<=T & T<=IF97.B23(p,false) & IF97.B23(T,true)<=p & p<=100*10^6) | ...
                                ((x0 | x1) & 623.15<T & T<=IF97.T_c) | ...
                                ((x0 | x1) & p_sat623<p & p<=IF97.p_c);
                regionmx(page(4))=(0<x & x<1 & 273.15<=T & T<=623.15) | ...
                                (0<x & x<1 & 611.213<=p & p<=p_sat623);
                regionmx(page(5))=1073.15<T & T<=2273.15 & 0<p & p<=50*10^6 & xexp;
                regionmx(page(6))=(0<x & x<1 & 623.15<T & T<=IF97.T_c) | ...
                                (0<x & x<1 & p_sat623<p & p<=IF97.p_c);
            else
                regionmx=false([sz,5]);
                
                regionmx(page(1))=273.15<=T & T<=623.15 & IF97.p_sat(T)<=p & p<=100*10^6 & xexp;
                regionmx(page(2))=(273.15<=T & T<=623.15 & 0<p & p<=IF97.p_sat(T)) | ...
                                (623.15<T & T<=863.15 & 0<p & p<=IF97.B23(T,true)) | ...
                                (863.15<T & T<=1073.15 & 0<p & p<=100*10^6) & xexp;
                regionmx(page(3))=623.15<=T & T<=IF97.B23(p,false) & IF97.B23(T,true)<=p & p<=100*10^6 & xexp;
                regionmx(page(4))=(0<x & x<1 & 273.15<=T & T<=IF97.T_c) | ...
                                (0<x & x<1 & 611.213<=p & p<=IF97.p_c);
                regionmx(page(5))=1073.15<T & T<=2273.15 & 0<p & p<=50*10^6 & xexp;
            end
                            
            
            regions=@(x) regionmx(page(x));
        end
        
        
        function [p,T]=pT2norm(p,T)
            psat=isnan(p);
            Tsat=isnan(T);
            
            p(psat)=IF97.p_sat(T(psat));
            T(Tsat)=IF97.T_sat(p(Tsat));
        end
        
        
        function [pi,tau]=pT2pitau(p,T,regions)
            sz=size(p);
            
            pi=NaN(sz);
            tau=NaN(sz);
            
            [pi(regions(1)),tau(regions(1))]=IF97.pitau1(p(regions(1)),T(regions(1)));
            [pi(regions(2)),tau(regions(2))]=IF97.pitau2(p(regions(2)),T(regions(2)));
            [pi(regions(5)),tau(regions(5))]=IF97.pitau5(p(regions(5)),T(regions(5)));
        end
        
        
        function [pi,tau]=pitau1(p,T)
            pi=p./(16.53*10^6);
            tau=1386./T;
        end
        
        
        function [pi,tau]=pitau2(p,T)
            pi=p./10^6;
            tau=540./T;
        end
        
        
        function [pi,tau]=pitau5(p,T)
            pi=p./10^6;
            tau=1000./T;
        end
    end
    
    
    %% Internal: Input Analyzing Functions (rho,T)
    methods(Static, Access=private)
        function [regions,region2,rho0,rho1]=rhoTregions(rho,T,sz,sat2regions)
            val=273.15<=T & T<=IF97.T_c;
            szT=size(T);
            
            rho0=NaN(szT);
            rho0(val)=IF97.rho_satLintern(T(val));
            
            rho1=NaN(szT);
            rho1(val)=IF97.rho_satVintern(T(val));
            
            
            page=IF97.pagefx(sz);
            
            
            region2mx=false([sz,3]);
            region2mx(page(2))=623.15<T & T<=863.15 & 0<rho & rho<=IF97.v(IF97.B23(T,true),T).^-1;
            region2mx(page(3))=863.15<T & T<=1073.15 & 0<rho & rho<=IF97.v(100*10^6,T).^-1;
            
            
            if sat2regions
                regionmx=false([sz,6]);
                
                regionmx(page(1))=273.15<=T & T<=623.15 & rho0<=rho & rho<=IF97.v(100*10^6,T).^-1;
                regionmx(page(3))=623.15<=T & T<=863.15 & IF97.v(IF97.B23(T,true),T).^-1<=rho & rho<=IF97.v(100*10^6,T).^-1;
                regionmx(page(4))=273.15<=T & T<=623.15 & rho1<rho & rho<rho0;
                regionmx(page(5))=1073.15<=T & T<=2273.15 & 0<rho & rho<=IF97.v(50*10^6,T).^-1;
                regionmx(page(6))=623.15<T & T<=IF97.T_c & rho1<rho & rho<rho0;
                
                regionmx(page(3))=regionmx(page(3)) & ~regionmx(page(6));
                
                region2mx(page(1))=273.15<=T & T<=623.15 & 0<rho & rho<=rho1;
            else
                regionmx=false([sz,5]);
                
                regionmx(page(1))=273.15<=T & T<=623.15 & rho0<rho & rho<=IF97.v(100*10^6,T).^-1;
                regionmx(page(3))=623.15<=T & T<=863.15 & IF97.v(IF97.B23(T,true),T).^-1<=rho & rho<=IF97.v(100*10^6,T).^-1;
                regionmx(page(4))=273.15<=T & T<=IF97.T_c & rho1<=rho & rho<=rho0;
                regionmx(page(5))=1073.15<=T & T<=2273.15 & 0<rho & rho<=IF97.v(50*10^6,T).^-1;
                
                regionmx(page(3))=regionmx(page(3)) & ~regionmx(page(4));
                
                region2mx(page(1))=273.15<=T & T<=623.15 & 0<rho & rho<rho1;
            end
            
            regionmx(page(2))=any(region2mx,ndims(region2mx));
            
            
            regions=@(x) regionmx(page(x));
            region2=@(x) region2mx(page(x));
        end
        
        
        function [delta,tau]=rhoT2deltatau(rho,T,regions)
            sz=size(rho);
            
            delta=IF97.delta3(rho(regions(3)));
            
            reg=regions(3) | regions(6);
            tau=NaN(sz);
            tau(reg)=IF97.tau3(T(reg));
        end
        
        
        function delta=delta3(rho)
            delta=rho./IF97.rho_c;
        end
        
        
        function tau=tau3(T)
            tau=IF97.T_c./T;
        end
    end
    
    
    %% Internal: pTx / rhoT Transformation Functions
    methods(Static, Access=private)
        v=v3(p,T)
        p=p_gamma(rho,T,regions,region2,gammaregion)
        
        
        function [rho,rho0,rho1]=pTx2rhoT(p,T,regions,x0,x1)
            rho=NaN(size(p));
            
            satL=regions(3) & x0;
            satV=regions(3) & x1;
            
            rho(satL)=IF97.rho_satLintern(T(satL));
            rho(satV)=IF97.rho_satVintern(T(satV));
            
            % rho(regions(3) & ~(x0 | x1))=IF97.v3(p(regions(3)),T(regions(3))).^-1;
            ind=regions(3) & ~(x0 | x1);
            rho(ind)=IF97.v3(p(ind),T(ind)).^-1;
            
            
            if nargout>1
                rho0=NaN(size(p));
                rho0(regions(6))=IF97.rho_satLintern(T(regions(6)));

                rho1=NaN(size(p));
                rho1(regions(6))=IF97.rho_satVintern(T(regions(6)));
            end
        end
        
        
        function [p,x]=rhoT2pTx(rho,T,regions,region2,gammaregion,rho0,rho1)
            p=NaN(size(rho));
            p(gammaregion)=IF97.p_gamma(rho,T,regions,region2,gammaregion);
            
            if nargin>5
                p(regions(4))=IF97.p_satIntern(T(regions(4)));
                
                reg=regions(4) | regions(6);
                x=NaN(size(rho));
                x(reg)=IF97.x_intern(rho(reg),rho0(reg),rho1(reg));
            end
        end
    end
    
    
    %% Internal: Saturation Functions
    methods(Static, Access=private)
        function p=p_satIntern(T)
            persistent n p_star
            if isempty(p_star)
                folder=IF97.getFolder();
                nstruct=load([folder,filesep,'Region4Constants.mat'],'n');
                n=nstruct.n;
                p_star=10^6;
            end
            
            theta=T+n(9)./(T-n(10));
            A=theta.^2+n(1)*theta+n(2);
            B=n(3)*theta.^2+n(4)*theta+n(5);
            C=n(6)*theta.^2+n(7)*theta+n(8);
            
            p=p_star*(2*C./(-B+sqrt(B.^2-4*A.*C))).^4;
        end
        
        
        function T=T_satIntern(p)
            persistent n p_star
            if isempty(p_star)
                folder=IF97.getFolder();
                nstruct=load([folder,filesep,'Region4Constants.mat'],'n');
                n=nstruct.n;
                p_star=10^6;
            end
                   
            beta=(p/p_star).^(1/4);
            E=beta.^2+n(3)*beta+n(6);
            F=n(1)*beta.^2+n(4)*beta+n(7);
            G=n(2)*beta.^2+n(5)*beta+n(8);
            D=2*G./(-F-sqrt(F.^2-4*E.*G));
            
            T=(n(10)+D-sqrt((n(10)+D).^2-4*(n(9)+n(10)*D)))/2;
        end
        
        
        function rho=rho_satLintern(T)
            persistent b e
            if isempty(b)
                folder=IF97.getFolder();
                vars=load([folder,filesep,'Region4Constants.mat'],'coeff');
                b=vars.coeff(:,2);
                e=[1/3;2/3;5/3;16/3;43/3;110/3];
            end
            
            tau=reshape(1-T./IF97.T_c,1,numel(T));
            rho=IF97.rho_c*(1+sum(b.*tau.^e,1));
        end
        
        
        function rho=rho_satVintern(T)
            persistent c e
            if isempty(c)
                folder=IF97.getFolder();
                vars=load([folder,filesep,'Region4Constants.mat'],'coeff');
                c=vars.coeff(:,3);
                e=[2/6;4/6;8/6;18/6;37/6;71/6];
            end
            
            tau=reshape(1-T./IF97.T_c,1,numel(T));
            rho=IF97.rho_c*exp(sum(c.*tau.^e,1));
        end
        
        
        function x=x_intern(rho,rho0,rho1)
            x=(rho.^-1-rho0.^-1)./(rho1.^-1-rho0.^-1);
        end
    end
    
    
    %% Internal: Melting and Sublimation Functions
    methods(Static, Access=private)
        function p=meltfx1(T,a,b,p_t,T_t)
            if ~isempty(T)
                p=p_t.*(1+sum(a.*(1-(T./T_t).^b),1));
            else
                p=NaN(1,0);
            end
        end
        
        
        function p=meltfx2(T,a,b,p_t,T_t)
            if ~isempty(T)
                p=p_t.*exp(sum(a.*(1-(T./T_t).^b),1));
            else
                p=NaN(1,0);
            end
        end
    end
    
    
    %% Internal: Additional Functions
    methods(Static, Access=private)
        function [val,reg5]=checkAddFx(rho,T,sz)
            persistent p_m273
            if isempty(p_m273)
                p_m273=IF97.p_melt(273.15);
            end
            
            is2phase=IF97.is2phase(rho,T);
            
            p=IF97.p(rho,T);
            
            check=IF97.p_t<=p & p<=p_m273;
            T_min=NaN(sz);
            T_min(check)=IF97.T_melt(p(check));
            T_min(p_m273<p)=273.15;
            
            % Validity region fitted to regular IF97 validity
            reg5=0<p & p<=50e6 & 1073.15<T & T<=1173.15;            
            val=(0<p & p<IF97.p_t & IF97.T_t<=T & T<=1073.15) | ...
                (p<=100e6 & T_min<=T & T<=1073.15) | ...
                reg5 & ~is2phase;
        end
        
        
        function my=myintern(rho,T,sz,val)
            persistent H0 H1 T_ast rho_ast my_ast
            if isempty(H0)
                folder=IF97.getFolder();
                vars=load([folder,filesep,'ViscosityConstants.mat'],'H0','H1');
                H0=vars.H0;
                H1=vars.H1;
                
                T_ast=647.096;
                rho_ast=322;
                my_ast=1e-6;
            end
            
            if any(val)
                rho=rho(val)./rho_ast;
                T=T(val)./T_ast;
                
                i=(0:3)';
                my0=NaN(sz);
                my0(val)=100.*sqrt(T)./sum(H0./(T.^i),1);
                
                i=(0:5)';
                j=reshape(0:6,1,1,7);
                my1=NaN(sz);
                my1(val)=exp(rho.*sum((T.^-1-1).^i.*sum(H1.*(rho-1).^j,3),1));

                my=my_ast.*my0.*my1;
            else
                my=NaN(1,0);
            end
        end
        
        
        function [lambda,c_p,my]=lambdaintern(rho,T,sz,val,reg5)
            persistent L0 L1 A lambda_ast my_ast R Lambda q_D ny gamma xi_0 Gamma_0 T_R
            if isempty(L0)
                folder=IF97.getFolder();
                vars=load([folder,filesep,'ThermCondConstants.mat'],'L0','L1','A');
                L0=vars.L0;
                L1=vars.L1;
                A=vars.A;
                
                lambda_ast=1e-3;
                my_ast=1e-6;
                R=0.46151805e3;
                
                Lambda=177.8514;
                q_D=(0.4e-9).^-1;
                ny=0.63;
                gamma=1.239;
                xi_0=0.13e-9;
                Gamma_0=0.06;
                T_R=1.5;
            end
            
            if any(val)
                n=nnz(val);
                
                my=IF97.myintern(rho,T,sz,val);
                my_ref=reshape(my(val)./my_ast,1,n);
                
                chi=IF97.chiintern(rho,T,sz,val);
                chi=reshape(chi(val),1,n);
                
                
                rho=rho(val);
                T=T(val);
                
                c_p=NaN(sz);
                [kappainv,c_p(val)]=IF97.kappa(rho,T);
                kappainv=kappainv.^-1;
                c_p_ref=reshape(c_p(val)./R,1,n);
                c_p_ref(c_p_ref<0 | 1e13<c_p_ref)=1e13;
                
                
                rho=rho./IF97.rho_c;
                T=T./IF97.T_c;
                
                zeta=chi.*rho.*IF97.p_c;
                zeta(zeta<0 | 1e13<zeta)=1e13;
                
                
                lambda0=NaN(sz);
                i=(0:4)';
                lambda0(val)=sqrt(T)./sum(L0./(T.^i),1);

                lambda1=NaN(sz);
                i=(0:4)';
                j=reshape(0:5,1,1,6);
                lambda1(val)=exp(rho.*sum((T.^-1-1).^i.*sum(L1.*(rho-1).^j,3),1));
                
                
                lambda2=NaN(sz);
                
                j=NaN(size(rho));
                j(rho<=0.310559006)=1;
                j(0.310559006<rho & rho<=0.776397516)=2;
                j(0.776397516<rho & rho<=1.242236025)=3;
                j(1.242236025<rho & rho<=1.863354037)=4;
                j(1.863354037<rho)=5;
                i=(0:5)';
                zeta_R=1./sum(A(:,j).*rho.^i,1);
                
                DeltaChi=rho.*(zeta-zeta_R.*T_R./T);
                DeltaChi(DeltaChi<0)=0;
                xi=xi_0.*(DeltaChi./Gamma_0).^(ny./gamma);
                y=q_D.*xi;
                
                Z=2./(pi.*y).*((1-kappainv).*atan(y)+kappainv.*y- ...
                    (1-exp(-1./(y.^-1+y.^2./(3*rho.^2)))));
                Z(y<1.2e-7)=0;
                
                
                lambda2(val)=Lambda.*rho.*c_p_ref.*T.*Z./my_ref;
                lambda2(reg5)=0;
                
                
                lambda=lambda_ast.*(lambda0.*lambda1+lambda2);
            else
                lambda=NaN(1,0);
            end
        end
        
        
        function chi=chiintern(rho,T,sz,val)
            persistent drho
            if isempty(drho)
                drho=1e-5;
            end
            
            if any(val)
                rho=rho(val);
                dp=diff(IF97.p([rho;rho+drho],T(val)),1,1);
                
                chi=NaN(sz);
                chi(val)=drho./dp./rho;
            else
                chi=NaN(1,0);
            end
        end
    end
    
    
    %% Internal: Region 1, 2, 4, 5 Property Equations
    methods(Static, Access=private)
        function v=v_pTx(p,T,pi,gamma)
            v=IF97.R.*T./p.*pi.*gamma(1,0);
        end
        
        
        function u=u_pTx(p,T,x,regions,gammaregion)
            [pi,tau]=IF97.pT2pitau(p,T,regions);
            
            u=NaN(size(p));
            
            gamma=@(dpi,dtau) IF97.gammaSwitch(pi,tau,dpi,dtau,regions);
            u(gammaregion)=ufx(pi(gammaregion),tau(gammaregion),gamma,gammaregion);
            
            u(regions(4))=IF97.gamma2phase(p,T,x,@ufx,regions(4));
            
            u=u(~isnan(u));
            
            
            function u=ufx(pi,tau,gamma,reg)
                u=IF97.R.*T(reg).*(tau.*gamma(0,1)-pi.*gamma(1,0));
            end
        end
        
        
        function s=s_pTx(p,T,x,regions,gammaregion)
            [pi,tau]=IF97.pT2pitau(p,T,regions);
            
            s=NaN(size(p));
            
            gamma=@(dpi,dtau) IF97.gammaSwitch(pi,tau,dpi,dtau,regions);
            s(gammaregion)=sfx(pi(gammaregion),tau(gammaregion),gamma,gammaregion);
            
            s(regions(4))=IF97.gamma2phase(p,T,x,@sfx,regions(4));
            
            s=s(~isnan(s));
            
            
            function s=sfx(~,tau,gamma,~)
                s=IF97.R.*(tau.*gamma(0,1)-gamma(0,0));
            end
        end
        
        
        function h=h_pTx(p,T,x,regions,gammaregion)
            [pi,tau]=IF97.pT2pitau(p,T,regions);
            
            h=NaN(size(p));
            
            gamma=@(dpi,dtau) IF97.gammaSwitch(pi,tau,dpi,dtau,regions);
            h(gammaregion)=hfx(pi(gammaregion),tau(gammaregion),gamma,gammaregion);
            
            h(regions(4))=IF97.gamma2phase(p,T,x,@hfx,regions(4));
            
            h=h(~isnan(h));
            
            
            function h=hfx(~,tau,gamma,reg)
                h=IF97.R.*T(reg).*(tau.*gamma(0,1));
            end
        end
        
        
        function c_p=c_p_pTx(p,T,regions,gammaregion)
            [pi,tau]=IF97.pT2pitau(p,T,regions);
            
            gamma=@(dpi,dtau) IF97.gammaSwitch(pi,tau,dpi,dtau,regions);
            c_p=IF97.R.*(-tau(gammaregion).^2.*gamma(0,2));
        end
        
        
        function c_v=c_v_pTx(p,T,regions,gammaregion)
            [pi,tau]=IF97.pT2pitau(p,T,regions);
            
            gamma=@(dpi,dtau) IF97.gammaSwitch(pi,tau,dpi,dtau,regions);
            c_v=IF97.R.*(-tau(gammaregion).^2.*gamma(0,2)+(gamma(1,0)-tau(gammaregion).*gamma(1,1)).^2./gamma(2,0));
        end
        
        
        function w=w_pTx(p,T,regions,gammaregion)
            [pi,tau]=IF97.pT2pitau(p,T,regions);
            
            gamma=@(dpi,dtau) IF97.gammaSwitch(pi,tau,dpi,dtau,regions);
            gamma1_0=gamma(1,0);
            w=sqrt(IF97.R.*T(gammaregion).*gamma1_0.^2./((gamma1_0-tau(gammaregion).*gamma(1,1)).^2./(tau(gammaregion).^2.*gamma(0,2))-gamma(2,0)));
        end
    end
    
    
    %% Internal: Region 3, 6 Property Equations
    methods(Static, Access=private)
        function u=u_rhoTintern(rho,T,regions,rho0,rho1,x)
            [delta,tau]=IF97.rhoT2deltatau(rho,T,regions);
            
            u=NaN(size(rho));
            
            phi=@(ddelta,dtau) IF97.phi3(delta,tau(regions(3)),ddelta,dtau);
            u(regions(3))=ufx(delta,phi,regions(3));
            
            u(regions(6))=IF97.phi2phase(rho0,rho1,tau,x,@ufx,regions(6));
            
            u=u(~isnan(u));
            
            function u=ufx(~,phi,reg)
                u=IF97.R.*T(reg).*tau(reg).*phi(0,1);
            end
        end
        
        
        function s=s_rhoTintern(rho,T,regions,rho0,rho1,x)
            [delta,tau]=IF97.rhoT2deltatau(rho,T,regions);
            
            s=NaN(size(rho));
            
            phi=@(ddelta,dtau) IF97.phi3(delta,tau(regions(3)),ddelta,dtau);
            s(regions(3))=sfx(delta,phi,regions(3));
            
            s(regions(6))=IF97.phi2phase(rho0,rho1,tau,x,@sfx,regions(6));
            
            s=s(~isnan(s));
            
            function s=sfx(~,phi,reg)
                s=IF97.R.*(tau(reg).*phi(0,1)-phi(0,0));
            end
        end
        
        
        function h=h_rhoTintern(rho,T,regions,rho0,rho1,x)
            [delta,tau]=IF97.rhoT2deltatau(rho,T,regions);
            
            h=NaN(size(rho));
            
            phi=@(ddelta,dtau) IF97.phi3(delta,tau(regions(3)),ddelta,dtau);
            h(regions(3))=hfx(delta,phi,regions(3));
            
            h(regions(6))=IF97.phi2phase(rho0,rho1,tau,x,@hfx,regions(6));
            
            h=h(~isnan(h));
            
            function h=hfx(delta,phi,reg)
                h=IF97.R.*T(reg).*(tau(reg).*phi(0,1)+delta.*phi(1,0));
            end
        end
        
        
        function c_p=c_p_rhoTintern(rho,T,regions)
            [delta,tau]=IF97.rhoT2deltatau(rho,T,regions);
            tau=tau(regions(3));
            
            phi=@(ddelta,dtau) IF97.phi3(delta,tau,ddelta,dtau);
            c_p=IF97.R.*(-tau.^2.*phi(0,2)+(delta.*phi(1,0)- ...
                            delta.*tau.*phi(1,1)).^2./(2*delta.*phi(1,0)+delta.^2.*phi(2,0)));
        end
        
        
        function c_v=c_v_rhoTintern(rho,T,regions)
            [delta,tau]=IF97.rhoT2deltatau(rho,T,regions);
            tau=tau(regions(3));
            
            phi=@(ddelta,dtau) IF97.phi3(delta,tau,ddelta,dtau);
            c_v=IF97.R.*-tau.^2.*phi(0,2);
        end
        
        
        function w=w_rhoTintern(rho,T,regions)
            [delta,tau]=IF97.rhoT2deltatau(rho,T,regions);
            tau=tau(regions(3));
            
            phi=@(ddelta,dtau) IF97.phi3(delta,tau,ddelta,dtau);
            phi1_0=phi(1,0);
            w=sqrt(IF97.R.*T(regions(3)).*(2*delta.*phi1_0+delta.^2.*phi(2,0)-(delta.*phi1_0-delta.*tau.*phi(1,1)).^2./(tau.^2.*phi(0,2))));
        end
    end
    
    
    %% Internal: gamma-functions
    methods(Static, Access=private)
        function gamma=gammaSwitch(pi,tau,dpi,dtau,regions)
            gamma=NaN(size(pi));
            gamma(regions(1))=IF97.gamma1(pi(regions(1)),tau(regions(1)),dpi,dtau);
            gamma(regions(2))=IF97.gamma2(pi(regions(2)),tau(regions(2)),dpi,dtau);
            gamma(regions(5))=IF97.gamma5(pi(regions(5)),tau(regions(5)),dpi,dtau);
            gamma=gamma(~isnan(gamma));
        end
        
        
        function gamma=gamma1(pi,tau,dpi,dtau)
            persistent I J n
            if any(isempty(I))
                folder=IF97.getFolder();
                vars=load([folder,filesep,'Region1constants.mat'],'I_pT','J_pT','n_pT');
                I=vars.I_pT;
                J=vars.J_pT;
                n=vars.n_pT;
            end
            
            
            if ~isempty(pi)
                if dpi==0 && dtau==0
                    gamma=sum(n.*(7.1-pi).^I.*(tau-1.222).^J,1);
                elseif dpi==1 && dtau==0
                    gamma=sum(-n.*I.*(7.1-pi).^(I-1).*(tau-1.222).^J,1);
                elseif dpi==2 && dtau==0
                    gamma=sum(n.*I.*(I-1).*(7.1-pi).^(I-2).*(tau-1.222).^J,1);
                elseif dpi==0 && dtau==1
                    gamma=sum(n.*(7.1-pi).^I.*J.*(tau-1.222).^(J-1),1);
                elseif dpi==0 && dtau==2
                    gamma=sum(n.*(7.1-pi).^I.*J.*(J-1).*(tau-1.222).^(J-2),1);
                elseif dpi==1 && dtau==1
                    gamma=sum(-n.*I.*(7.1-pi).^(I-1).*J.*(tau-1.222).^(J-1),1);
                end
            else
                gamma=NaN(1,0);
            end

        end


        function gamma=gamma2(pi,tau,dpi,dtau)
            persistent J_0 n_0 I_r J_r n_r
            if any(isempty(J_0))
                folder=IF97.getFolder();
                vars=load([folder,filesep,'Region2constants.mat'],'J_pT0','n_pT0','I_pTr','J_pTr','n_pTr');
                J_0=vars.J_pT0;
                n_0=vars.n_pT0;
                I_r=vars.I_pTr;
                J_r=vars.J_pTr;
                n_r=vars.n_pTr;
            end


            if ~isempty(pi)
                if dpi==0 && dtau==0
                    gamma_0=log(pi)+sum(n_0.*tau.^J_0,1);
                    gamma_r=sum(n_r.*pi.^I_r.*(tau-0.5).^J_r,1);
                elseif dpi==1 && dtau==0
                    gamma_0=pi.^-1;
                    gamma_r=sum(n_r.*I_r.*pi.^(I_r-1).*(tau-0.5).^J_r,1);
                elseif dpi==2 && dtau==0
                    gamma_0=-pi.^-2;
                    gamma_r=sum(n_r.*I_r.*(I_r-1).*pi.^(I_r-2).*(tau-0.5).^J_r,1);
                elseif dpi==0 && dtau==1
                    gamma_0=sum(n_0.*J_0.*tau.^(J_0-1),1);
                    gamma_r=sum(n_r.*pi.^I_r.*J_r.*(tau-0.5).^(J_r-1),1);
                elseif dpi==0 && dtau==2
                    gamma_0=sum(n_0.*J_0.*(J_0-1).*tau.^(J_0-2),1);
                    gamma_r=sum(n_r.*pi.^I_r.*J_r.*(J_r-1).*(tau-0.5).^(J_r-2),1);
                elseif dpi==1 && dtau==1
                    gamma_0=0;
                    gamma_r=sum(n_r.*I_r.*pi.^(I_r-1).*J_r.*(tau-0.5).^(J_r-1),1);
                end

                gamma=gamma_0+gamma_r;
            else
                gamma=NaN(1,0);
            end

        end


        function gamma=gamma5(pi,tau,dpi,dtau)
            persistent J_0 n_0 I_r J_r n_r
            if any(isempty(J_0))
                folder=IF97.getFolder();
                vars=load([folder,filesep,'Region5constants.mat'],'J_0','n_0','I_r','J_r','n_r');
                J_0=vars.J_0;
                n_0=vars.n_0;
                I_r=vars.I_r;
                J_r=vars.J_r;
                n_r=vars.n_r;
            end


            if ~isempty(pi)
                if dpi==0 && dtau==0
                    gamma_0=log(pi)+sum(n_0.*tau.^J_0,1);
                    gamma_r=sum(n_r.*pi.^I_r.*tau.^J_r,1);
                elseif dpi==1 && dtau==0
                    gamma_0=pi.^-1;
                    gamma_r=sum(n_r.*I_r.*pi.^(I_r-1).*tau.^J_r,1);
                elseif dpi==2 && dtau==0
                    gamma_0=-pi.^-2;
                    gamma_r=sum(n_r.*I_r.*(I_r-1).*pi.^(I_r-2).*tau.^J_r,1);
                elseif dpi==0 && dtau==1
                    gamma_0=sum(n_0.*J_0.*tau.^(J_0-1),1);
                    gamma_r=sum(n_r.*pi.^I_r.*J_r.*tau.^(J_r-1),1);
                elseif dpi==0 && dtau==2
                    gamma_0=sum(n_0.*J_0.*(J_0-1).*tau.^(J_0-2),1);
                    gamma_r=sum(n_r.*pi.^I_r.*J_r.*(J_r-1).*tau.^(J_r-2),1);
                elseif dpi==1 && dtau==1
                    gamma_0=0;
                    gamma_r=sum(n_r.*I_r.*pi.^(I_r-1).*J_r.*tau.^(J_r-1),1);
                end

                gamma=gamma_0+gamma_r;
            else
                gamma=NaN(1,0);
            end

        end
        
        
        function prop=gamma2phase(p,T,x,propfx,reg4)
            sat=NaN(2,nnz(reg4));
            
            [pi,tau]=IF97.pitau1(p(reg4),T(reg4));
            gamma=@(dpi,dtau) IF97.gamma1(pi,tau,dpi,dtau);
            sat(1,:)=propfx(pi,tau,gamma,reg4);
            
            [pi,tau]=IF97.pitau2(p(reg4),T(reg4));
            gamma=@(dpi,dtau) IF97.gamma2(pi,tau,dpi,dtau);
            sat(2,:)=propfx(pi,tau,gamma,reg4);
            
            prop=x(reg4).*sat(2,:)+(1-x(reg4)).*sat(1,:);
        end
    end
    
    
    %% Internal: phi-functions
    methods(Static, Access=private)
        function phi=phi3(delta,tau,ddelta,dtau)
            persistent I J n n1
            if isempty(n1)
                folder=IF97.getFolder();
                vars=load([folder,filesep,'Region3constants.mat'],'I','J','n');
                I=vars.I(2:end);
                J=vars.J(2:end);
                n=vars.n(2:end);
                n1=vars.n(1);
            end


            if ~isempty(delta)
                if ddelta==0 && dtau==0
                    phi=n1*log(delta)+sum(n.*delta.^I.*tau.^J,1);
                elseif ddelta==1 && dtau==0
                    phi=n1*delta.^-1+sum(n.*I.*delta.^(I-1).*tau.^J,1);
                elseif ddelta==2 && dtau==0
                    phi=-n1*delta.^-2+sum(n.*I.*(I-1).*delta.^(I-2).*tau.^J,1);
                elseif ddelta==0 && dtau==1
                    phi=sum(n.*delta.^I.*J.*tau.^(J-1),1);
                elseif ddelta==0 && dtau==2
                    phi=sum(n.*delta.^I.*J.*(J-1).*tau.^(J-2),1);
                elseif ddelta==1 && dtau==1
                    phi=sum(n.*I.*delta.^(I-1).*J.*tau.^(J-1),1);
                end
            else
                phi=NaN(1,0);
            end

        end
        
        
        function prop=phi2phase(rho0,rho1,tau,x,propfx,reg6)
            delta0=IF97.delta3(rho0(reg6));
            delta1=IF97.delta3(rho1(reg6));
            
            sat=NaN(2,nnz(reg6));
            phi=@(ddelta,dtau) IF97.phi3(delta0,tau(reg6),ddelta,dtau);
            sat(1,:)=propfx(delta0,phi,reg6);
            
            phi=@(ddelta,dtau) IF97.phi3(delta1,tau(reg6),ddelta,dtau);
            sat(2,:)=propfx(delta1,phi,reg6);
            
            prop=x(reg6).*sat(2,:)+(1-x(reg6)).*sat(1,:);
        end
    end
    
    
    %% Internal: Backward Equations
    methods(Static, Access=private)
        [regions,region2,region3]=psregions(p,s,sz)
        T=T_psIntern(p,s,regions,region2,region3)
        
        [regions,region2,region3]=hsregions(h,s,sz)
        p=p_hsIntern(h,s,regions,region2,region3)
        T=T_sathsIntern(h,s,regions);
    end
    
        
    %% Auxillary functions
    methods(Static, Access=private)
        function folder=getFolder()
            folder=mfilename('fullpath');
            foldercell=strsplit(folder,filesep);
            folder=strjoin(foldercell(1:end-1),filesep);
        end
        
        
        function padsz=padsz(dims,x)
            padsz=[size(x),ones(1,dims-ndims(x))];
        end
        
        
        function sz=implExpSize(varargin)
            dims=max(cellfun(@(x) ndims(x),varargin));
            szCell=cellfun(@(x) IF97.padsz(dims,x),varargin,'UniformOutput',false);
            sz=cell2mat(szCell');
            
            %all dimensions not equal to 1 must be the same
            not1=sz~=1;
            check=arrayfun(@(x) all(sz(not1(:,x),x)==sz(find(not1(:,x),1),x)),1:size(sz,2));
            
            if any(~check)
                throw(MException('IF97:incompatibleArrays','Arrays have incompatible sizes'));
            else
                sz=max(sz,[],1);
            end
        end
        
        
        function varargout=normalize(sz,varargin)
            n=prod(sz);
            pad=@(x) IF97.padsz(numel(sz),x);
            
            varargout=cellfun(@(x) reshape(repmat(x,sz./pad(x)),1,n),varargin,'UniformOutput',false);
        end
        
        
        function page=pagefx(sz)
            n=prod(sz);
            page=@(x) (x-1)*n+1:x*n;
        end
    end
end




