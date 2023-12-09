%% Silos
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
%This class contains functions to calculate loads, stresses, and required 
%wall thicknesses of silos based on the various norms that cover it.
%
%
%Requires all files packaged in the class folder and on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%   - Parallel Computing Toolbox, version 7.8
%Additional classes:
%   - implExp


classdef Silo
    properties(Constant)
        gamma_F=1.35;
    end


    methods(Static)
        %Cylindrical part LS1
        function [t,p,n,coords]=cylinder(h_c,d_c,T,steel,wall,particles,critProp,nz,z)
            if nargin<8
                nz=1e2;
            end
            E=0;
        
            sz=implExp.size(h_c,d_c);
            [h_c,d_c]=implExp.normalize(sz,h_c,d_c);
            
            part=Silo.getPartProp(particles,wall,critProp);
            
            % %Optional
            isSlender=false;     %§11, p. 30
            
            
            %% Geometry
            A=d_c.^2*pi/4;
            U=d_c*pi;
            r=d_c/2;
            
            slenderness=h_c./d_c;
            squat=slenderness<=1;
            inter=1<slenderness & slenderness<2;
            slender=2<=slenderness | isSlender;
        
            squatInter=squat | inter;
            h_o=NaN(size(d_c));
            h_o(squatInter)=d_c(squatInter)/6.*tan(part.phi_r);
            s=d_c*pi/16;
        
            if nargin<9 || isscalar(z)
                z=NaN(nz,length(h_c));
                z(:,~squatInter)=cell2mat(arrayfun(@(x) linspace(0,h_c(x),nz)',find(~squatInter),'UniformOutput',false));
                z(:,squatInter)=cell2mat(arrayfun(@(x) linspace(h_o(x),h_c(x),nz)',find(squatInter),'UniformOutput',false));
            end
            z_p=cell2mat(arrayfun(@(x) linspace(z(1,x)+s(x)/2,h_c(x)-s(x)/2,nz)',1:length(h_c),'UniformOutput',false));
            z_p=permute(z_p,[3,2,1]);
            
            
            %% Filling, symmetrical
            z_o=A./(U.*part.K.*part.my);
            z_oSquat=z_o(squatInter);
            p_ho=part.gamma.*part.K.*z_o;
        
        
            if any(slender)
                Y_J=1-exp(-z(:,slender)./z_o(slender));
            else
                Y_J=NaN;
            end
        
            
            if any(squatInter)
                n=-(1+tan(part.phi_r)).*(1-h_o./z_oSquat);
                Y_R=1-((z(:,squatInter)-h_o)./(z_oSquat-h_o)+1).^n;
                Y_R(Y_R<0)=0;
            
                z_v=h_o-1./(n+1).*(z_oSquat-h_o-(z(:,squatInter)+z_oSquat-2*h_o).^(n+1)./(z_oSquat-h_o).^n);
            else
                Y_R=NaN;
                z_v=NaN;
            end
        
        
            p_hf=NaN(size(z));
            p_hf(:,slender)=p_ho(slender).*Y_J;
            p_hf(:,squatInter)=p_ho(squatInter).*Y_R;
        
            p_wf=part.my.*p_hf;
        
            %p_vf only at z=h_c for hopper calculations
            p_vf=NaN(size(z));
            p_vf(:,slender)=p_hf(:,slender)./part.K;
            p_vf(:,squatInter)=part.gamma.*z_v;
        
        
            %% Discharge, symmetrical
            C_S=slenderness(inter)-1;
            
            C_h=NaN(size(slenderness));
            C_h(slender)=1.15;
            C_h(inter)=1+0.15*C_S;
            C_h(squat)=1;
            
            C_w=NaN(size(slenderness));
            C_w(slender)=1.1;
            C_w(inter)=1+0.1*C_S;
            C_w(squat)=1;
            
            p_he=C_h.*p_hf;
            p_we=C_w.*p_wf;
        
        
            %% Discharge, patch load
            C_pe=zeros(size(slenderness));
            idx=slenderness>1.2;
            C_pe(idx)=0.42*part.C_op.*(1+2*E.^2).*(1-exp(-1.5*(slenderness(idx)-1)));
            C_pe(~idx)=0.272*part.C_op*(slenderness(~idx)-1+E);
            C_pe(C_pe<0)=0;
        
            val=~isnan(h_c);
            idx=find(val);
            p_peCell=cell(1,length(d_c));
            p_peCell(~val)=repmat({NaN(1,1,nz)},1,nnz(~val));
            p_peCell(idx)=arrayfun(@(x) interp1(z(:,x),p_he(:,x),z_p(1,x,:)),idx,'UniformOutput',false);
            p_pe=C_pe.*cell2mat(p_peCell);
            p_pe(:,C_pe==0,:)=0;
            
            F_pe=pi/2*s.*d_c.*p_pe;
        
        
            %% Stresses
            idx=z_p-s/2<=z & z<=z_p+s/2;
        
            V=zeros(size(idx));
            M=zeros(size(idx));
        
            s1=z-z_p+s/2;
            Vloc=F_pe./s.*s1;
            V(idx)=Vloc(idx);
        
            Mloc=Vloc.*s1/2;
            M(idx)=Mloc(idx);
            
            idx=z_p+s/2<z;
            Vloc=repmat(F_pe,size(z,1),1,1);
            V(idx)=Vloc(idx);
        
            Mloc=F_pe.*(z-z_p);
            M(idx)=Mloc(idx);
        
        
            p_weInt=cell2mat(arrayfun(@(i) cumtrapz(z(:,i),p_we(:,i),1),1:length(h_c),'UniformOutput',false));
            n_x=-M./(pi*r.^2)-p_weInt;
            n_theta=p_he.*r;
            n_xTheta=V./(pi*r);
        
            t_minVec=Silo.vonMises(n_x,n_theta,n_xTheta,T,steel);
            t_min=reshape(max(t_minVec,[],1),sz);


            t=struct('min',t_min,'vec',t_minVec);
            p=struct('vf',p_vf,'he',p_he,'we',p_we);
            n=struct('x',n_x,'theta',n_theta,'xTheta',n_xTheta);
            coords=struct('z',z,'z_p',z_p);
        end


        %Hopper
        function [t,p,x]=hopper(d_c,d_h,beta,p_vf,T,steel,particles,wall)
            nx=60;
        
            sz=implExp.size(d_c,d_h,beta,p_vf);
            [d_c,d_h,beta,p_vf]=implExp.normalize(sz,d_c,d_h,beta,p_vf);
            
            part=Silo.getPartProp(particles,wall,3);


            %% Geometry
            steep=tan(beta)<(1-part.K)/(2*part.my);
            shallow=~steep;

            r1=d_h./2;
            r2=d_c./2;
            r=cell2mat(arrayfun(@(x) linspace(r1(x),r2(x),nx)',1:length(d_c),'UniformOutput',false));

            s1=r1./sin(beta);
            s2=r2./sin(beta);
            s=cell2mat(arrayfun(@(x) linspace(s1(x),s2(x),nx)',1:length(d_c),'UniformOutput',false));

            h_h=r2./tan(beta);
            h_tip=r1./tan(beta);
            x=cell2mat(arrayfun(@(x) linspace(h_tip(x),h_h(x),nx)',1:length(d_c),'UniformOutput',false));

            %% Loads
            %General
            C_b=1;
            p_vft=C_b*p_vf;


            %Filling
            b=0.2;
            my_heff=NaN(size(d_c));
            my_heff(steep)=part.my;
            my_heff(shallow)=(1-part.K)./(2*tan(beta(shallow)));

            F_f=1-b./(1+tan(beta)./my_heff);
            p_nf=F_f.*p_vFx(F_f,true(size(p_vft)));
            p_tf=my_heff.*p_nf;
            


            %Discharge
            phi_wh=min(atan(part.my),part.phi_i);
            epsilon=phi_wh+asin(sin(phi_wh)/sin(part.phi_i));
            F_e=(1+sin(part.phi_i).*cos(epsilon))./(1-sin(part.phi_i).*cos(2*beta(steep)+epsilon));

            p_ne=NaN(nx,length(d_c));
            p_te=NaN(size(p_ne));

            p_vf=p_vFx(F_e,steep);
            p_ne(:,steep)=F_e.*p_vf;
            p_te(:,steep)=my_heff(steep).*p_ne(:,steep);

            p_ne(:,shallow)=p_nf(:,shallow);
            p_te(:,shallow)=p_tf(:,shallow);


            %% Stresses
            t_minf=tmin(p_nf,p_tf);
            t_mine=tmin(p_ne,p_te);

            t_minVec=max(cat(3,t_minf,t_mine),[],3);
            t_min=max(t_minVec,[],1);
            % p_n=NaN(size(d_c));
            % p_n(idx==1)=p_nf(end,idx==1);
            % p_n(idx==2)=p_ne(end,idx==2);
            % 
            % p_t=NaN(size(d_c));
            % p_t(idx==1)=p_tf(end,idx==1);
            % p_t(idx==2)=p_te(end,idx==2);

            t_min=reshape(t_min,sz);
            % p_n=reshape(p_n,sz);
            % p_t=reshape(p_t,sz);

            t=struct('min',t_min,'vec',t_minVec,'minf',t_minf,'mine',t_mine);
            p=struct('nf',p_nf,'ne',p_ne,'tf',p_tf,'te',p_te,'vf',p_vf);


            %% Auxiliary functions
            function [p_v,n]=p_vFx(F,idx)
                S=2;
                n=S.*(F.*my_heff(idx).*cot(beta(idx))+F)-2;
                f1=x(:,idx)./h_h(idx);
                f2=f1.^n;
                p_v=part.gamma.*h_h(idx)./(n-1).*(f1-f2)+p_vft(idx).*f2;
            end


            function t_min=tmin(p_n,p_t)
                p_phi=cell2mat(arrayfun(@(n) interp1(x(:,n),p_t(:,n),s(:,n).*cos(beta(n)),'linear','extrap'),1:length(d_c),'UniformOutput',false));
                p_phiInt=cell2mat(arrayfun(@(i) cumtrapz(s(:,i),p_phi(:,i).*s(:,i),1),1:length(d_c),'UniformOutput',false));
                n_phi1=p_phiInt./s;
                n_phi1=cell2mat(arrayfun(@(n) interp1(s(:,n),n_phi1(:,n),r(:,n)./sin(beta(n)),'linear','extrap'),1:length(d_c),'UniformOutput',false));
    
                p_nloc=cell2mat(arrayfun(@(n) interp1(x(:,n),p_n(:,n),r(:,n)./tan(beta(n)),'linear','extrap'),1:length(d_c),'UniformOutput',false));
                n_phi2=p_nloc.*r.*((r2./r).^2-1)./(2*cos(beta));
                n_phi=n_phi1+n_phi2;
    
                n_theta=p_nloc.*r./cos(beta);
    
                t_min=Silo.vonMises(n_phi,n_theta,0,T,steel);
            end
        end


        %Cylindrical part LS3
        function [t_min,t_vec]=buckling(d_c,t,p_hMin,p_hMax,n_x,T,steel)
            sz=implExp.size(d_c,t);
            [d_c,t]=implExp.normalize(sz,d_c,t);

            r=(d_c+t)./2;

            ny=0.3;
            [f_pTheta,E_aTheta]=Silo.yieldElast(T,steel);
            sigma_xRcr=E_aTheta*t./(sqrt(3*(1-ny^2))*r);


            %% Non-uniformity parameter
            s=0.5;
            psi_b=0.4;

            j=0.25*sqrt(r./t).*acos(s);
            b1=0.5*sqrt(t./r);
            b2=(1-b1)./psi_b-1;

            idx=j>1./b1;
            j(idx)=1./b1(idx);

            psi=(1-b1.*j)./(1+b2.*j);


            %% Elastic imperfection reduction factor
            Q=25;
            alpha=NaN([size(p_hMin),2]);

            w_ok=t./Q.*sqrt(r./t);
            alpha_0=0.62./(1+1.91*psi.*(w_ok./t).^1.44);


            p_s=pfx(p_hMin);
            alpha(:,:,1)=alpha_0+(1-alpha_0).*p_s./(p_s+0.3./sqrt(alpha_0));


            p_g=pfx(p_hMax);
            s=r./(400*t);
            lambda_x=f_pTheta./sigma_xRcr;
            alpha(:,:,2)=(1-(p_g./lambda_x).^2).*(1-(1.12+s.^(3/2)).^-1).*...
                        (s.^2+1.21*lambda_x)./(s.*(s+1));


            alpha=min(alpha,[],3);


            %% Stresses
            beta=0.6;
            eta=1;

            lambda_x=repmat(sqrt(lambda_x),size(p_hMin,1),1);
            lambda_0=repmat(0.2,size(p_hMin));
            lambda_p=sqrt(alpha./(1-beta));

            chi_x=ones(size(lambda_p));
            idx=lambda_0<lambda_x & lambda_x<lambda_p;
            chi_x(idx)=1-beta.*((lambda_x(idx)-lambda_0(idx))./...
                        (lambda_p(idx)-lambda_0(idx))).^eta;
            idx=lambda_p<=lambda_x;
            chi_x(idx)=alpha(idx)./lambda_x(idx).^2;

            sigma_xRk=chi_x.*f_pTheta;

            gamma_M1=1.1;
            sigma_xRd=sigma_xRk./gamma_M1;
            t_vec=n_x.*Silo.gamma_F./sigma_xRd;

            t_min=reshape(max(t_vec,[],1),sz);


            %% Auxiliary function
            function p=pfx(p_h)
                p=p_h.*r./(t.*sigma_xRcr);
            end
        end


        %Silo geometry to achieve minimal surface area
        function [A,d_cyl,h_cyl,h_hop,h_cylEff,h_heap,h_r]=minSurf(V,d_h,beta,phi_r,gamma,eta_free)
            if isnan(V)
                V=1e6;
            end

            sz=implExp.size(V,d_h,beta,phi_r,gamma,eta_free);
            [V,d_h,beta,phi_r,gamma,eta_free]=implExp.normalize(sz,V,d_h,beta,phi_r,gamma,eta_free);

            [d_cyl,h_cyl,h_cylEff,h_hop,h_heap,h_r]=Silo.getSizes(V,d_h,phi_r,beta,gamma,1e5,[0.2,2]);
            
        
            A_h=pi.*(d_cyl+d_h)./2.*h_hop./cos(beta);
            A_cEff=d_cyl.*pi.*h_cylEff;
            A_cFree=d_cyl.*pi.*(h_cyl-h_cylEff);
            A_r=d_cyl.^2.*pi./(4.*cos(gamma));
            
            A_eff=A_h+A_cEff+(A_cFree+A_r).*eta_free;


            [~,idx]=min(A_eff,[],1);
            getVal=@(mat) arrayfun(@(x) mat(idx(x),x),1:length(d_h));

            d_cyl=reshape(getVal(d_cyl),sz);
            h_cyl=reshape(getVal(h_cyl),sz);
            h_hop=reshape(getVal(h_hop),sz);
            h_cylEff=reshape(getVal(h_cylEff),sz);
            h_heap=reshape(getVal(h_heap),sz);
            h_r=reshape(getVal(h_r),sz);

            A_h=reshape(getVal(A_h),sz);
            A_cEff=reshape(getVal(A_cEff),sz);
            A_cFree=reshape(getVal(A_cFree),sz);
            A_r=reshape(getVal(A_r),sz);
            A_eff=reshape(getVal(A_eff),sz);

            A=struct('A_h',A_h,'A_cEff',A_cEff,'A_cFree',A_cFree,'A_r',A_r,'A_eff',A_eff);

        end


        %Battery geometry to achieve minimal surface area
        function [bat,d_cyl,h_cyl,h_hop,h_cylEff,h_heap,h_r]=minSurfBat(V,d_h,beta,phi_r,gamma,nHop)
            sz=implExp.size(V,d_h,beta,phi_r,gamma,nHop);
            [V,d_h,beta,phi_r,gamma,nHop]=implExp.normalize(sz,V,d_h,beta,phi_r,gamma,nHop);

            
            %% Ideal hopper arrangement
            l_bat=NaN(size(V));
            w_bat=NaN(size(V));
            U=NaN(size(V));
            A_roof=NaN(size(V));
            for i=1:numel(U)
                f=factor(nHop(i));
                switch length(f)
                    case 1
                        n=[1,f];
                    case 2
                        n=f;
                    otherwise
                        nPicks=floor(length(f)/2);
                        n=NaN(nPicks,2);
                
                        for j=1:nPicks
                            c=unique(nchoosek(f,j),'rows');
                    
                            nLoc=NaN(size(c,1),2);
                            nLoc(:,1)=prod(c,2);
                            nLoc(:,2)=nHop(i)./nLoc(:,1);
                
                            [~,idx]=min(sum(nLoc,2));
                            n(j,:)=nLoc(idx,:);
                        end
                
                        [~,idx]=min(sum(n,2));
                        n=n(idx,:);
                end
                l_bat(i)=n(1);
                w_bat(i)=n(2);
                U(i)=2.*sum(n,2);
                A_roof(i)=prod(n,2);
            end
            
            
            %%
            Vi=V./nHop;
            
            [d_cyl,h_cyl,h_cylEff,h_hop,h_heap,h_r]=Silo.getSizes(Vi,d_h,phi_r,beta,gamma,1e4,[0.2,20]);
            h_bat=h_cyl+h_r;
            A_bat=U.*d_cyl.*h_bat+A_roof.*d_cyl.^2;
            
            A_hop=nHop.*pi.*(d_cyl+d_h)./2.*h_hop./cos(beta);
            
            A_tot=A_bat+A_hop;
            
            
            [A_tot,idx]=min(A_tot,[],1);
            getVal=@(mat) arrayfun(@(x) mat(idx(x),x),1:length(V));

            A_tot=reshape(A_tot,sz);
            d_cyl=reshape(getVal(d_cyl),sz);
            h_cyl=reshape(getVal(h_cyl),sz);
            h_cylEff=reshape(getVal(h_cylEff),sz);
            h_hop=reshape(getVal(h_hop),sz);
            h_heap=reshape(getVal(h_heap),sz);
            h_r=reshape(getVal(h_r),sz);
            h_bat=reshape(getVal(h_bat),sz);

            bat=struct('l',reshape(l_bat,sz),'w',reshape(w_bat,sz),'h',h_bat,'A_tot',A_tot);
        end


        % function [d_cyl,h_cyl,h_hop,h_cylEff,h_heap,h_r]=minHeight(V,d_h,beta,phi_r,gamma)
        %     sz=implExp.size(V,d_h,beta,phi_r,gamma);
        %     [V,d_h,beta,phi_r,gamma]=implExp.normalize(sz,V,d_h,beta,phi_r,gamma);
        %     l=1:length(d_h);
        % 
        %     [d_cyl,h_cyl,h_cylEff,h_hop,h_heap,h_r]=Silo.getSizes(V,d_h,phi_r,beta,gamma,1e5);
        % 
        % 
        %     H=h_hop+h_cyl+h_r;
        % 
        % 
        %     [~,idx]=min(H,[],1);
        %     getVal=@(mat) arrayfun(@(x) mat(idx(x),x),l);
        % 
        %     d_cyl=reshape(getVal(d_cyl),sz);
        %     h_cyl=reshape(getVal(h_cyl),sz);
        %     h_hop=reshape(getVal(h_hop),sz);
        %     h_cylEff=reshape(getVal(h_cylEff),sz);
        %     h_heap=reshape(getVal(h_heap),sz);
        %     h_r=reshape(getVal(h_r),sz);
        % end


        %Silo geometry to achieve minimal steel mass
        function [d_cyl,h_cyl,h_hop,h_cylEff,h_heap,h_r]=minMass(V,d_h,beta,phi_r,gamma,T,steel,wall,particles,nsz,slenderRange)
            if ischar(steel)
                steel={steel};
            end
            if nargin<10
                nsz=1e2;
            end
            if nargin<11
                slenderRange=[0.4,2];
            end
            nz=120;
            
            sz=implExp.size(V,d_h,beta,phi_r,gamma,T,steel,wall);
            [V,d_h,beta,phi_r,gamma,T,steel,wall]=implExp.normalize(sz,V,d_h,beta,phi_r,gamma,T,steel,wall);

            

            [d_c,h_cyl,h_cylEff,h_hop,h_heap,h_r]=Silo.getSizes(V,d_h,phi_r,beta,gamma,nsz,slenderRange);
            
            h_o=d_c./(6*tan(phi_r));
            h_c=h_cylEff+h_o;
            t_abr=2e-3;


            c=0;
            qCount=parallel.pool.DataQueue;
            afterEach(qCount,@cfx);

            qWarn=parallel.pool.DataQueue;
            afterEach(qWarn,@(i) warning('Limit value reached at Index %d \n',i));


            idx=NaN(size(V));
            parfor i=1:length(V)
                tic;
                if any(isnan(d_c(:,i)))
                    send(qCount,toc);
                    continue
                end


                [~,coords,~,LS1LS3]=Silo.getTcyl(h_c(:,i),d_c(:,i),T(i),steel{i},wall(i),particles,24);

                dz=diff(coords.z(1:2,:));
                n=zeros(4,length(dz));
                nAdapt=round(nz/4);


                zBound=NaN(5,length(dz));
                zBound(1,:)=Silo.getZ(coords.z,0,'last',LS1LS3);
                zBound(2,:)=Silo.getZ(coords.z,1,'first',LS1LS3);
                zBound(3,:)=Silo.getZ(coords.z,1,'last',LS1LS3);
                zBound(4,:)=Silo.getZ(coords.z,2,'first',LS1LS3);
                zBound(5,:)=coords.z(end,:);


                val=isnan(zBound(1,:));
                zBound(1,val)=zBound(2,val);
                zBound(1,:)=zBound(1,:)-dz./2;
                zBound(1,zBound(1,:)<0)=0;


                n(1,~isnan(zBound(2,:)))=nAdapt;
                n(2,~(isnan(zBound(3,:)) | zBound(2,:)==zBound(3,:)))=nAdapt;
                n(3,~isnan(zBound(4,:)))=nAdapt;
                n(4,~isnan(zBound(4,:)) & zBound(4,:)~=zBound(5,:))=nAdapt;

                n=round(n.*nz./sum(n,1));
                [~,idx2]=max(n,[],1);
                for j=1:length(idx2)
                    n(idx2(j),j)=n(idx2(j),j)+(nz-sum(n(:,j),1));
                end
                n=cumsum(n,1);
                n=[ones(1,nsz);n];


                val=isnan(zBound(3,:)) & ~isnan(zBound(4,:));
                zBound(3,val)=zBound(4,val)-dz(val)./2;
                zBound(3,zBound(3,:)<0)=0;

                zBound(isnan(zBound))=0;

                
                z=NaN(nz,nsz);
                for j=1:nsz
                    for k=1:size(n,1)-1
                        s=n(k,j);
                        e=n(k+1,j);
                        z(s:e,j)=linspace(zBound(k,j),zBound(k+1,j),e-s+1);
                    end
                end


                [tCylVec,coords,p_vf]=Silo.getTcyl(h_c(:,i),d_c(:,i),T(i),steel{i},wall(i),particles,nz,z);



                tCylVec=tCylVec+t_abr;

                
                %Hopper
                p_vf=p_vf(end,:)';
                [t,~,x]=Silo.hopper(d_c(:,i),d_h(i),beta(i),p_vf,T(i),steel{i},particles,wall(i));
                tHopVec=t.vec;
                
                
                %Geometry
                
                tHopVec=tHopVec+t_abr;
                tRoof=2e-3;
                
                d_oc=d_c(:,i)'+2.*tCylVec;
                d_cLoc=d_c(:,i);
                d_ih=cell2mat(arrayfun(@(j) linspace(d_h(i),d_cLoc(j),size(x,1))',1:nsz,'UniformOutput',false));
                d_oh=d_ih+2.*tHopVec./cos(beta(i));
                

                if gamma(i)>0
                    d_or=d_c(:,i)'+2.*tRoof./sin(gamma(i));
                    Vroof=h_r(:,i).*pi./12.*(d_or'.^2-d_c(:,i).^2);
                else
                    d_or=d_oc;
                    Vroof=d_or.^2.*pi./4.*tRoof;
                end


                f=(d_oc.^2-d_c(:,i)'.^2).*pi./4;
                Vcyl=arrayfun(@(j) trapz(coords.z(:,j),f(:,j),1),1:nsz)';

                f=(d_oh.^2-d_ih.^2).*pi./4;
                Vhop=arrayfun(@(j) trapz(x(:,j),f(:,j),1),1:nsz)';
                
                
                Vtot=Vcyl+Vhop+Vroof;
                [~,idx(i)]=min(Vtot);
                if idx(i)==1 || idx(i)==nsz
                    send(qWarn,i);
                end

                send(qCount,toc);
            end
            delete([qCount,qWarn]);
            

            val=~isnan(idx);
            getVal=@(mat) arrayfun(@(x) mat(idx(x),x),find(val));

            d_cyl=NaN(size(V));
            d_cyl(val)=getVal(d_c);
            h_cyl(1,val)=getVal(h_cyl);
            h_hop(1,val)=getVal(h_hop);
            h_cylEff(1,val)=getVal(h_cylEff);
            h_heap(1,val)=getVal(h_heap);
            h_r(1,val)=getVal(h_r);

            d_cyl=reshape(d_cyl,sz);
            h_cyl=reshape(h_cyl(1,:),sz);
            h_hop=reshape(h_hop(1,:),sz);
            h_cylEff=reshape(h_cylEff(1,:),sz);
            h_heap=reshape(h_heap(1,:),sz);
            h_r=reshape(h_r(1,:),sz);


            function cfx(t)
                c=c+1;
                disp(['i=',num2str(c),'/',num2str(length(V)),', Time=',num2str(t)]);
            end

        end


        %Read constants from Excel table and write it into a matlab file
        function loadConstants()
            tab=readtable(['@Silo',filesep,'constants.xls'],'Sheet','thermalFactors');
            tab.T=tab.T+273.15;
            k_p=@(T) interp1(tab.T,tab.k_p,T);
            k_E=@(T) interp1(tab.T,tab.k_E,T);

            material=readtable(['@Silo',filesep,'constants.xls'],'Sheet','materials');
            material.Properties.RowNames=material.Material;
            material.Material=[];
            material{:,:}=material{:,:}.*10^6;

            particles=readtable(['@Silo',filesep,'constants.xls'],'Sheet','particles');
            particles.gamma_l=particles.gamma_l.*10^3;
            particles.gamma_u=particles.gamma_u.*10^3;
            particles.phi_r=deg2rad(particles.phi_r);
            particles.phi_im=deg2rad(particles.phi_im);

            save(['@Silo',filesep,'constants.mat'],'k_p','k_E','material','particles');
        end
    end


    methods(Static,Access=private)
        function [d_cyl,h_cyl,h_cylEff,h_hop,h_heap,h_r]=getSizes(V,d_h,phi_r,beta,gamma,n,slenderInterval)
            sz=implExp.size(V,d_h,phi_r,beta,gamma);
            [V,d_h,phi_r,beta,gamma]=implExp.normalize(sz,V,d_h,phi_r,beta,gamma);
            l=1:length(V);

            d_cApprox=(V.*4./pi).^(1/3);
            
            h_heapfx=@(d_c,i) d_c.*tan(phi_r(i))./2;
            V_heapfx=@(d_c,i) d_c.^2.*pi.*h_heapfx(d_c,i)./12;
        
            h_hfx=@(d_c,d_h,i) (d_c-d_h)./(2*tan(beta(i)));
            V_hfx=@(d_c,d_h,i) h_hfx(d_c,d_h,i).*pi.*(d_c.^2+d_c.*d_h+d_h.^2)./12;

            h_cylEffFx=@(d_c,d_h,i) (V(i)-V_heapfx(d_c,i)-V_hfx(d_c,d_h,i))./(d_c.^2.*pi./4);
            h_rfx=@(d_c,i) d_c./2.*tan(gamma(i));
            h_cylfx=@(d_c,d_h,i) h_cylEffFx(d_c,d_h,i)-h_rfx(d_c,i)+h_heapfx(d_c,i);


            d_cMax=arrayfun(@(i) fzero(@(d_c) V_hfx(d_c,0,i)+V_heapfx(d_c,i)-V(i),d_cApprox(i)),l);

            slenderMin=h_cylfx(d_cMax,d_h,l)./d_cMax;
            idx=slenderInterval(1)>slenderMin;
            d_cMax(idx)=arrayfun(@(i) fzero(@(d_c) h_cylfx(d_c,d_h(i),i)./d_c-slenderInterval(1),d_cMax(i)),find(idx));


            d_cMin=arrayfun(@(i) fzero(@(d_c) h_cylfx(d_c,d_h(i),i)./d_c-slenderInterval(2),d_cMax(i)),l);

            idx=d_cMin<d_h;
            d_cMin(idx)=d_h(idx);


            d_cyl=NaN(n,length(V));
            idx=d_h<d_cMax;
            d_cyl(:,idx)=cell2mat(arrayfun(@(x) linspace(d_cMin(x)+1e-6,d_cMax(x),n)',find(idx),'UniformOutput',false));
            % d_cyl(d_cyl<d_h)=NaN;


            h_heap=h_heapfx(d_cyl,l);
            h_hop=h_hfx(d_cyl,d_h,l);
            h_cylEff=h_cylEffFx(d_cyl,d_h,l);
            h_r=h_rfx(d_cyl,l);
            h_cyl=h_cylfx(d_cyl,d_h,l);



            % V_heap=V_heapfx(d_cyl,l);
            % V_h=V_hfx(d_cyl,d_h,l);
            % idx=h_cylEff<0;
            % h_cylEff(idx)=NaN;
            % d_cyl(idx)=NaN;
            % h_heap(idx)=NaN;
            % h_hop(idx)=NaN;

            
        end


        function [tCylVec,coords,p_vf,LS1LS3]=getTcyl(h_c,d_c,T,steel,wall,particles,nz,z)
            if nargin<8
                z=NaN;
            end

            nsz=length(h_c);

            %Cylinder
                t_minLS1=NaN(nsz,2);
                t_minLS1vec=NaN(nz,nsz,2);
                p_hMin=NaN(nz,nsz,2);
                
                silfx=@(critProp) Silo.cylinder(h_c,d_c,T,steel,wall,particles,critProp,nz,z);
                
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
                t_minLS3=NaN(size(t_minLS1));
                t_minLS3vec=NaN(size(t_minLS1vec));
                for i=1:nsz
                    b1=1;
                    count=0;
                    isInBound=false;
                    while ~isInBound && count<4
                        try
                            t_minLS3(i)=fzero(@(t) Silo.buckling(d_c(i),t,p_hMin(:,i),p_hMax(:,i),n_xMax(:,i),T,steel)-t,...
                                                        [t_minLS1(i)/3,t_minLS1(i)*3*b1]);
                            isInBound=true;
                        catch ME
                            if strcmp(ME.identifier,'MATLAB:fzero:ValuesAtEndPtsSameSign')
                                b1=b1*10;
                            else
                                rethrow(ME);
                            end
                        end
                        count=count+1;
                    end
                    if ~isInBound
                        t_minLS3(i)=t_minLS1(i)*3*b1;
                        warning('Could not find t_minLS3 at index %d',i)
                    end
                    [~,t_minLS3vec(:,i)]=Silo.buckling(d_c(i),t_minLS3(i),p_hMin(:,i),p_hMax(:,i),n_xMax(:,i),T,steel);
                end


                %Critical wall thickness
                [tCylVec,LS1LS3]=max(cat(3,t_minLS1vec,t_minLS3vec),[],3);
                LS1LS3(tCylVec==0)=0;
        end


        function z=getZ(zmat,i,order,LS1LS3)
            z=NaN(1,size(LS1LS3,2));

            val=any(LS1LS3==i,1);
            idx=NaN(1,length(z));
            idx(val)=arrayfun(@(j) find(LS1LS3(:,j)==i,1,order),find(val));
            z(val)=arrayfun(@(j) zmat(idx(j),j),find(val));
        end


        function part=getPartProp(name,wall,critProp)
            persistent tab
            if isempty(tab)
                partStruct=load(['@Silo',filesep,'constants.mat'],'particles');
                tab=partStruct.particles;
            end

            tabloc=tab(strcmp(tab.Name,name),:);


            gamma=tabloc.gamma_u;
            phi_r=tabloc.phi_r;
            C_op=tabloc.C_op;
            
            phi_i=tabloc.phi_im*[1/tabloc.a_phi,tabloc.a_phi];   

            K=tabloc.K_m*[1/tabloc.a_K,tabloc.a_K];
            
            my_m=[tabloc.my_D1,tabloc.my_D2,tabloc.my_D3]';
            my=my_m*[1/tabloc.a_my,tabloc.a_my];  %my<=tan(phi_i) always
            my=my(wall,:);
            

            switch critProp
                case 1  %Maximum normal pressure on vertical wall
                    phi_i=phi_i(1);
                    my=min(my(1),tan(phi_i));
                    K=K(2);
                case 2  %Maximum frictional traction on vertical wall
                    phi_i=phi_i(1);
                    my=min(my(2),tan(phi_i));
                    K=K(2);
                case 3  %Maximum vertical load on hopper or silo bottom
                    phi_i=phi_i(2);
                    my=min(my(1),tan(phi_i));
                    K=K(1);
            end


            part=struct('gamma',gamma,'phi_r',phi_r,'C_op',C_op,'phi_i',phi_i,...
                        'K',K,'my',my);
        end


        function t_min=vonMises(n_x,n_theta,n_xTheta,T,steel)
            eqMax=max(sqrt(n_x.^2+n_theta.^2-...
                                n_x.*n_theta+...
                                3*n_xTheta.^2),[],3);
            
            f_pTheta=Silo.yieldElast(T,steel);
            
            gamma_M0=1;
            f_eRD=f_pTheta./gamma_M0;
        
            t_min=eqMax.*Silo.gamma_F./f_eRD;
        end


        function [f_pTheta,E_aTheta]=yieldElast(T,steel)
            E=210e9;

            matStruct=load(['@Silo',filesep,'constants.mat'],'material','k_p','k_E');
            mat=matStruct.material;
            k_p=matStruct.k_p;
            k_E=matStruct.k_E;
        
            f_y=mat.f_y40(steel);
            f_pTheta=f_y.*k_p(T);

            E_aTheta=E.*k_E(T);
        end
    end
end




