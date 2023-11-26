function T=T_ph(p,h)
    sz=IF97.implExpSize(p,h);
    [regions,region2,region3]=phregions(p,h,sz);
    [p,h]=IF97.normalize(sz,p,h);
    
    T=NaN(sz);
    T(regions(1))=T_ph1(p(regions(1)),h(regions(1)));
    
    T(region2(1))=T_ph2a(p(region2(1)),h(region2(1)));
    T(region2(2))=T_ph2b(p(region2(2)),h(region2(2)));
    T(region2(3))=T_ph2c(p(region2(3)),h(region2(3)));
    
    T(region3(1))=T_ph3a(p(region3(1)),h(region3(1)));
    T(region3(2))=T_ph3b(p(region3(2)),h(region3(2)));
    
    T(regions(4))=IF97.T_satIntern(p(regions(4)));
end


function T=T_ph1(p,h)
    persistent I J n p_star h_star
    if isempty(p_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'Region1constants.mat'],'I_ph','J_ph','n_ph');
        I=vars.I_ph;
        J=vars.J_ph;
        n=vars.n_ph;
        p_star=10^6;
        h_star=2500*10^3;
    end
    
    
    if ~isempty(p)
        T=sum(n.*(p./p_star).^I.*(h./h_star+1).^J,1);
    else
        T=NaN(1,0);
    end
end


function T=T_ph2a(p,h)
    persistent I J n p_star h_star
    if isempty(p_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'Region2constants.mat'],'I_pha','J_pha','n_pha');
        I=vars.I_pha;
        J=vars.J_pha;
        n=vars.n_pha;
        p_star=10^6;
        h_star=2000*10^3;
    end
    
    
    if ~isempty(p)
        T=sum(n.*(p./p_star).^I.*(h./h_star-2.1).^J,1);
    else
        T=NaN(1,0);
    end
end


function T=T_ph2b(p,h)
    persistent I J n p_star h_star
    if isempty(p_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'Region2constants.mat'],'I_phb','J_phb','n_phb');
        I=vars.I_phb;
        J=vars.J_phb;
        n=vars.n_phb;
        p_star=10^6;
        h_star=2000*10^3;
    end
    
    
    if ~isempty(p)
        T=sum(n.*(p./p_star-2).^I.*(h./h_star-2.6).^J,1);
    else
        T=NaN(1,0);
    end
end


function T=T_ph2c(p,h)
    persistent I J n p_star h_star
    if isempty(p_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'Region2constants.mat'],'I_phc','J_phc','n_phc');
        I=vars.I_phc;
        J=vars.J_phc;
        n=vars.n_phc;
        p_star=10^6;
        h_star=2000*10^3;
    end
    
    
    if ~isempty(p)
        T=sum(n.*(p./p_star+25).^I.*(h./h_star-1.8).^J,1);
    else
        T=NaN(1,0);
    end
end


function T=T_ph3a(p,h)
    persistent I J n T_star p_star h_star
    if isempty(p_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'Region3constants.mat'],'I_pha','J_pha','n_pha');
        I=vars.I_pha;
        J=vars.J_pha;
        n=vars.n_pha;
        T_star=760;
        p_star=100*10^6;
        h_star=2300*10^3;
    end
    
    
    if ~isempty(p)
        T=T_star.*sum(n.*(p./p_star+0.24).^I.*(h./h_star-0.615).^J,1);
    else
        T=NaN(1,0);
    end
end


function T=T_ph3b(p,h)
    persistent I J n T_star p_star h_star
    if isempty(p_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'Region3constants.mat'],'I_phb','J_phb','n_phb');
        I=vars.I_phb;
        J=vars.J_phb;
        n=vars.n_phb;
        T_star=860;
        p_star=100*10^6;
        h_star=2800*10^3;
    end
    
    
    if ~isempty(p)
        T=T_star.*sum(n.*(p./p_star+0.298).^I.*(h./h_star-0.72).^J,1);
    else
        T=NaN(1,0);
    end
end


function [regions,region2,region3]=phregions(p,h,sz)
    persistent p_sat273 p_sat623
    if isempty(p_sat273)
        p_sat273=IF97.p_satIntern(273.15);
        p_sat623=IF97.p_satIntern(623.15);
    end
    
    
%     T_sat=IF97.T_sat(p);
    h0=IF97.h(p,NaN,0);
    h1=IF97.h(p,NaN,1);
    
    
    szp=size(p);
    p=reshape(p,1,numel(p));
    
    
    hmin=NaN(szp);
    reg=p_sat273<=p & p<=100*10^6;
    hmin(reg)=hreg(273.15,1);
    reg=0<p & p<p_sat273;
    hmin(reg)=hreg(273.15,2);
    
    
    h13=NaN(szp);
    reg=p_sat623<=p & p<=100*10^6;
    h13(reg)=hreg(623.15,1);
    
    
    h23=NaN(szp);
    h23(reg)=hreg(IF97.B23(p(reg),false),2);
    
    
    hmax2=NaN(szp);
    reg=0<p & p<=100*10^6;
    hmax2(reg)=hreg(1073.15,2);
    
    
    hmax5=NaN(szp);
    reg=0<p & p<=50*10^6;
    hmax5(reg)=hreg(2273.15,5);


    p=reshape(p,szp);
    page=IF97.pagefx(sz);
    regionmx=false([sz,5]);

    regionmx(page(1))=(p_sat273<=p & p<=p_sat623 & hmin<=h & h<h0) | ... 
                        (p_sat623<p & p<=100*10^6 & hmin<=h & h<=h13);
    regionmx(page(2))=(0<p & p<p_sat273 & hmin<=h & h<=hmax2) | ...
                    (p_sat273<=p & p<=p_sat623 & h1<h & h<=hmax2) | ...
                    (p_sat623<p & p<=100*10^6 & h23<=h & h<=hmax2);
    regionmx(page(3))=(p_sat623<p & p<=IF97.p_c & h13<h & h<h0) | ...
                    (p_sat623<p & p<=IF97.p_c & h1<h & h<h23) | ...
                    (IF97.p_c<p & p<=100*10^6 & h13<h & h<h23);
    regionmx(page(4))=p_sat273<=p & p<=IF97.p_c & h0<=h & h<=h1;
    regionmx(page(5))=0<p & p<=50*10^6 & hmax2<h & h<=hmax5;
    
    regions=@(x) regionmx(page(x));
    
    
    region2mx=false([sz,3]);
    region2mx(page(1))=reshape(regions(2),sz) & p<=4*10^6 & true(size(h));
    region2mx(page(3))=reshape(regions(2),sz) & 6.5467*10^6<p & h<B2bc(p);
    region2mx(page(2))=regions(2) & ~region2mx(page(1)) & ~region2mx(page(3));
    
    region2=@(x) region2mx(page(x));
    
    
    region3mx=false([sz,2]);
    region3mx(page(1))=(reshape(regions(3),sz) & p<=IF97.p_c & h<h0) | ...
                        (reshape(regions(3),sz) & IF97.p_c<p & h<=B3ab(p));
    region3mx(page(2))=regions(3) & ~region3mx(page(1));
    
    region3=@(x) region3mx(page(x));
    
    
    function h=hreg(T,i)
        if ~isempty(T)
            n=nnz(reg);
            exp=@(x) repmat(x,1,n./numel(x));
            h=IF97.h_pTx(p(reg),exp(T),NaN(1,n),@(x) exp(x==i),true(1,n));
        else
            h=NaN(1,0);
        end
    end
end


function h2bc=B2bc(p)
    persistent n p_star h_star
    if isempty(p_star)
        folder=IF97.getFolder();
        nstruct=load([folder,filesep,'BoundaryConstants.mat'],'n_B2bc');
        n=nstruct.n_B2bc;
        p_star=10^6;
        h_star=10^3;
    end
    
    h2bc=h_star*(n(4)+sqrt((p/p_star-n(5))/n(3)));
end


function h3ab=B3ab(p)
    persistent n p_star h_star
    if isempty(p_star)
        folder=IF97.getFolder();
        nstruct=load([folder,filesep,'BoundaryConstants.mat'],'n_B3ab');
        n=nstruct.n_B3ab;
        p_star=10^6;
        h_star=10^3;
    end
    
    h3ab=h_star*sum(n.*(p/p_star).^((0:3)'),1);
end




