function [regions,region2,region3]=hsregions(h,s,sz)
    persistent s_satL273 s_satV273 s_satL623 p_sat273 s_B2abmin s_B2abmax s_B13min s_B23min s_B23max h_B23min h_B23max s_pmax h_pmax smax hmax2
    if isempty(s_satL273)
        s_satL273=-1.545495919*10^-1;
        s_satV273=9.155759395*10^3;
        s_satL623=3.77828134*10^3;
        
        p_sat273=IF97.p_satIntern(273.15);
        
        s_B2abmin=IF97.s(4*10^6,NaN,1);
        s_B2abmax=IF97.s(4*10^6,1073.15);
        
        s_B13min=3.397782955*10^3;
        
        s_B23min=5.048096828*10^3;
        s_B23max=5.260578707*10^3;
        h_B23min=2.563592004*10^6;
        h_B23max=2.812942061*10^6;
        
        s_pmax=IF97.s(100*10^6,1073.15);
        h_pmax=IF97.h(100*10^6,1073.15);
        
        smax=IF97.s(p_sat273,1073.15);
        hmax2=IF97.h(p_sat273,IF97.T_ps(p_sat273,smax));
    end
    
    
    hnorm=h;
    snorm=s;
    
    reg=h_B23min<=h & h<=h_B23max;
    hnorm(~reg)=NaN;
    
    reg=s_B23min<=s & s<=s_B23max;
    snorm(~reg)=NaN;
    
    [hnorm,snorm]=IF97.normalize(sz,hnorm,snorm);
    
    reg=~isnan(hnorm) & ~isnan(snorm);
    hnorm=hnorm(reg);
    snorm=snorm(reg);
    
    p23=NaN(sz);
    p23(reg)=IF97.B23(B23(hnorm,snorm),true);
    
    p2c=NaN(sz);
    p2c(reg)=p_hsreg2(hnorm,snorm,3);
    
    
    szs=size(s);
    s=reshape(s,1,numel(s));
    
    
    hmax=NaN(szs);
    reg=s_satL273<=s & s<=s_pmax;
    if any(reg)
        hmax(reg)=IF97.h(100*10^6,IF97.T_ps(100*10^6,s(reg)));
    end
    
    reg=true;
    for i=find(s_pmax<s & s<=smax)
        hmax(i)=fzero(@(h) T_psreg2(p_hsreg2(h,s(i),1),s(i),1)-1073.15,[h_pmax,hmax2]);
    end
    
    
    hmin=NaN(szs);
    reg=s_satL273<=s & s<=s_satV273;
    [exp,n]=expfx();
    hminfx(exp(273.15),(s(reg)-s_satL273)./(s_satV273-s_satL273),4);
    
    reg=s_satV273<s & s<=10*10^3;
    [exp,n]=expfx();
    hminfx(T_psreg2(p_sat273,s(reg),1),NaN(1,n),2);
    
    
    h1=NaN(szs);
    reg=s_satL273<=s & s<=s_satL623;
    h1(reg)=B1(s(reg));
    
    
    h3a=NaN(szs);
    reg=s_satL623<s & s<=IF97.s_c;
    h3a(reg)=B3a(s(reg));
    
    
    h2ab_satV=NaN(szs);
    reg=5.85*10^3<s & s<=s_satV273;
    h2ab_satV(reg)=B2ab_satV(s(reg));
    
    
    h2c3b=NaN(szs);
    reg=IF97.s_c<s & s<=5.85*10^3;
    h2c3b(reg)=B2c3b(s(reg));
    
    
    h13=NaN(szs);
    reg=s_B13min<=s & s<=s_satL623;
    h13(reg)=B13(s(reg));
    
    
    h2ab=NaN(szs);
    reg=s_B2abmin<s & s<=s_B2abmax;
    h2ab(reg)=B2ab(s(reg));
    
    
    
    s=reshape(s,szs);
    page=IF97.pagefx(sz);
    regionmx=false([sz,5]);

    regionmx(page(1))=(s<s_B13min & h1<h & h<=hmax) | ... 
                    (h1<h & h<=h13);
    regionmx(page(2))=(s_B23max<s & h2c3b<h & h<=hmax) | ...
                    (h2ab_satV<h & h<=hmax) | ...
                    (s_satV273<s & s<=10*10^3 & hmin<=h & h<=hmax) | ...
                    (s_B23min<=s & s<=s_B23max & h_B23max<h & h<=hmax) | ...
                    (p2c<=p23);
    regionmx(page(3))=(h13<h & h<=hmax) | ...
                    (h3a<h & h<=hmax) | ...
                    (s<s_B23min & h2c3b<h & h<=hmax) | ...
                    (s_B23min<=s & s<=s_B23max & h2c3b<h & h<h_B23min) | ...
                    (p23<p2c);
    regionmx(page(4))=(hmin<=h & h<=h1) | ...
                    (hmin<=h & h<=h3a) | ...
                    (hmin<=h & h<=h2c3b) | ...
                    (hmin<=h & h<=h2ab_satV);
    
    regions=@(x) regionmx(page(x));
    
    
    region2mx=false([sz,3]);
    region2mx(page(1))=(reshape(regions(2),sz) & h<=h2ab) | ...
                        (reshape(regions(2),sz) & s_B2abmax<s & true(size(h)));
    region2mx(page(3))=reshape(regions(2),sz) & s<5.85*10^3 & true(size(h));
    region2mx(page(2))=regions(2) & ~region2mx(page(1)) & ~region2mx(page(3));
    
    region2=@(x) region2mx(page(x));
    
    
    region3mx=false([sz,2]);
    region3mx(page(1))=reshape(regions(3),sz) & s<=IF97.s_c & true(size(h));
    region3mx(page(2))=regions(3) & ~region3mx(page(1));
    
    region3=@(x) region3mx(page(x));
    
    
    function p=p_hsreg2(h,s,regi)
        [exp,n]=expfx();
        p=IF97.p_hsIntern(h,s,@(x) exp(x==2),@(x) exp(x==regi),@(x) false(1,n));
    end
    
    
    function [exp,n]=expfx()
        n=nnz(reg);
        exp=@(x) repmat(x,1,n);
    end


    function hminfx(T,x,regi)
        hmin(reg)=IF97.h_pTx(exp(p_sat273),T,x,@(x) exp(x==regi),exp(regi~=4));
    end
end


function T=T_psreg2(p,s,reg2)
    n=numel(s);
    exp=@(x) repmat(x,1,n);
    T=IF97.T_psIntern(exp(p),s,@(x) exp(x==2),@(x) exp(x==reg2),@(x) false(1,n));
end


function h1=B1(s)
    persistent I J n h_star s_star
    if isempty(h_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'BoundaryConstants.mat'],'I_h1','J_h1','n_h1');
        I=vars.I_h1;
        J=vars.J_h1;
        n=vars.n_h1;
        h_star=1700*10^3;
        s_star=3.8*10^3;
    end
    
    if ~isempty(s)
        s=s./s_star;
        h1=h_star.*sum(n.*(s-1.09).^I.*(s+0.366*10^-4).^J,1);
    else
        h1=NaN(1,0);
    end
end


function h3a=B3a(s)
    persistent I J n h_star s_star
    if isempty(h_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'BoundaryConstants.mat'],'I_h3a','J_h3a','n_h3a');
        I=vars.I_h3a;
        J=vars.J_h3a;
        n=vars.n_h3a;
        h_star=1700*10^3;
        s_star=3.8*10^3;
    end
    
    if ~isempty(s)
        s=s./s_star;
        h3a=h_star.*sum(n.*(s-1.09).^I.*(s+0.366*10^-4).^J,1);
    else
        h3a=NaN(1,0);
    end
end


function h2ab_satV=B2ab_satV(s)
    persistent I J n h_star s1_star s2_star
    if isempty(h_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'BoundaryConstants.mat'],'I_h2ab_satV','J_h2ab_satV','n_h2ab_satV');
        I=vars.I_h2ab_satV;
        J=vars.J_h2ab_satV;
        n=vars.n_h2ab_satV;
        h_star=2800*10^3;
        s1_star=5.21*10^3;
        s2_star=9.2*10^3;
    end
    
    if ~isempty(s)
        h2ab_satV=h_star.*exp(sum(n.*(s1_star./s-0.513).^I.*(s./s2_star-0.524).^J,1));
    else
        h2ab_satV=NaN(1,0);
    end
end


function h2c3b=B2c3b(s)
    persistent I J n h_star s_star
    if isempty(h_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'BoundaryConstants.mat'],'I_h2c3b','J_h2c3b','n_h2c3b');
        I=vars.I_h2c3b;
        J=vars.J_h2c3b;
        n=vars.n_h2c3b;
        h_star=2800*10^3;
        s_star=5.9*10^3;
    end
    
    if ~isempty(s)
        s=s./s_star;
        h2c3b=h_star.*sum(n.*(s-1.02).^I.*(s-0.726).^J,1).^4;
    else
        h2c3b=NaN(1,0);
    end
end


function h13=B13(s)
    persistent I J n h_star s_star
    if isempty(h_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'BoundaryConstants.mat'],'I_hB13','J_hB13','n_hB13');
        I=vars.I_hB13;
        J=vars.J_hB13;
        n=vars.n_hB13;
        h_star=1700*10^3;
        s_star=3.8*10^3;
    end
    
    if ~isempty(s)
    s=s./s_star;
        h13=h_star.*sum(n.*(s-0.884).^I.*(s-0.864).^J,1);
    else
        h13=NaN(1,0);
    end
end


function h2ab=B2ab(s)
    persistent n h_star s_star
    if isempty(h_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'BoundaryConstants.mat'],'n_2ab');
        n=vars.n_2ab;
        h_star=10^3;
        s_star=10^3;
    end
    
    if ~isempty(s)
        h2ab=h_star*sum(n.*(s./s_star).^((0:3)'),1);
    else
        h2ab=NaN(1,0);
    end
end


function T=B23(h,s)
    persistent I J n T_star h_star s_star
    if isempty(T_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'BoundaryConstants.mat'],'I_TB23','J_TB23','n_TB23');
        I=vars.I_TB23;
        J=vars.J_TB23;
        n=vars.n_TB23;
        T_star=900;
        h_star=3000*10^3;
        s_star=5.3*10^3;
    end
    
    if ~isempty(h)
        T=T_star.*sum(n.*(h./h_star-0.727).^I.*(s./s_star-0.864).^J,1);
    else
        T=NaN(1,0);
    end
end




