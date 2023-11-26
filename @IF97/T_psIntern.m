function T=T_psIntern(p,s,regions,region2,region3)
    T=NaN(size(p));
    T(regions(1))=T_ps1(p(regions(1)),s(regions(1)));
    
    T(region2(1))=T_ps2a(p(region2(1)),s(region2(1)));
    T(region2(2))=T_ps2b(p(region2(2)),s(region2(2)));
    T(region2(3))=T_ps2c(p(region2(3)),s(region2(3)));
    
    T(region3(1))=T_ps3a(p(region3(1)),s(region3(1)));
    T(region3(2))=T_ps3b(p(region3(2)),s(region3(2)));
    
    T(regions(4))=IF97.T_satIntern(p(regions(4)));
end


function T=T_ps1(p,s)
    persistent I J n p_star s_star
    if isempty(p_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'Region1constants.mat'],'I_ps','J_ps','n_ps');
        I=vars.I_ps;
        J=vars.J_ps;
        n=vars.n_ps;
        p_star=10^6;
        s_star=10^3;
    end
    
    
    if ~isempty(p)
        T=sum(n.*(p./p_star).^I.*(s./s_star+2).^J,1);
    else
        T=NaN(1,0);
    end
end


function T=T_ps2a(p,s)
    persistent I J n p_star s_star
    if isempty(p_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'Region2constants.mat'],'I_psa','J_psa','n_psa');
        I=vars.I_psa;
        J=vars.J_psa;
        n=vars.n_psa;
        p_star=10^6;
        s_star=2*10^3;
    end
    
    
    if ~isempty(p)
        T=sum(n.*(p./p_star).^I.*(s./s_star-2).^J,1);
    else
        T=NaN(1,0);
    end
end


function T=T_ps2b(p,s)
    persistent I J n p_star s_star
    if isempty(p_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'Region2constants.mat'],'I_psb','J_psb','n_psb');
        I=vars.I_psb;
        J=vars.J_psb;
        n=vars.n_psb;
        p_star=10^6;
        s_star=0.7853*10^3;
    end
    
    
    if ~isempty(p)
        T=sum(n.*(p./p_star).^I.*(10-s./s_star).^J,1);
    else
        T=NaN(1,0);
    end
end


function T=T_ps2c(p,s)
    persistent I J n p_star s_star
    if isempty(p_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'Region2constants.mat'],'I_psc','J_psc','n_psc');
        I=vars.I_psc;
        J=vars.J_psc;
        n=vars.n_psc;
        p_star=10^6;
        s_star=2.9251*10^3;
    end
    
    
    if ~isempty(p)
        T=sum(n.*(p./p_star).^I.*(2-s./s_star).^J,1);
    else
        T=NaN(1,0);
    end
end


function T=T_ps3a(p,s)
    persistent I J n T_star p_star s_star
    if isempty(p_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'Region3constants.mat'],'I_psa','J_psa','n_psa');
        I=vars.I_psa;
        J=vars.J_psa;
        n=vars.n_psa;
        T_star=760;
        p_star=100*10^6;
        s_star=4.4*10^3;
    end
    
    
    if ~isempty(p)
        T=T_star.*sum(n.*(p./p_star+0.24).^I.*(s./s_star-0.703).^J,1);
    else
        T=NaN(1,0);
    end
end


function T=T_ps3b(p,s)
    persistent I J n T_star p_star s_star
    if isempty(p_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'Region3constants.mat'],'I_psb','J_psb','n_psb');
        I=vars.I_psb;
        J=vars.J_psb;
        n=vars.n_psb;
        T_star=860;
        p_star=100*10^6;
        s_star=5.3*10^3;
    end
    
    
    if ~isempty(p)
        T=T_star.*sum(n.*(p./p_star+0.76).^I.*(s./s_star-0.818).^J,1);
    else
        T=NaN(1,0);
    end
end


