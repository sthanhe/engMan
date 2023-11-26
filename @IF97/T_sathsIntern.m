function T=T_sathsIntern(h,s,regions)
    persistent s_satV623
    if isempty(s_satV623)
        s_satV623=5.210887825*10^3;
    end
    
    
    reg41=regions(4) & s_satV623<=s;
    
    
    T=NaN(size(h));
    T(reg41)=T_sat41(h(reg41),s(reg41));
    
    for i=find(~reg41)
        T(i)=fzero(@(T) deltax(h(i),s(i),T),[273.15,IF97.T_c-10^-10]);
    end
end


function T=T_sat41(h,s)
    persistent I J n T_star h_star s_star
    if isempty(T_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'Region4constants.mat'],'I_hs','J_hs','n_hs');
        I=vars.I_hs;
        J=vars.J_hs;
        n=vars.n_hs;
        T_star=550;
        h_star=2800*10^3;
        s_star=9.2*10^3;
    end
    
    
    if ~isempty(h)
        T=T_star.*sum(n.*(h./h_star-0.119).^I.*(s./s_star-1.07).^J,1);
    else
        T=NaN(1,0);
    end
end


function deltax=deltax(h,s,T)
    if T<=623.15
        p=IF97.p_satIntern(T);
        
        h0=IF97.h_pTx(p,T,0,@(x) x==1,true);
        h1=IF97.h_pTx(p,T,1,@(x) x==2,true);

        s0=IF97.s_pTx(p,T,0,@(x) x==1,true);
        s1=IF97.s_pTx(p,T,1,@(x) x==2,true);
    else
        rho0=IF97.rho_satLintern(T);
        rho1=IF97.rho_satVintern(T);
        
        h0=IF97.h_rhoTintern(rho0,T,@(x) x==3,rho0,rho1,0);
        h1=IF97.h_rhoTintern(rho1,T,@(x) x==3,rho0,rho1,1);

        s0=IF97.s_rhoTintern(rho0,T,@(x) x==3,rho0,rho1,0);
        s1=IF97.s_rhoTintern(rho1,T,@(x) x==3,rho0,rho1,1);
    end

    deltax=(h-h0)/(h1-h0)-(s-s0)/(s1-s0);
end




