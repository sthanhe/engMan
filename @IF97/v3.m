function v=v3(p,T)
    persistent I_v3 J_v3 n_v3 v_star p_star T_star a b c d e
    if isempty(v_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'Region3BackwardsConstants.mat'],'I_v3','J_v3','n_v3','v_star','p_star','T_star','a','b','c','d','e');
        I_v3=vars.I_v3;
        J_v3=vars.J_v3;
        n_v3=vars.n_v3;
        v_star=vars.v_star;
        p_star=vars.p_star;
        T_star=vars.T_star;
        a=vars.a;
        b=vars.b;
        c=vars.c;
        d=vars.d;
        e=vars.e;
    end
    
    
    if ~isempty(p)
        v=NaN(size(p));

        [regions,region14]=getT3regions(p,T);
        p=reshape(p,numel(p),1);
        T=reshape(T,numel(T),1);


        pi=p./p_star(regions)-a(regions);
        theta=T./T_star(regions)-b(regions);

        
        pi14=pi(region14)';
        theta14=theta(region14)';

        pi=(pi(~region14).^c(regions(~region14)))';
        theta=(theta(~region14).^d(regions(~region14)))';


        if any(region14)
            Iloc14=I_v3(:,regions(region14));
            Jloc14=J_v3(:,regions(region14));
            nloc14=n_v3(:,regions(region14));
            
            v(region14)=v_star(regions(region14))'.*exp(sum(nloc14.*pi14.^Iloc14.*theta14.^Jloc14,1));
        end
        
        if any(~region14)
            Iloc=I_v3(:,regions(~region14));
            Jloc=J_v3(:,regions(~region14));
            nloc=n_v3(:,regions(~region14));
            
            v(~region14)=v_star(regions(~region14))'.*sum(nloc.*pi.^Iloc.*theta.^Jloc,1).^(e(regions(~region14))');
        end
        

    else
        v=NaN(1,0);
    end
    
end


function [regions,region3n]=getT3regions(p,T)
    persistent I_T3 n_T3 p_star p_sat643 p_sat623
    if isempty(p_star)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'Region3BackwardsConstants.mat'],'I_T3','n_T3');
        I_T3=vars.I_T3;
        n_T3=vars.n_T3;
        p_star=10^6;
        
        p_sat643=IF97.p_satIntern(643.15);
        p_sat623=IF97.p_satIntern(623.15);
    end
    
    
    plim1=40*10^6<p;
    plim2=25*10^6<p;
    plim3=23.5*10^6<p;
    plim4=23*10^6<p;
    plim5=22.5*10^6<p;
    plim6=p_sat643<p;
    plim7=20.5*10^6<p;
    plim8=1.900881189173929*10^7<p;
    plim9=22.11*10^6<p;
    plim10=22.064*10^6<p;
    plim11=2.193161551*10^7<p;
    plim12=p_sat643<p;
    plim13=2.190096265*10^7<p;
    
    
    prange1=plim1 & p<=100*10^6;
    prange2=plim2 & ~plim1;
    prange3=plim3 & ~plim2;
    prange4=plim4 & ~plim3;
    prange5=plim5 & ~plim4;
    prange6=plim6 & ~plim5;
    prange7=plim7 & ~plim6;
    prange8=plim8 & ~plim7;
    prange9=p_sat623<p & ~plim8;    
    prange10=plim9 & ~plim5;
    prange11=plim10 & ~plim9;
    prange12=plim11 & ~plim10;
    prange14=plim13 & ~plim10;
    
    prange3_4=prange3 | prange4;
    prange4_5=prange4 | prange5;
    prange3_4_5=prange3 | prange4_5;
    prange10_11=prange10 | prange11;
    

    T_sat=IF97.T_satIntern(p);
    p=p./p_star;
    
    Trange3gh=T<=T_3gh;
    Trange3ij=T<=T_3ij;
    Trange3mn=T<=T_3mn;
    Trange3op=T<=T_3op;
    Trange3qu=T_3qu<T;
    Trange3uv=T_3uv<T;
    Trange3rx=T<=T_3rx;
    Trange3wx=T<=T_3wx;
    
    Trange3a=T<=T_3ab;
    Trange3c1=T<=T_3cd;
    Trange3c2=T<=T_sat;
    Trange3f=T_3ef<T;
    Trange3k=T_3jk<T;
    Trange3t=T_sat<=T;
    Trange3u=Trange3qu & ~Trange3uv;
    Trange3v=Trange3uv & ~Trange3f;
    Trange3w=Trange3f & Trange3wx;
    Trange3x=~Trange3wx & Trange3rx;
    
     
    region3a=prange1 & Trange3a;
    region3b=prange1 & ~Trange3a;
    region3c=(prange2 & Trange3c1) | ...
            (prange3 & Trange3c1) | ...
            (prange4 & Trange3c1) | ...
            (prange5 & Trange3c1) | ...
            (prange6 & Trange3c1) | ...
            (prange7 & Trange3c1) | ...
            (prange8 & Trange3c1) | ...
            (prange9 & Trange3c2);
    region3d=prange2 & ~Trange3c1 & Trange3a;
    region3e=prange2 & ~Trange3a & ~Trange3f;
    region3f=prange2 & Trange3f;
    region3g=prange3 & ~Trange3c1 & Trange3gh;
    region3h=prange3_4 & ~Trange3gh & ~Trange3f;
    region3i=prange3_4 & Trange3f & Trange3ij;
    region3j=prange3_4_5 & ~Trange3ij & ~Trange3k;
    region3k=(prange3_4_5 | prange6 | prange7) & Trange3k;
    region3l=prange4_5 & ~Trange3c1 & Trange3gh;
    region3m=prange5 & ~Trange3gh & Trange3mn;
    region3n=prange5 & ~Trange3mn & ~Trange3f;
    region3o=prange5 & Trange3f & Trange3op;
    region3p=prange5 & ~Trange3op & Trange3ij;
    region3q=prange6 & ~Trange3c1 & ~Trange3qu;
    region3r=(prange6 & ~Trange3rx & ~Trange3k) | (prange7 & Trange3t & ~Trange3k);
    region3s=(prange7 | prange8) & ~Trange3c1 & Trange3c2;
    region3t=(prange8 | prange9) & Trange3t;
    region3u=(prange10_11 & Trange3u) | (prange12 & Trange3u & Trange3c2) | (plim12 & ~plim11 & Trange3qu & Trange3c2);
    region3v=prange10 & Trange3v;
    region3w=prange10 & Trange3w;
    region3x=(prange10_11 & Trange3x) | (prange14 & Trange3x & Trange3t) | (plim12 & ~plim13 & Trange3rx & Trange3t);
    region3y=(prange11 & Trange3v) | (prange12 & Trange3uv & Trange3c2);
    region3z=(prange11 & Trange3w) | (prange14 & Trange3wx & Trange3t);
    
      
    regions=NaN(size(p));
    regions(region3a)=1;
    regions(region3b)=2;
    regions(region3c)=3;
    regions(region3d)=4;
    regions(region3e)=5;
    regions(region3f)=6;
    regions(region3g)=7;
    regions(region3h)=8;
    regions(region3i)=9;
    regions(region3j)=10;
    regions(region3k)=11;
    regions(region3l)=12;
    regions(region3m)=13;
    regions(region3n)=14;
    regions(region3o)=15;
    regions(region3p)=16;
    regions(region3q)=17;
    regions(region3r)=18;
    regions(region3s)=19;
    regions(region3t)=20;
    regions(region3u)=21;
    regions(region3v)=22;
    regions(region3w)=23;
    regions(region3x)=24;
    regions(region3y)=25;
    regions(region3z)=26;
    

    function T=T_3ab()
        T=T3fx2(I_T3(:,1),n_T3(:,1));
    end

    function T=T_3cd()
        T=T3fx1(I_T3(:,2),n_T3(:,2));
    end

    function T=T_3ef()
        T=3.727888004*(p-22.064)+647.096;
    end

    function T=T_3gh()
        T=T3fx1(I_T3(:,3),n_T3(:,3));
    end

    function T=T_3ij()
        T=T3fx1(I_T3(:,4),n_T3(:,4));
    end

    function T=T_3jk()
        T=T3fx1(I_T3(:,5),n_T3(:,5));
    end

    function T=T_3mn()
        T=T3fx1(I_T3(:,6),n_T3(:,6));
    end

    function T=T_3op()
        T=T3fx2(I_T3(:,7),n_T3(:,7));
    end

    function T=T_3qu()
        T=T3fx1(I_T3(:,8),n_T3(:,8));
    end

    function T=T_3rx()
        T=T3fx1(I_T3(:,9),n_T3(:,9));
    end

    function T=T_3uv()
        T=T3fx1(I_T3(:,10),n_T3(:,10));
    end

    function T=T_3wx()
        T=T3fx2(I_T3(:,11),n_T3(:,11));
    end


    function T=T3fx1(I,n)
        T=sum(n.*p.^I,1);
    end

    function T=T3fx2(I,n)
        T=sum(n.*log(p).^I,1);
    end

end

%     prange1=40*10^6<p & p<=100*10^6;
%     prange2=25*10^6<p & p<=40*10^6;
%     prange3=23.5*10^6<p & p<=25*10^6;
%     prange4=23*10^6<p & p<=23.5*10^6;
%     prange5=22.5*10^6<p & p<=23*10^6;
%     prange6=p_sat643<p & p<=22.5*10^6;
%     prange7=20.5*10^6<p & p<=p_sat643;
%     prange8=1.900881189173929*10^7<p & p<=20.5*10^6;
%     prange9=p_sat623<p & p<=1.900881189173929*10^7;
%     
%     prange10=22.11*10^6<p & p<=22.5*10^6;
%     prange11=22.064*10^6<p & p<=22.11*10^6;
%     prange12=2.193161551*10^7<p & p<=22.064*10^6;
%     prange13=p_sat643<p & p<=2.193161551*10^7;
%     prange14=2.190096265*10^7<p & p<=22.064*10^6;
%     prange15=p_sat643<p & p<=2.190096265*10^7;
% 
% %     1.Iteration:
% %     prange13=plim12 & ~plim11;
% %     prange15=plim12 & ~plim13;
% 
% 
% 
%     Trange3a=T<=T_3ab;
%     Trange3b=T_3ab<T;
%     Trange3c1=T<=T_3cd;
%     Trange3c2=T<=IF97.T_sat(p);
%     Trange3d=T_3cd<T & T<=T_3ab;
%     Trange3e=T_3ab<T & T<=T_3ef;
%     Trange3f=T_3ef<T;
%     Trange3g=T_3cd<T & T<=T_3gh;
%     Trange3h=T_3gh<T & T<=T_3ef;
%     Trange3i=T_3ef<T & T<=T_3ij;
%     Trange3j=T_3ij<T & T<=T_3jk;
%     Trange3k=T_3jk<T;
%     Trange3l=T_3cd<T & T<=T_3gh;
%     Trange3m=T_3gh<T & T<=T_3mn;
%     Trange3n=T_3mn<T & T<=T_3ef;
%     Trange3o=T_3ef<T & T<=T_3op;
%     Trange3p=T_3op<T & T<=T_3ij;
%     Trange3q=T_3cd<T & T<=T_3qu;
%     Trange3r1=T_3rx<T & T<=T_3jk;
%     Trange3r2=IF97.T_sat(p)<=T & T<=T_3jk;
%     Trange3s=T_3cd<T & T<=IF97.T_sat(p);
%     Trange3t=IF97.T_sat(p)<=T;
%     
%     Trange3uSup=T_3qu<T & T<=T_3uv;
%     Trange3vSup=T_3uv<T & T<=T_3ef;
%     Trange3wSup=T_3ef<T & T<=T_3wx;
%     Trange3xSup=T_3wx<T & T<=T_3rx;
%     Trange3ySup=T_3uv<T & T<=T_3ef;
%     Trange3zSup=T_3ef<T & T<=T_3wx;
%     
%     Trange3uSub1=T_3qu<T & T<=T_3uv & T<=IF97.T_sat(p);
%     Trange3uSub2=T_3qu<T & T<=IF97.T_sat(p);
%     Trange3xSub1=T_3wx<T & T<=T_3rx & IF97.T_sat(p)<=T;
%     Trange3xSub2=T<=T_3rx & IF97.T_sat(p)<=T;
%     Trange3ySub=T_3uv<T & T<=IF97.T_sat(p);
%     Trange3zSub=T<=T_3wx & IF97.T_sat(p)<=T;
% 
% 
% %     1. Iteration:
% % 
% %     Trange3b=~Trange3a;
% %     Trange3d=~Trange3c1 & Trange3a;
% %     Trange3e=~Trange3a & ~Trange3f;
% %     Trange3g=~Trange3c1 & Trange3gh;
% %     Trange3h=~Trange3gh & ~Trange3f;
% %     Trange3i=Trange3f & Trange3ij;
% %     Trange3j=~Trange3ij & ~Trange3k;
% %     Trange3l=~Trange3c1 & Trange3gh;
% %     Trange3m=~Trange3gh & Trange3mn;
% %     Trange3n=~Trange3mn & ~Trange3f;
% %     Trange3o=Trange3f & Trange3op;
% %     Trange3p=~Trange3op & Trange3ij;
% %     Trange3q=~Trange3c1 & T<=T_3qu;
% %     Trange3r1=T_3rx<T & ~Trange3k;
% %     Trange3r2=Trange3t & ~Trange3k;
% %     Trange3s=~Trange3c1 & Trange3c2;
% % 
% %     Trange3ySup=Trange3vSup;
% %     Trange3zSup=Trange3wSup;
% %     
% %     Trange3uSub1=Trange3uSup & Trange3c2;
% %     Trange3uSub2=Trange3qu & Trange3c2;
% %     Trange3xSub1=Trange3xSup & Trange3t;
% %     Trange3xSub2=Trange3rx & Trange3t;
% %     Trange3ySub=Trange3uv & Trange3c2;
% %     Trange3zSub=Trange3wx & Trange3t;
% 
% 
%     region3a=prange1 & Trange3a;
%     region3b=prange1 & Trange3b;
%     region3c=(prange2 & Trange3c1) | ...
%             (prange3 & Trange3c1) | ...
%             (prange4 & Trange3c1) | ...
%             (prange5 & Trange3c1) | ...
%             (prange6 & Trange3c1) | ...
%             (prange7 & Trange3c1) | ...
%             (prange8 & Trange3c1) | ...
%             (prange9 & Trange3c2);
%     region3d=prange2 & Trange3d;
%     region3e=prange2 & Trange3e;
%     region3f=prange2 & Trange3f;
%     region3g=prange3 & Trange3g;
%     region3h=(prange3 | prange4) & Trange3h;
%     region3i=(prange3 | prange4) & Trange3i;
%     region3j=(prange3 | prange4 | prange5) & Trange3j;
%     region3k=(prange3 | prange4 | prange5 | prange6 | prange7) & Trange3k;
%     region3l=(prange4 | prange5) & Trange3l;
%     region3m=prange5 & Trange3m;
%     region3n=prange5 & Trange3n;
%     region3o=prange5 & Trange3o;
%     region3p=prange5 & Trange3p;
%     region3q=prange6 & Trange3q;
%     region3r=(prange6 & Trange3r1) | (prange7 & Trange3r2);
%     region3s=(prange7 | prange8) & Trange3s;
%     region3t=(prange8 | prange9) & Trange3t;
%     region3u=((prange10 | prange11) & Trange3uSup) | (prange12 & Trange3uSub1) | (prange13 & Trange3uSub2);
%     region3v=prange10 & Trange3vSup;
%     region3w=prange10 & Trange3wSup;
%     region3x=((prange10 | prange11) & Trange3xSup) | (prange14 & Trange3xSub1) | (prange15 & Trange3xSub2);
%     region3y=(prange11 & Trange3ySup) | (prange12 & Trange3ySub);
%     region3z=(prange11 & Trange3zSup) | (prange14 & Trange3zSub);


