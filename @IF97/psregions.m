function [regions,region2,region3]=psregions(p,s,sz)
    persistent p_sat273 p_sat623
    if isempty(p_sat273)
        p_sat273=IF97.p_satIntern(273.15);
        p_sat623=IF97.p_satIntern(623.15);
    end
    
    
    T_sat=IF97.T_sat(p);
    s0=IF97.s(p,T_sat,0);
    s1=IF97.s(p,T_sat,1);
    
    
    szp=size(p);
    p=reshape(p,1,numel(p));
    
    
    smin=NaN(szp);
    reg=p_sat273<=p & p<=100*10^6;
    smin(reg)=sreg(273.15,1);
    reg=0<p & p<p_sat273;
    smin(reg)=sreg(273.15,2);
    
    
    s13=NaN(szp);
    reg=p_sat623<=p & p<=100*10^6;
    s13(reg)=sreg(623.15,1);
    
    
    s23=NaN(szp);
    s23(reg)=sreg(IF97.B23(p(reg),false),2);
    
    
    smax2=NaN(szp);
    reg=0<p & p<=100*10^6;
    smax2(reg)=sreg(1073.15,2);
    
    
    smax5=NaN(szp);
    reg=0<p & p<=50*10^6;
    smax5(reg)=sreg(2273.15,5);


    p=reshape(p,szp);
    page=IF97.pagefx(sz);
    regionmx=false([sz,5]);

    regionmx(page(1))=(p_sat273<=p & p<=p_sat623 & smin<=s & s<s0) | ... 
                        (p_sat623<p & p<=100*10^6 & smin<=s & s<=s13);
    regionmx(page(2))=(0<p & p<p_sat273 & smin<=s & s<=smax2) | ...
                    (p_sat273<=p & p<=p_sat623 & s1<s & s<=smax2) | ...
                    (p_sat623<p & p<=100*10^6 & s23<=s & s<=smax2);
    regionmx(page(3))=(p_sat623<p & p<=IF97.p_c & s13<s & s<s0) | ...
                    (p_sat623<p & p<=IF97.p_c & s1<s & s<s23) | ...
                    (IF97.p_c<p & p<=100*10^6 & s13<s & s<s23);
    regionmx(page(4))=p_sat273<=p & p<=IF97.p_c & s0<=s & s<=s1;
    regionmx(page(5))=0<p & p<=50*10^6 & smax2<s & s<=smax5;
    
    regions=@(x) regionmx(page(x));
    
    
    region2mx=false([sz,3]);
    region2mx(page(1))=reshape(regions(2),sz) & p<=4*10^6 & true(size(s));
    region2mx(page(3))=reshape(regions(2),sz) & 6.5467*10^6<p & s<5.85*10^3;
    region2mx(page(2))=regions(2) & ~region2mx(page(1)) & ~region2mx(page(3));
    
    region2=@(x) region2mx(page(x));
    
    
    region3mx=false([sz,2]);
    region3mx(page(1))=(reshape(regions(3),sz) & p<=IF97.p_c & s<s0) | ...
                        (reshape(regions(3),sz) & IF97.p_c<p & s<=IF97.s_c);
    region3mx(page(2))=regions(3) & ~region3mx(page(1));
    
    region3=@(x) region3mx(page(x));
    
    
    function s=sreg(T,i)
        if ~isempty(T)
            n=nnz(reg);
            exp=@(x) repmat(x,1,n./numel(x));
            s=IF97.s_pTx(p(reg),exp(T),NaN(1,n),@(x) exp(x==i),true(1,n));
        else
            s=NaN(1,0);
        end
    end
end




