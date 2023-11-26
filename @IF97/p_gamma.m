function p=p_gamma(rho,T,regions,region2,gammaregion)
    sz=size(rho);


    regind=NaN(sz);
    regind(regions(1))=1;
    regind(regions(2))=2;
    regind(regions(5))=5;
    regind=regind(gammaregion);


    p_lowBound=NaN(sz);
    p_upBound=NaN(sz);

    p_lowBound(regions(1))=IF97.p_satIntern(T(regions(1)));
    p_upBound(regions(1))=repmat(100*10^6,1,nnz(regions(1)));

    p_lowBound(regions(2))=realmin;
    p_upBound(region2(1))=IF97.p_satIntern(T(region2(1)));
    p_upBound(region2(2))=IF97.B23(T(region2(2)),true);
    p_upBound(region2(3))=repmat(100*10^6,1,nnz(region2(3)));
    
    p_lowBound(regions(5))=realmin;
    p_upBound(regions(5))=repmat(50*10^6,1,nnz(regions(5)));


    p_lowBound=p_lowBound(gammaregion);
    p_upBound=p_upBound(gammaregion);
    rho=rho(gammaregion);
    T=T(gammaregion);
    
    
    p=NaN(size(rho));
    for i=1:length(p)
        p(i)=fzero(@(p) vSwitch(p,T(i),regind(i)).^-1-rho(i),[p_lowBound(i),p_upBound(i)]);
    end
end


function v=vSwitch(p,T,reg)
    switch reg
        case 1
            [pi,tau]=IF97.pitau1(p,T);
            gamma=@(dpi,dtau) IF97.gamma1(pi,tau,dpi,dtau);
        case 2
            [pi,tau]=IF97.pitau2(p,T);
            gamma=@(dpi,dtau) IF97.gamma2(pi,tau,dpi,dtau);
        case 5
            [pi,tau]=IF97.pitau5(p,T);
            gamma=@(dpi,dtau) IF97.gamma5(pi,tau,dpi,dtau);
    end
    
    v=IF97.v_pTx(p,T,pi,gamma);
end




