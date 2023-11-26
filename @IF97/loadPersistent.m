function loadPersistent(n)
    persistent Tp Trho Tph Tps p_satT T_satp vpT
    if isempty(Tp)
        folder=IF97.getFolder();
        vars=load([folder,filesep,'CompVerifConstants.mat'],'Tp','Trho','Tph','Tps','p_satT','T_satp','vpT');
        Tp=vars.Tp;
        Trho=vars.Trho;
        Tph=vars.Tph;
        Tps=vars.Tps;
        p_satT=vars.p_satT;
        T_satp=vars.T_satp;
        vpT=vars.vpT;
    end
    
    if nargin<1
        n=1;
    else
        tic;
    end
    
    ME=@(x) MException('IF97:verficationFailed','Test of function(s) %s failed',x);
    
    
    p=repmat([3,80,3,0.0035,0.0035,30,0.5,30,30]*10^6,1,1,n);
    T=repmat([300,300,500,300,700,700,1500,1500,2000],1,1,n);
    test=norm(vertcat(IF97.v(p,T),IF97.h(p,T),IF97.u(p,T),IF97.s(p,T),IF97.c_p(p,T),IF97.w(p,T)));
    err=abs(test-repmat(norm(Tp),1,1,n));
    if any(err>10^-9)
        if n>1
            toc;
        end
        throw(ME('v(p,T), h(p,T), u(p,T), s(p,T), c_p(p,T) or w(p,T)'));
    end
    
    
    if n>1
        T=repmat(p_satT(:,1),1,n);
        test=norm(IF97.p_sat(T));
        err=abs(test-repmat(norm(p_satT(:,2)),1,n));
        if any(err>10^-9)
            toc;
            throw(ME('p_sat(T)'));
        end


        p=repmat(T_satp(:,1),1,n);
        test=norm(IF97.T_sat(p));
        err=abs(test-repmat(norm(T_satp(:,2)),1,n));
        if any(err>10^-9)
            toc;
            throw(ME('T_sat(p)'));
        end
        
        
        p=repmat(vpT(:,1),1,n);
        T=repmat(vpT(:,2),1,n);
        test=norm(IF97.v(p,T,NaN));
        err=abs(test-repmat(norm(vpT(:,3)),1,n));
        if any(err>10^-9)
            toc;
            throw(ME('v3(p,T)'));
        end
    end
    
    
    if n>1
        toc;
    end
    

end


function X=norm(X)
    X=X.*10.^(-floor(log10(X))-1);
end




