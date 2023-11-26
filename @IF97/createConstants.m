function createConstants()
    clear('functions'); %#ok<CLFUNC>
    folder=IF97.getFolder();
    file=[folder,filesep,'constants.xls'];

    
    %% Boundaries
    boundaries=xlsread(file,'Boundaries','B2:C6');
    n_B23=boundaries(:,1);
    n_B2bc=boundaries(:,2);
    save([folder,filesep,'BoundaryConstants.mat'],'n_B23','n_B2bc');
    
    boundaries=xlsread(file,'Boundaries','D2:F5');
    n_B3ab=boundaries(:,1);
    n_2ab=boundaries(:,2);
    n_3ab=boundaries(:,3);
    save([folder,filesep,'BoundaryConstants.mat'],'n_B3ab','n_2ab','n_3ab','-append');
    
    boundaries=xlsread(file,'Boundaries','G2:I28');
    I_h1=boundaries(:,1);
    J_h1=boundaries(:,2);
    n_h1=boundaries(:,3);
    save([folder,filesep,'BoundaryConstants.mat'],'I_h1','J_h1','n_h1','-append');
    
    boundaries=xlsread(file,'Boundaries','J2:L20');
    I_h3a=boundaries(:,1);
    J_h3a=boundaries(:,2);
    n_h3a=boundaries(:,3);
    save([folder,filesep,'BoundaryConstants.mat'],'I_h3a','J_h3a','n_h3a','-append');
    
    boundaries=xlsread(file,'Boundaries','M2:O31');
    I_h2ab_satV=boundaries(:,1);
    J_h2ab_satV=boundaries(:,2);
    n_h2ab_satV=boundaries(:,3);
    save([folder,filesep,'BoundaryConstants.mat'],'I_h2ab_satV','J_h2ab_satV','n_h2ab_satV','-append');
    
    boundaries=xlsread(file,'Boundaries','P2:R17');
    I_h2c3b=boundaries(:,1);
    J_h2c3b=boundaries(:,2);
    n_h2c3b=boundaries(:,3);
    save([folder,filesep,'BoundaryConstants.mat'],'I_h2c3b','J_h2c3b','n_h2c3b','-append');
    
    boundaries=xlsread(file,'Boundaries','S2:U7');
    I_hB13=boundaries(:,1);
    J_hB13=boundaries(:,2);
    n_hB13=boundaries(:,3);
    save([folder,filesep,'BoundaryConstants.mat'],'I_hB13','J_hB13','n_hB13','-append');
    
    boundaries=xlsread(file,'Boundaries','V2:X26');
    I_TB23=boundaries(:,1);
    J_TB23=boundaries(:,2);
    n_TB23=boundaries(:,3);
    save([folder,filesep,'BoundaryConstants.mat'],'I_TB23','J_TB23','n_TB23','-append');
    

    %% Region1
    region1_pT=xlsread(file,'Region1','B2:D35');
    I_pT=region1_pT(:,1);
    J_pT=region1_pT(:,2);
    n_pT=region1_pT(:,3);
    save([folder,filesep,'Region1constants.mat'],'I_pT','J_pT','n_pT');

    region1_ph=xlsread(file,'Region1','E2:G21');
    I_ph=region1_ph(:,1);
    J_ph=region1_ph(:,2);
    n_ph=region1_ph(:,3);
    save([folder,filesep,'Region1constants.mat'],'I_ph','J_ph','n_ph','-append');

    region1_ps=xlsread(file,'Region1','H2:J21');
    I_ps=region1_ps(:,1);
    J_ps=region1_ps(:,2);
    n_ps=region1_ps(:,3);
    save([folder,filesep,'Region1constants.mat'],'I_ps','J_ps','n_ps','-append');
    
    region1_hs=xlsread(file,'Region1','K2:M20');
    I_hs=region1_hs(:,1);
    J_hs=region1_hs(:,2);
    n_hs=region1_hs(:,3);
    save([folder,filesep,'Region1constants.mat'],'I_hs','J_hs','n_hs','-append');


    %% Region2
    region2_pT0=xlsread(file,'Region2','B2:C10');
    J_pT0=region2_pT0(:,1);
    n_pT0=region2_pT0(:,2);
    save([folder,filesep,'Region2constants.mat'],'J_pT0','n_pT0');

    region2_pTr=xlsread(file,'Region2','D2:F44');
    I_pTr=region2_pTr(:,1);
    J_pTr=region2_pTr(:,2);
    n_pTr=region2_pTr(:,3);
    save([folder,filesep,'Region2constants.mat'],'I_pTr','J_pTr','n_pTr','-append');

    region2_pha=xlsread(file,'Region2','G2:I35');
    I_pha=region2_pha(:,1);
    J_pha=region2_pha(:,2);
    n_pha=region2_pha(:,3);
    save([folder,filesep,'Region2constants.mat'],'I_pha','J_pha','n_pha','-append');

    region2_phb=xlsread(file,'Region2','J2:L39');
    I_phb=region2_phb(:,1);
    J_phb=region2_phb(:,2);
    n_phb=region2_phb(:,3);
    save([folder,filesep,'Region2constants.mat'],'I_phb','J_phb','n_phb','-append');

    region2_phc=xlsread(file,'Region2','M2:O24');
    I_phc=region2_phc(:,1);
    J_phc=region2_phc(:,2);
    n_phc=region2_phc(:,3);
    save([folder,filesep,'Region2constants.mat'],'I_phc','J_phc','n_phc','-append');

    region2_psa=xlsread(file,'Region2','P2:R47');
    I_psa=region2_psa(:,1);
    J_psa=region2_psa(:,2);
    n_psa=region2_psa(:,3);
    save([folder,filesep,'Region2constants.mat'],'I_psa','J_psa','n_psa','-append');

    region2_psb=xlsread(file,'Region2','S2:U45');
    I_psb=region2_psb(:,1);
    J_psb=region2_psb(:,2);
    n_psb=region2_psb(:,3);
    save([folder,filesep,'Region2constants.mat'],'I_psb','J_psb','n_psb','-append');

    region2_psc=xlsread(file,'Region2','V2:X31');
    I_psc=region2_psc(:,1);
    J_psc=region2_psc(:,2);
    n_psc=region2_psc(:,3);
    save([folder,filesep,'Region2constants.mat'],'I_psc','J_psc','n_psc','-append');
    
    region2_hsa=xlsread(file,'Region2','Y2:AA30');
    I_hsa=region2_hsa(:,1);
    J_hsa=region2_hsa(:,2);
    n_hsa=region2_hsa(:,3);
    save([folder,filesep,'Region2constants.mat'],'I_hsa','J_hsa','n_hsa','-append');

    region2_hsb=xlsread(file,'Region2','AB2:AD34');
    I_hsb=region2_hsb(:,1);
    J_hsb=region2_hsb(:,2);
    n_hsb=region2_hsb(:,3);
    save([folder,filesep,'Region2constants.mat'],'I_hsb','J_hsb','n_hsb','-append');

    region2_hsc=xlsread(file,'Region2','AE2:AG32');
    I_hsc=region2_hsc(:,1);
    J_hsc=region2_hsc(:,2);
    n_hsc=region2_hsc(:,3);
    save([folder,filesep,'Region2constants.mat'],'I_hsc','J_hsc','n_hsc','-append');


    %% Region3
    region3=xlsread(file,'Region3','B2:D41');
    I=region3(:,1);
    J=region3(:,2);
    n=region3(:,3);
    save([folder,filesep,'Region3constants.mat'],'I','J','n');
    
    region3_pha=xlsread(file,'Region3','E2:G32');
    I_pha=region3_pha(:,1);
    J_pha=region3_pha(:,2);
    n_pha=region3_pha(:,3);
    save([folder,filesep,'Region3constants.mat'],'I_pha','J_pha','n_pha','-append');

    region3_phb=xlsread(file,'Region3','H2:J34');
    I_phb=region3_phb(:,1);
    J_phb=region3_phb(:,2);
    n_phb=region3_phb(:,3);
    save([folder,filesep,'Region3constants.mat'],'I_phb','J_phb','n_phb','-append');
    
    region3_psa=xlsread(file,'Region3','K2:M34');
    I_psa=region3_psa(:,1);
    J_psa=region3_psa(:,2);
    n_psa=region3_psa(:,3);
    save([folder,filesep,'Region3constants.mat'],'I_psa','J_psa','n_psa','-append');

    region3_psb=xlsread(file,'Region3','N2:P29');
    I_psb=region3_psb(:,1);
    J_psb=region3_psb(:,2);
    n_psb=region3_psb(:,3);
    save([folder,filesep,'Region3constants.mat'],'I_psb','J_psb','n_psb','-append');
    
    region3_hsa=xlsread(file,'Region3','Q2:S34');
    I_hsa=region3_hsa(:,1);
    J_hsa=region3_hsa(:,2);
    n_hsa=region3_hsa(:,3);
    save([folder,filesep,'Region3constants.mat'],'I_hsa','J_hsa','n_hsa','-append');

    region3_hsb=xlsread(file,'Region3','T2:V36');
    I_hsb=region3_hsb(:,1);
    J_hsb=region3_hsb(:,2);
    n_hsb=region3_hsb(:,3);
    save([folder,filesep,'Region3constants.mat'],'I_hsb','J_hsb','n_hsb','-append');
    
    
    %% Region3Backwards
    region3_redQuan=xlsread(file,'Region3Backwards','B2:J27');
    v_star=region3_redQuan(:,1);
    p_star=region3_redQuan(:,2);
    T_star=region3_redQuan(:,3);
    a=region3_redQuan(:,5);
    b=region3_redQuan(:,6);
    c=region3_redQuan(:,7);
    d=region3_redQuan(:,8);
    e=region3_redQuan(:,9);
    save([folder,filesep,'Region3BackwardsConstants.mat'],'v_star','p_star','T_star','a','b','c','d','e');
    
    I_T3=xlsread(file,'Region3Backwards','B31:L35');
    save([folder,filesep,'Region3BackwardsConstants.mat'],'I_T3','-append');
    
    n_T3=xlsread(file,'Region3Backwards','B39:L43');
    save([folder,filesep,'Region3BackwardsConstants.mat'],'n_T3','-append');
    
    I_v3=xlsread(file,'Region3Backwards','B47:AA89');
    save([folder,filesep,'Region3BackwardsConstants.mat'],'I_v3','-append');
    
    J_v3=xlsread(file,'Region3Backwards','B93:AA135');
    save([folder,filesep,'Region3BackwardsConstants.mat'],'J_v3','-append');
    
    n_v3=xlsread(file,'Region3Backwards','B139:AA181');
    save([folder,filesep,'Region3BackwardsConstants.mat'],'n_v3','-append');
    

    %% Region4
    n=xlsread(file,'Region4','B2:B11');
    save([folder,filesep,'Region4constants.mat'],'n');
    
    coeff=xlsread(file,'Region4','C2:F7');
    save([folder,filesep,'Region4constants.mat'],'coeff','-append');
    
    region4_hs=xlsread(file,'Region4','G2:I37');
    I_hs=region4_hs(:,1);
    J_hs=region4_hs(:,2);
    n_hs=region4_hs(:,3);
    save([folder,filesep,'Region4constants.mat'],'I_hs','J_hs','n_hs','-append');


    %% Region5
    region5=xlsread(file,'Region5','B2:F7');
    J_0=region5(:,1);
    n_0=region5(:,2);
    I_r=region5(:,3);
    J_r=region5(:,4);
    n_r=region5(:,5);
    save([folder,filesep,'Region5constants.mat'],'J_0','n_0','I_r','J_r','n_r');
    
    
    %% Melting and Sublimation
    Ih=xlsread(file,'MeltSubl','B2:C4');
    III=xlsread(file,'MeltSubl','B8:C8');
    V=xlsread(file,'MeltSubl','B12:C12');
    VI=xlsread(file,'MeltSubl','B16:C16');
    VII=xlsread(file,'MeltSubl','B20:C22');
    triple=xlsread(file,'MeltSubl','B26:C29');
    
    a=struct('Ih',Ih(:,1),'III',III(:,1),'V',V(:,1),'VI',VI(:,1),'VII',VII(:,1));
    b=struct('Ih',Ih(:,2),'III',III(:,2),'V',V(:,2),'VI',VI(:,2),'VII',VII(:,2));
    T_t=triple(:,1);
    p_t=triple(:,2);
    save([folder,filesep,'MeltSublConstants.mat'],'a','b','T_t','p_t');
    
    
    sub=xlsread(file,'MeltSubl','B33:C35');
    aSub=sub(:,1);
    bSub=sub(:,2);
    save([folder,filesep,'MeltSublConstants.mat'],'aSub','bSub','-append');
    
    
    %% Viscosity
    H0=xlsread(file,'Viscosity','B2:B5');
    H1=xlsread(file,'Viscosity','B10:H15');
    H1=reshape(H1,6,1,7);
    
    save([folder,filesep,'ViscosityConstants.mat'],'H0','H1');
    
    
    %% Thermal Conductivity
    L0=xlsread(file,'ThermCond','B2:B6');
    L1=xlsread(file,'ThermCond','B11:G15');
    L1=reshape(L1,5,1,6);
    A=xlsread(file,'ThermCond','B19:F24');
    
    save([folder,filesep,'ThermCondConstants.mat'],'L0','L1','A');
    
    
    %% CompVerif
    Tp=xlsread(file,'CompVerif','B2:J7');
    save([folder,filesep,'CompVerifConstants.mat'],'Tp');
    
    Trho=xlsread(file,'CompVerif','B11:D16');
    save([folder,filesep,'CompVerifConstants.mat'],'Trho','-append');
    
    Tph=xlsread(file,'CompVerif','A20:C31');
    save([folder,filesep,'CompVerifConstants.mat'],'Tph','-append');
    
    Tps=xlsread(file,'CompVerif','A35:C46');
    save([folder,filesep,'CompVerifConstants.mat'],'Tps','-append');
    
    p_satT=xlsread(file,'CompVerif','A50:B52');
    save([folder,filesep,'CompVerifConstants.mat'],'p_satT','-append');
    
    T_satp=xlsread(file,'CompVerif','A56:B58');
    save([folder,filesep,'CompVerifConstants.mat'],'T_satp','-append');
    
    vpT=xlsread(file,'CompVerif','B62:D113');
    save([folder,filesep,'CompVerifConstants.mat'],'vpT','-append');
end




