%% Inner heat transfer coefficient
%By Theresa Brunauer
%
%Part of the sandTES Engineering Manual
%
%All required files for this function can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.xxxxxxx
% 
%All parameters and results are in SI base units.
%
%
%
%This function calculate the inner heat transfer coefficient according to
%Gnielinski. It is a slightly modified copy of the same function in the
%heat exchanger calculation class (see the documentation there):
%
% Link Github
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%Necessary files, classes, functions, and scripts:
%   - @IF97


function alpha_1=alpha_in(Rho1,T1m,mdot1,D,N,x)
    persistent g 
    if isempty(g)
        g = 9.80665; %acceleration of gravity
    end
    len=numel(Rho1);
    
    mdot1=repmat(mdot1,1,len);
    D=repmat(D,1,len);
    N=repmat(N,1,len);

    x(x<D) = D(x<D); 

    is2phase = IF97.is2phase(Rho1,T1m);
    alphainside = NaN(1,len); 

    a = find(is2phase==0);

    A = D.^2.*pi/4;
    mDot = mdot1./(A.*N);

    if a > 0
        eta = IF97.my(Rho1,T1m); % viscosity
        lambda1 = IF97.lambda(Rho1,T1m); %heat conductivity
        u = mDot./Rho1; %velocity
        Pr = IF97.Pr(Rho1,T1m); % Pr-number
        Re = u.*D.*Rho1./eta; % Re-number
    
        %%% SINGL PHASE
        % G1-3 Wärmeübertragung bei laminarer Strömung durch Rohre (S.787)
        % G1-4 Wärmeübertragung bei turbulenter Strömung durch Rohre (S.788)

        Nu_x = NaN(1, len); %local Nu-number
        lam = Re < 2300; % laminar flow
        Nu_x2 = 1.302 .* (Re.*Pr.*D./x).^(1/3);
        Nu_x(lam) = (4.364.^3 + 1 + (Nu_x2(lam) - 1).^3).^(1./3);
        tur = Re > 10^4; % turulent flow
        Xi = (1.8 .* log10(Re)-1.5).^(-2);
        Nu_x(tur) = (Xi(tur)./8).*(Re(tur)-10^3).*Pr(tur)./(1+12.7.*(Xi(tur)./8).^0.5.*(Pr(tur)-1)).*(1+1/3.*(D(tur)./x(tur)).^(2/3)); %neue Formel

        tra = 2300 < Re & Re < 10^4; %transition region
        Nu_x2L = 1.302 .* (2300.*Pr(tra).*D(tra)./x(tra)).^(1./3);
        Nu_x3L = 0.462*Pr(tra).^(1/3).*(2300.*D(tra)./x(tra)).^(1/2);
        Nu_xL = (4.364^3+(Nu_x2L-0.6).^3+(Nu_x3L).^3).^(1/3);
        Nu_xT = (0.0308./8).*10^4.*Pr(tra)./(1+12.7.*(0.0308./8).^(1/2).*(Pr(tra).^(2./3)-1)).*(1+(1/3).*(D(tra)./x(tra)).^(2/3));
        Gamma = (Re(tra) - 2300)./(10.^4 -2300);
        Nu_x(tra) = (1-Gamma) .* Nu_xL + Gamma .* Nu_xT;

        alphainside(~is2phase) = Nu_x(~is2phase) .* lambda1(~is2phase) ./ D(~is2phase); % heat transfer coefficient
    end

    % % H3.4.2 Blasensieden reiner Stoffe in durchströmten Rohren (S.919)
    % persistent cf q0 d0 Ra0 pc alpha_0                    
    % if isempty(cf)
    %     cf = IF97.cf;    %properties of the fluid                                                    
    %     q0 = IF97.q0;  %substance-IF97ific values
    %     d0 = IF97.d0;
    %     Ra0 = IF97.Ra0;
    %     pc = IF97.p_c;
    %     alpha_0 = IF97.alpha_0;
    % end
    % Ra_in = 0.5 * 10^-6; %Angenommen
    % p_st = P1./pc;
    % n = 0.8 - 0.1.*10.^(0.76 .* p_st);
    % 
    % q = abs(QDot1(2:end)./(D.^2*pi./4));
    % 
    % A = NaN(1,len);
    % B = NaN(1,len);
    % C = NaN(1,len);
    % 
    % A(is2phase) = cf.*(q(is2phase)./q0).^n(is2phase);
    % B(is2phase) = 2.816 .* p_st(is2phase) .^0.45 + (3.4 + 1.7./(1-p_st(is2phase).^7)).*p_st(is2phase).^3.7;
    % C(is2phase) = (d0./D(is2phase)).^0.4 .* (Ra_in./Ra0).^0.133;
    % 
    % alphainside(is2phase) = A(is2phase).*B(is2phase).*C(is2phase).*alpha_0; % heat transfer coefficient 
    alphainside(is2phase)=1e6;
    alpha_1 = alphainside;
end




