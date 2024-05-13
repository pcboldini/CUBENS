function [mu,mu_T,mu_TT] = Sutherland(T)

global visc

    BL_PARA=Set_BL_Para_IG;

    if strcmp(visc,'Constant')
    
        mu=1;
        mu_T=0;
        mu_TT=0;
    
    elseif strcmp(visc,'PowerLaw')

        powerexp=BL_PARA.powerexp;
        mu=T^powerexp;
        mu_T=powerexp*T^(powerexp-1);
        mu_TT=powerexp*(powerexp-1)*T^(powerexp-2);

    elseif strcmp(visc,'Sutherland')
        
        aa=BL_PARA.S_mu_ref/BL_PARA.T_inf;

        mu = (1.0+aa)*T.^1.5./(T+aa);
        mu_T = (1.0+aa)*sqrt(T).*(0.5*T+1.5*aa)./((T+aa).^2.0);
        mu_TT= (0.75*(T.^0.5+aa*T.^(-0.5)).*(T+aa).^2.0-...
        (T.^1.5+3.0*aa*T.^0.5).*(T+aa)).*(1.0+aa)./(T+aa).^4.0;
    
    end
    
end

