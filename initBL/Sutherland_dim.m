function [mu,kappa] = Sutherland_dim(BL_PARA)

global visc
    
    T_inf=BL_PARA.T_inf;
    T_ref=BL_PARA.T_ref;
    mu_ref=BL_PARA.mu_ref;
    kappa_ref=BL_PARA.kappa_ref;
    
    if strcmp(visc,'PowerLaw')

        powerexp=BL_PARA.powerexp;
        mu_ref=BL_PARA.mu_ref;
        kappa_ref=BL_PARA.kappa_ref;

        mu = mu_ref * (T_inf/T_ref).^(powerexp);
        kappa= kappa_ref*(T_inf/T_ref).^(powerexp);

    elseif strcmp(visc,'Sutherland')

        S1=BL_PARA.S_mu_ref;
        S2=BL_PARA.S_kappa_ref;
        
        mu = mu_ref * (T_ref + S1) ./ (T_inf + S1) .* (T_inf / T_ref).^(3.0 / 2.0);
        kappa= kappa_ref*(T_ref + S2) ./ (T_inf + S2).* (T_inf / T_ref).^(3.0 / 2.0);

    end
    
end

