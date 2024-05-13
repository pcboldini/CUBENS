    function [LST_EIGS]=Postproc_LST_Spectrum(BL_PARA,LST_PARA,LST_MESH,LST_SPECTRUM,FLOW_LST,PRESSURE_LST,THERMO_LST,ENERGY_LST);

    global dimomega USE_EOS show_plot LST
    
    % Define parameters
    M_inf=BL_PARA.M_inf;
    if strcmp(USE_EOS,'IG')
        Rg=BL_PARA.Rg;
    end
    a_local=FLOW_LST(:,11);
    p=PRESSURE_LST(:,1);
    p_T=PRESSURE_LST(:,2);
    p_Rho=PRESSURE_LST(:,4);
    
    Spatial=LST_PARA.Spatial;
    Full = LST_PARA.Full;
    
    N_LST=LST_MESH.N_LST;
    y_LST=LST_MESH.y_LST;    
    Eval=LST_SPECTRUM.Eval;
    Evec=LST_SPECTRUM.Evec;
    Re0=LST_PARA.Re0;
    T_LST=FLOW_LST(:,7);
    Rho_LST=FLOW_LST(:,1);
    Rho_y_LST=FLOW_LST(:,2);
   
    U_LST=FLOW_LST(:,4);
    U_y_LST=FLOW_LST(:,5);
    U_yy_LST=FLOW_LST(:,6);
    
        
    if strcmp(Spatial,'true')
    if strcmp(dimomega,'yes')
        omega=LST_PARA.omega_0;  
    else
        omega     = LST_PARA.F_0*Re0;  
    end
    else 
        alpha=LST_PARA.alpha_0;
    end

    beta=LST_PARA.beta_0; 

    
%% Calculation

    cph_fast=1+1/M_inf;
    cph_slow=1-1/M_inf;
    
  % Inflection point
  
  IP_profile=Rho_y_LST.*U_y_LST+Rho_LST.*U_yy_LST;
     
    if strcmp(USE_EOS,'IG')
        count=5;   
        for i=1:N_LST
        rhoi_LST(i)=1+(i-1)*count;
        ui_LST(i)=2+(i-1)*count;
        vi_LST(i)=3+(i-1)*count;
        wi_LST(i)=4+(i-1)*count;
        Ti_LST(i)=5+(i-1)*count;
        end
    else
        count=6;   
        for i=1:N_LST
        rhoi_LST(i)=1+(i-1)*count;
        ui_LST(i)=2+(i-1)*count;
        vi_LST(i)=3+(i-1)*count;
        wi_LST(i)=4+(i-1)*count;
        Ti_LST(i)=5+(i-1)*count;
        pi_LST(i)=6+(i-1)*count;
        end
    end

    
    for index=1:numel(Eval)
        if strcmp(Spatial,'true')
        c_r(index)=omega./real(Eval(index));
        else
        c_r(index)=real(Eval(index))./real(alpha);
        end
    end
    
    % Entropy waves
    index_ew=find(abs(c_r-1)<0.001);
    ew=Eval(index_ew);
    if strcmp(Spatial,'true')
        c_r_ew=omega./real(ew);
    else
        c_r_ew=real(ew)./alpha;
    end

    % Slow waves
    index_sw=find(c_r<=cph_slow & c_r>0);
    sw=Eval(index_sw);
    if strcmp(Spatial,'true')
        c_r_sw=omega./real(sw);
    else
        c_r_sw=real(sw)./alpha;    
    end

    % Fast waves
    index_fw=find(c_r>=cph_fast);
    fw=Eval(index_fw);
    if strcmp(Spatial,'true')
        c_r_fw=omega./real(fw);
    else
        c_r_fw=real(fw)./alpha;    
    end
    
    if strcmp(Spatial,'true')               
    
   figure(1);
    
    plot(squeeze(real(Eval)), squeeze(imag(Eval)), 'kx','Linewidth',2)
    hold on
    plot(c_r_ew, imag(ew), 'yo','Linewidth',2); 
    plot(c_r_sw, imag(sw), 'ro','Linewidth',2);     
    plot(c_r_fw, imag(fw), 'bo','Linewidth',2); 
         if strcmp(USE_EOS,'IG')
            for j=1:count*N_LST, text(real(Eval(j)), imag(Eval(j)), [' ' num2str(j)]), end
         else
            for j=1:(2*count)*N_LST, text(real(Eval(j)), imag(Eval(j)), [' ' num2str(j)]), end
         end
    xlim([omega-0.3 omega+0.3]), ylim([-0.02 0.015])
    ylabel('$\alpha_i$'); xlabel('$\alpha_r$');
    reply = 'C';
    selectedEVi = 0;
    disp(' ');
    disp(' ');
    j = input('Index of eigenvector to plot (i.e. "59"+ENTER): ');

        selectedEVi = j;
        disp(['ALPHA(' num2str(selectedEVi) ') = ' num2str(Eval(selectedEVi)),', C(' num2str(selectedEVi) ') = ' num2str(omega./Eval(selectedEVi))]);
        eigV = Eval(selectedEVi);
        eigVec = Evec(:,selectedEVi);
        
        k_ph=sqrt(real(eigV)^2+beta^2);
        cr_eigV=omega/k_ph;
        disp(['c_r(' num2str(selectedEVi) ') = ' num2str(cr_eigV)]);

        
    else
        
    figure(1);
    
    plot(squeeze(real(Eval)), squeeze(imag(Eval)), 'kx','Linewidth',2)
    hold on
     plot(c_r_ew, imag(ew), 'yo','Linewidth',2); 
     plot(c_r_sw, imag(sw), 'ro','Linewidth',2);     
     plot(c_r_fw, imag(fw), 'bo','Linewidth',2); 
    
                 for j=1:5*N_LST, text(real(Eval(j)), imag(Eval(j)), [' ' num2str(j)]), end
                 
    xlim([1-0.8 1+0.8]), ylim([-0.02 0.005])
    ylabel('\omega_i'); xlabel('\omega_r');
    reply = 'C';
    selectedEVi = 0;
    disp(' ');
    disp(' ');
    j = input('Index of eigenvector to plot (i.e. "59"+ENTER): ');
        
        selectedEVi = j;
        disp(['OMEGA(' num2str(selectedEVi) ') = ' num2str(Eval(selectedEVi)),', C(' num2str(selectedEVi) ') = ' num2str(Eval(selectedEVi)/alpha)]);
        eigVec = Evec(:,j);
    end

        if strcmp(Spatial,'true')
        cr_new=omega./real(Eval(selectedEVi));
        else
        cr_new=real(Eval(selectedEVi))./alpha;
        end
        
    close;     
           
% Eigenfunctions        
        
        rho_EF=eigVec(rhoi_LST);
        u_EF=eigVec(ui_LST);
        v_EF=eigVec(vi_LST);
        w_EF=eigVec(wi_LST);
        T_EF=eigVec(Ti_LST);
        if strcmp(USE_EOS,'IG')
        p_EF = Rg*(T_LST.*rho_EF + Rho_LST.*T_EF);
        else
        p_EF=eigVec(pi_LST);
        end
        
        %% Phase
        phaseu_EF=angle(u_EF)-angle(p_EF(1)); 
        phasev_EF=angle(v_EF)-angle(p_EF(1));
        %phasew2=angle(w1)-angle(p1(1));
        phaseT_EF=angle(T_EF)-angle(p_EF(1));
        phasep_EF=angle(p_EF)-angle(p_EF(1));
           
% Absolute values and ratios      
        normalize_EF='velocity'; % velocity/pressure
        
        rho_EFAbs       = abs(rho_EF);
        u_EFAbs       = abs(u_EF);
        v_EFAbs       = abs(v_EF);
        w_EFAbs       = abs(w_EF);
        T_EFAbs       = abs(T_EF);
        p_EFAbs       = abs(p_EF);
        
        u_EFMax        = max(u_EFAbs);
        p_EFAbs1       = abs(p_EF(1));
        
        if strcmp(normalize_EF,'velocity')
            den_EF=u_EFMax;
        elseif strcmp(normalize_EF,'pressure')
            den_EF=p_EFAbs1;
        end
        
        rho_EFAbsN      = rho_EFAbs./den_EF;
        u_EFAbsN      = u_EFAbs./den_EF; % p_EFAbs1; %
        v_EFAbsN      = v_EFAbs./den_EF;
        w_EFAbsN      = w_EFAbs./den_EF;
        T_EFAbsN      = T_EFAbs./den_EF; % p_EFAbs1; %
        p_EFAbsN      = p_EFAbs./den_EF;
        
        if strcmp(show_plot,"true")
        Ymax=20;
        fig     = figure(2);
        box
          plot(y_LST, u_EFAbsN,'LineWidth',2)
          hold on
          grid on
          plot(y_LST, T_EFAbsN ,'LineWidth',2)
          plot(y_LST, rho_EFAbsN,'LineWidth',2)
          plot(y_LST, v_EFAbsN,'LineWidth',2)
          plot(y_LST, p_EFAbsN,'LineWidth',2)
          
          ax=gca;
          ax.TickLabelInterpreter='latex';
          legend('$u''$','$T''$','$\rho''$','$v''$','$p''$','interpreter','latex');
          xlabel('$\eta$','interpreter','latex'); ylabel('$|q_{pert}|$','interpreter','latex');
          xlim([0,Ymax])
          set(gca,'XMinorTick','on');
          set(gca,'Fontsize',26,'fontWeight','normal');
          set(findall(gcf,'type','text'),'Fontsize',26,'fontWeight','normal');
        end
          
% Fluctuations
zz=sqrt(-1);
alpha=real(eigV)+zz*imag(eigV);
x=linspace(40.331622296443442,80,40);
for i=1:numel(x)
    for j=1:numel(y_LST)
        u_fluct(i,j)=real(u_EF(j)*exp(zz*(alpha*x(i)-omega*0)));
        v_fluct(i,j)=real(v_EF(j)*exp(zz*(alpha*x(i)-omega*0)));
        T_fluct(i,j)=real(T_EF(j)*exp(zz*(alpha*x(i)-omega*0)));
    end
end

if strcmp(show_plot,"true")
    [X,Y]=meshgrid(x,y_LST);
    figure(3)
    contourf(X,Y,u_fluct','LevelStep',0.2);
    ylim([0 10])
end

LST_EIGS.rho_EF=rho_EF;
LST_EIGS.u_EF=u_EF;
LST_EIGS.v_EF=v_EF;
LST_EIGS.T_EF=T_EF;
LST_EIGS.p_EF=p_EF;

end
    


