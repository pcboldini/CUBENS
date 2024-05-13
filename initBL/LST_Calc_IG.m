function [LST_SPECTRUM] = LST_Calc(BL_PARA,LST_PARA,LST_MESH,FLOW_IG_LST,PRESSURE_IG_LST,THERMO_IG_LST,ENERGY_IG_LST);
    Pr=BL_PARA.Pr_inf;
    Ec=BL_PARA.Ec_inf;    
    M=BL_PARA.M_inf;

    Spatial=LST_PARA.Spatial;
    Full = LST_PARA.Full;
    
    Re0=LST_PARA.Re0;
    beta    = LST_PARA.beta_0;
    if strcmp(Spatial,'true')
    omega   = LST_PARA.omega_0;
    elseif strcmp(Spatial,'false')
    alpha   = LST_PARA.alpha_0;
    end

    dzdy_LST=LST_MESH.dzdy_LST;
    d2zdy2_LST=LST_MESH.d2zdy2_LST;
    N_LST=LST_MESH.N_LST;
    Dist_wall_temp=LST_PARA.Dist_wall_temp;
    Dist_wall_rho=LST_PARA.Dist_wall_rho;
    
    %Parameters
    zz = sqrt(-1);
    Dn = chebDiff(N_LST); %1st order
    Dn2 = Dn * Dn; %2nd order
    II = eye(5);
    Total_Dn = kron(Dn,II);
    Total_DDn= kron(Dn2,II);
        
    AAL1 = zeros([5*N_LST,5*N_LST]); % R
    AAL2 = zeros([5*N_LST,5*N_LST]); % S'+T'
    AAL3 = zeros([5*N_LST,5*N_LST]); % T''

    BBR1 = zeros([5*N_LST,5*N_LST]); % M
    BBR2 = zeros([5*N_LST,5*N_LST]); % N'
    BBR3 = zeros([5*N_LST,5*N_LST]); % P

    for i=1:N_LST
       
        [Lt,Lx,Ly,Lz,Lq,Vxx,Vxy,Vyy,Vxz,Vyz,Vzz]=LST_Calc_Matrix_IG(Re0,Pr,Ec,M,FLOW_IG_LST(i,:),PRESSURE_IG_LST(i,:),THERMO_IG_LST(i,:),ENERGY_IG_LST(i,:));

       if strcmp(Spatial,'true')       
        A1 = -zz*omega*Lt + zz * beta * Lz + Lq - beta * beta * Vzz;
        A2 = Ly + zz * beta * Vyz;
        A3 = Vyy;
        B1 = -zz * Lx + beta * Vxz;
        B2 = -zz * Vxy;
        B3 = Vxx;
        
        AA1 = A1;
        AA2 = A2 * dzdy_LST(i) + A3 * d2zdy2_LST(i);
        AA3 = A3 * dzdy_LST(i) * dzdy_LST(i);
        BB1 = B1;
        BB2 = B2 * dzdy_LST(i);
        BB3 = B3;
        
        icc = (i - 1) * 5;
        
        AAL1(icc+1:icc+5,icc+1:icc+5)=AA1;
        AAL2(icc+1:icc+5,icc+1:icc+5)=AA2;
        AAL3(icc+1:icc+5,icc+1:icc+5)=AA3;
        BBR1(icc+1:icc+5,icc+1:icc+5)=BB1;
        BBR2(icc+1:icc+5,icc+1:icc+5)=BB2;
        BBR3(icc+1:icc+5,icc+1:icc+5)=BB3;     
        
       else
           
        A1 = alpha * Lx + beta * Lz - zz * Lq + zz * alpha * alpha * Vxx + zz * beta * beta * Vzz + zz * alpha * beta * Vxz;
        A2 = -zz * Ly + alpha * Vxy + beta * Vyz;
        A3 = -zz * Vyy; 
        B1 = Lt;
        
        AA1 = A1;
        AA2 = A2 * dzdy_LST(i) + A3 * d2zdy2_LST(i);
        AA3 = A3 * dzdy_LST(i) * dzdy_LST(i);
        BB1 = B1;
        
        icc = (i - 1) * 5;
        
        AAL1(icc+1:icc+5,icc+1:icc+5)=AA1;
        AAL2(icc+1:icc+5,icc+1:icc+5)=AA2;
        AAL3(icc+1:icc+5,icc+1:icc+5)=AA3;
        BBR1(icc+1:icc+5,icc+1:icc+5)=BB1; 
        
       end
    end
    
    AA = AAL1 + AAL2 * Total_Dn + AAL3 * Total_DDn;
    if strcmp(Spatial,'true') 
	BB = BBR1 + BBR2 * Total_Dn;
    else
    BB = BBR1;
    end
    
    if strcmp(Spatial,'true') 
        if strcmp(Full,'false')        
        else
            A=[zeros([5*N_LST,5*N_LST]),eye(5*N_LST);
            AA,-BB]; 
% 
            B=[eye(5*N_LST),zeros([5*N_LST,5*N_LST]);
            zeros([5*N_LST,5*N_LST]),BBR3];   
        end
    end

 %% BCs

         if strcmp(Full,'false') || strcmp(Spatial,'false')
             % Wall
             if strcmp(Dist_wall_rho,'adiabatic') 
                AA(1:5,:)=0;
                AA(1,1)=1.0;
                AA(1,5+1)=-1.0;
             else
                AA(2:5,:)=0;
             end
                
                AA(2,2) = 1.0; % u
                AA(3,3) = 1.0; % v
                AA(4,4) = 1.0; % w
                AA(5,5) = 1.0; % T (isothermal) 
       
                if strcmp(Dist_wall_temp,'adiabatic') 
                AA(5,5+5) = -1.0;  
                end
                if strcmp(Dist_wall_rho,'adiabatic')
                BB(1:5,:)=0;
                else
                BB(2:5,:)=0;
                end
                
              % Farfield    
                AA(5*N_LST-3:5*N_LST,:) = 0;
                AA(5*N_LST-3,5*N_LST-3) = 1.0; % u
                AA(5*N_LST-2,5*N_LST-2) = 1.0; % v
                AA(5*N_LST-1,5*N_LST-1) = 1.0; % w
                AA(5*N_LST,5*N_LST) = 1.0; % T
                
                BB(5*N_LST-3:5*N_LST,:)=0;
                
         elseif strcmp(Full,'true') && strcmp(LST_PARA.Spatial,'true')
                
                % Wall
                
                if strcmp(Dist_wall_rho,'adiabatic')                                
                A(1+5*N_LST:5+5*N_LST,:)=0;
                A(1+5*N_LST,1) = 1.0; % rho
                A(1+5*N_LST,1+5) = -1.0; % rho
                else
                A(2+5*N_LST:5+5*N_LST,:)=0;    
                end
                
                A(2+5*N_LST,2) = 1.0; % u
                A(3+5*N_LST,3) = 1.0; % v
                A(4+5*N_LST,4) = 1.0; % w
                A(5+5*N_LST,5) = 1.0; % T
                
                if strcmp(Dist_wall_temp,'adiabatic')
                A(5+5*N_LST,5+5) = -1.0; % T    
                end               

                if strcmp(Dist_wall_rho,'adiabatic')
                B(1+5*N_LST:5+5*N_LST,:)=0;
                else
                B(2+5*N_LST:5+5*N_LST,:)=0;    
                end
                  
                 % farfield    
                A(5*N_LST-3:5*N_LST,:)=0;
                A(5*N_LST-3,5*N_LST-3) = 1.0;
                A(5*N_LST-2,5*N_LST-2) = 1.0;
                A(5*N_LST-1,5*N_LST-1) = 1.0;
                A(5*N_LST,5*N_LST) = 1.0;
 
                B(5*N_LST-3:5*N_LST,:)=0;
    
         end
% 
 %% EP
         if strcmp(Full,'false') || strcmp(Spatial,'false')
             [Evec,P]=eig(AA,BB); 
         elseif strcmp(Full,'true') && strcmp(Spatial,'true')
             [Evec,P]=eig(A,B);            
         end
 
             Eval=diag(P);
              
LST_SPECTRUM.Evec=Evec;
LST_SPECTRUM.Eval=Eval;
% 
%     figure(1)
%     plot(real(LST_SPECTRUM.Eval),imag(LST_SPECTRUM.Eval),'*');
%     hold on
%  %   for j=1:(5)*LST_MESH.N_LST, text(real(LST_SPECTRUM.Eval(j)), imag(LST_SPECTRUM.Eval(j)), [' ' num2str(j)]), end
%     grid on;
%     if strcmp(LST_PARA.Spatial,'false')
%     axis([-0.05 0.15 -0.10 0.02])
%     else
%     axis([-0.05 0.15 -0.10 0.02])
%     end

end