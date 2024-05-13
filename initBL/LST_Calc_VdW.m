function [LST_SPECTRUM] = LST_Calc_VdW(BL_PARA,LST_PARA,LST_MESH,FLOW_VDW_LST,PRESSURE_VDW_LST,THERMO_VDW_LST,ENERGY_VDW_LST)
    Pr=BL_PARA.Pr_inf;
    Ec=BL_PARA.Ec_inf;    

    Re0=LST_PARA.Re0;
    beta    = LST_PARA.beta_0;
    omega   = LST_PARA.omega_0;
    
    Full = LST_PARA.Full;
    dzdy_LST=LST_MESH.dzdy_LST;
    d2zdy2_LST=LST_MESH.d2zdy2_LST;
    N_LST=LST_MESH.N_LST;
    Dist_wall_temp=LST_PARA.Dist_wall_temp;
    Dist_wall_rho=LST_PARA.Dist_wall_rho;
    Dist_wall_p=LST_PARA.Dist_wall_p;
    
    %Parameters
    zz = sqrt(-1);
    Dn = chebDiff(N_LST); %1st order
    Dn2 = Dn * Dn; %2nd order
    II = eye(6);
    Total_Dn = kron(Dn,II);
    Total_DDn= kron(Dn2,II);
        
    AAL1 = zeros([6*N_LST,6*N_LST]); % R
    AAL2 = zeros([6*N_LST,6*N_LST]); % S'+T'
    AAL3 = zeros([6*N_LST,6*N_LST]); % T''

    BBR1 = zeros([6*N_LST,6*N_LST]); % M
    BBR2 = zeros([6*N_LST,6*N_LST]); % N'
    BBR3 = zeros([6*N_LST,6*N_LST]); % P

    for i=1:N_LST
       
        [Lt,Lx,Ly,Lz,Lq,Vxx,Vxy,Vyy,Vxz,Vyz,Vzz]=LST_Calc_Matrix_VdW(Re0,Pr,Ec,FLOW_VDW_LST(i,:),PRESSURE_VDW_LST(i,:),THERMO_VDW_LST(i,:),ENERGY_VDW_LST(i,:));

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
        
        icc = (i - 1) * 6;
        
        AAL1(icc+1:icc+6,icc+1:icc+6)=AA1;
        AAL2(icc+1:icc+6,icc+1:icc+6)=AA2;
        AAL3(icc+1:icc+6,icc+1:icc+6)=AA3;
        BBR1(icc+1:icc+6,icc+1:icc+6)=BB1;
        BBR2(icc+1:icc+6,icc+1:icc+6)=BB2;
        BBR3(icc+1:icc+6,icc+1:icc+6)=BB3;
    end
    
    AA = AAL1 + AAL2 * Total_Dn + AAL3 * Total_DDn;
	BB = BBR1 + BBR2 * Total_Dn;
   
       A=[zeros([6*N_LST,6*N_LST]),eye(6*N_LST);
            AA,-BB]; 
% 
       B=[eye(6*N_LST),zeros([6*N_LST,6*N_LST]);
            zeros([6*N_LST,6*N_LST]),BBR3];   
   
       II = eye(12);
       Total_Dn = kron(Dn,II);

 %% BCs                           
                % Wall
                
                if strcmp(Dist_wall_rho,'adiabatic')                                
                A(1+6*N_LST:5+6*N_LST,:)=0;
                A(1+6*N_LST,1) = 1.0; % rho
                A(1+6*N_LST,1+6) = -1.0; % rho
                else
                A(2+6*N_LST:5+6*N_LST,:)=0;    
                end
                
                if strcmp(Dist_wall_p,'adiabatic')
                A(6+6*N_LST:6+6*N_LST,:)=0;
                A(6+6*N_LST,6) = 1.0; % p
                A(6+6*N_LST,6+6) = -1.0; % p
                end 
                
                A(2+6*N_LST,2) = 1.0; % u
                A(3+6*N_LST,3) = 1.0; % v
                A(4+6*N_LST,4) = 1.0; % w
                A(5+6*N_LST,5) = 1.0; % T
                
                if strcmp(Dist_wall_temp,'adiabatic')
                A(5+6*N_LST,5+6) = -1.0; % T    
                end               

                if strcmp(Dist_wall_rho,'adiabatic')
                B(1+6*N_LST:5+6*N_LST,:)=0;
                else
                B(2+6*N_LST:5+6*N_LST,:)=0;    
                end
                
                if strcmp(Dist_wall_p,'adiabatic')
                B(6+6*N_LST:6+6*N_LST,:)=0;
                end  
                
                 % farfield    
                A(6*N_LST-4:6*N_LST-1,:)=0;
                A(6*N_LST-4,6*N_LST-4) = 1.0;
                A(6*N_LST-3,6*N_LST-3) = 1.0;
                A(6*N_LST-2,6*N_LST-2) = 1.0;
                A(6*N_LST-1,6*N_LST-1) = 1.0;
 
                B(6*N_LST-4:6*N_LST-1,:)=0;
                
   


 %% EP
             [Evec,P]=eig(A,B);            
 
             Eval=diag(P);
              
LST_SPECTRUM.Evec=Evec;
LST_SPECTRUM.Eval=Eval;

end