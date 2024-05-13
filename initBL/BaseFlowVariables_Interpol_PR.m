function [FLOW_PR_LST,PRESSURE_PR_LST,THERMO_PR_LST,ENERGY_PR_LST] = BaseFlowVariables_Interpol_PR(Neta,Ny_ext,y,ym,BL_PARA,FLOW_PR,PRESSURE_PR,THERMO_PR,ENERGY_PR)
    
    FLOW_PR_inf = [1,0,0,1,0,0,1,0,0,BL_PARA.Pr_inf,BL_PARA.a_inf,BL_PARA.V0_ns_inf];
    PRESSURE_PR_inf=PRESSURE_PR(Neta,:);
    THERMO_PR_inf=THERMO_PR(Neta,:);
    ENERGY_PR_inf=ENERGY_PR(Neta,:);
    
    %==========================================================================
    % Interpolation from eta to y
    y_ext = linspace(y(Neta),500,Ny_ext+1)'; 
    y_all = [y;y_ext(2:Ny_ext+1)]; 
    
    
    FLOW_PR_LST = zeros(size(ym,1),size(FLOW_PR,2));
    for j = 1:size(FLOW_PR_LST,2)
        FLOW_PR_LST(:,j) = interp1(y_all,[FLOW_PR(:,j);FLOW_PR_inf(j)*ones(Ny_ext,1)],ym,'spline');
    end
    
    PRESSURE_PR_LST = zeros(size(ym,1),size(PRESSURE_PR,2));
    for j = 1:size(PRESSURE_PR_LST,2)
        PRESSURE_PR_LST(:,j) = interp1(y_all,[PRESSURE_PR(:,j);PRESSURE_PR_inf(j)*ones(Ny_ext,1)],ym,'spline');
    end
    
    THERMO_PR_LST = zeros(size(ym,1),size(THERMO_PR,2));
    for j = 1:size(THERMO_PR_LST,2)
        THERMO_PR_LST(:,j) = interp1(y_all,[THERMO_PR(:,j);THERMO_PR_inf(j)*ones(Ny_ext,1)],ym,'spline');
    end
    
    ENERGY_PR_LST = zeros(size(ym,1),size(ENERGY_PR,2));
    for j = 1:size(ENERGY_PR_LST,2)
        ENERGY_PR_LST(:,j) = interp1(y_all,[ENERGY_PR(:,j);ENERGY_PR_inf(j)*ones(Ny_ext,1)],ym,'spline');
    end
    
    %==========================================================================

end