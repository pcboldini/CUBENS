function [FLOW_IG_LST,PRESSURE_IG_LST,THERMO_IG_LST,ENERGY_IG_LST] = BaseFlowVariables_interpol_IG(Neta,Ny_ext,y,ym,BL_PARA,FLOW_IG,PRESSURE_IG,THERMO_IG,ENERGY_IG)
    
    FLOW_IG_inf = [1,0,0,1,0,0,1,0,0,BL_PARA.Pr_inf,BL_PARA.a_inf,BL_PARA.V0_ns_inf]; % ,
    PRESSURE_IG_inf=PRESSURE_IG(Neta,:);
    THERMO_IG_inf=THERMO_IG(Neta,:);
    ENERGY_IG_inf=ENERGY_IG(Neta,:);
    
    %==========================================================================
    % Interpolation from eta to y
    y_ext = linspace(y(Neta),500,Ny_ext+1)'; 
    y_all = [y;y_ext(2:Ny_ext+1)]; 
    
    
    FLOW_IG_LST = zeros(size(ym,1),size(FLOW_IG,2));
    for j = 1:size(FLOW_IG_LST,2)
        FLOW_IG_LST(:,j) = interp1(y_all,[FLOW_IG(:,j);FLOW_IG_inf(j)*ones(Ny_ext,1)],ym,'spline');
    end
    
    PRESSURE_IG_LST = zeros(size(ym,1),size(PRESSURE_IG,2));
    for j = 1:size(PRESSURE_IG_LST,2)
        PRESSURE_IG_LST(:,j) = interp1(y_all,[PRESSURE_IG(:,j);PRESSURE_IG_inf(j)*ones(Ny_ext,1)],ym,'spline');
    end
    
    THERMO_IG_LST = zeros(size(ym,1),size(THERMO_IG,2));
    for j = 1:size(THERMO_IG_LST,2)
        THERMO_IG_LST(:,j) = interp1(y_all,[THERMO_IG(:,j);THERMO_IG_inf(j)*ones(Ny_ext,1)],ym,'spline');
    end
    
    ENERGY_IG_LST = zeros(size(ym,1),size(ENERGY_IG,2));
    for j = 1:size(ENERGY_IG_LST,2)
        ENERGY_IG_LST(:,j) = interp1(y_all,[ENERGY_IG(:,j);ENERGY_IG_inf(j)*ones(Ny_ext,1)],ym,'spline');
    end
    
    %==========================================================================

end