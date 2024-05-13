function [FLOW_VDW_LST,PRESSURE_VDW_LST,THERMO_VDW_LST,ENERGY_VDW_LST] = BaseFlowVariables_interpol_VdW(Neta,Ny_ext,y,ym,BL_PARA,FLOW_VDW,PRESSURE_VDW,THERMO_VDW,ENERGY_VDW)
    
    FLOW_VDW_inf = [1,0,0,1,0,0,1,0,0,BL_PARA.Pr_inf,BL_PARA.a_inf,BL_PARA.V0_ns_inf];
    PRESSURE_VDW_inf=PRESSURE_VDW(Neta,:);
    THERMO_VDW_inf=THERMO_VDW(Neta,:);
    ENERGY_VDW_inf=ENERGY_VDW(Neta,:);
    
    %==========================================================================
    % Interpolation from eta to y
    y_ext = linspace(y(Neta),500,Ny_ext+1)'; 
    y_all = [y;y_ext(2:Ny_ext+1)]; 
    
    
    FLOW_VDW_LST = zeros(size(ym,1),size(FLOW_VDW,2));
    for j = 1:size(FLOW_VDW_LST,2)
        FLOW_VDW_LST(:,j) = interp1(y_all,[FLOW_VDW(:,j);FLOW_VDW_inf(j)*ones(Ny_ext,1)],ym,'spline');
    end
    
    PRESSURE_VDW_LST = zeros(size(ym,1),size(PRESSURE_VDW,2));
    for j = 1:size(PRESSURE_VDW_LST,2)
        PRESSURE_VDW_LST(:,j) = interp1(y_all,[PRESSURE_VDW(:,j);PRESSURE_VDW_inf(j)*ones(Ny_ext,1)],ym,'spline');
    end
    
    THERMO_VDW_LST = zeros(size(ym,1),size(THERMO_VDW,2));
    for j = 1:size(THERMO_VDW_LST,2)
        THERMO_VDW_LST(:,j) = interp1(y_all,[THERMO_VDW(:,j);THERMO_VDW_inf(j)*ones(Ny_ext,1)],ym,'spline');
    end
    
    ENERGY_VDW_LST = zeros(size(ym,1),size(ENERGY_VDW,2));
    for j = 1:size(ENERGY_VDW_LST,2)
        ENERGY_VDW_LST(:,j) = interp1(y_all,[ENERGY_VDW(:,j);ENERGY_VDW_inf(j)*ones(Ny_ext,1)],ym,'spline');
    end
    
    %==========================================================================

end