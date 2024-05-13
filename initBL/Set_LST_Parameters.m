function [LST_PARA]=Set_LST_Parameters

global threed dimomega viscous_LST local_LST inviscid_LST

LST_PARA.Spatial   = 'true';
LST_PARA.Full   = 'false';    % Spatial/Temporal problem, full matrix, if false alfa^2/omega^2 dropped
LST_PARA.NN        = 0;       % numbers of eigenvalues to be solved 
LST_PARA.Re0       = 250; %2000; 

threed     = 'false'; % true/false
dimomega   = 'false'; % true/false
viscous_LST    = 'true'; % true/false
local_LST      = 'false'; % true/false
inviscid_LST   = 'spectrum'; % shooting/spectrum

LST_PARA.sigma= 0.245577+0.0134408i; %0.220195-0.00309131i; %0.056853-0.00031235i %0.22719+0.00283982i; %  0.220195-0.00309131i %0.22719+0.00283982i; %0.240404-0.00345523i;% 0.471596-0.0146847i %;%;%0.22719+0.00283982i; %0.226338-0.0102011i;%0.22719+0.00283982i 0.2171660-0.0059525i 0.220157-0.0034701i

if strcmp(LST_PARA.Spatial,'true')
    if strcmp(dimomega,'true')
        LST_PARA.omega_0     = 0.05;   %0.2032 0.08
        LST_PARA.F_0=LST_PARA.omega_0/LST_PARA.Re0;
    elseif strcmp(dimomega,'false')
        LST_PARA.F_0=110*1e-6;
        LST_PARA.omega_0     = LST_PARA.F_0*LST_PARA.Re0;  
    end
else
LST_PARA.alpha_0     = 0;    
end
if strcmp(threed,'false')
LST_PARA.beta_0      = 0;
else
LST_PARA.beta_0      = 0.1;   
end

% Disturbance boundary conditions

LST_PARA.Dist_wall_temp = 'constant'; %constant/adiabatic
LST_PARA.Dist_wall_rho = 'constant'; %constant/adiabatic 
LST_PARA.Dist_wall_p = 'constant'; %constant/adiabatic

end

