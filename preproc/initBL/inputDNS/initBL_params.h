!------------ use equation of state ------------------ 
USE_EOS  = "PR" 
!------------ set non-dimensional free-stream values for computation ------------------ 
delta99   = 5.5648011057 
Pra   = 1.0000000000 
Ec   = 1.0000000000e-02 
Ma   = 1.7834290694e-01 
eos_dof   = 9.000000 
eos_ac   = 0.224000 
!------------ set dimensional free-stream values for computation ------------------ 
Tcrit   = 304.128200 
Pcrit   = 7377300.000000 
Vcrit   = 2.138580e-03 
eos_Rgas   = 188.924058 
! ------------ set wall BC for computation ------------------ 
wall_bc  = "isoth" 
Twall_bot   = 0.900000 
!------------ set reference free-stream values for computation ------------------ 
Tref   = 0.9000000000 
Pref   = 1.1000000000 
Rhoref   = 2.3366345038 
Cpref   = 13.8693915475 
SOSref   = 3.5511933911 
Rref   = 3.2133676093 
! ----------- set viscosity and conductivity ---------------- 
USE_VISC  = "JST" 
Muref   = 3.2420139420e-03 
Kref   = 9.0490527358e-02 
