!------------ use equation of state ------------------ 
USE_EOS  = "PR" 
!------------ set non-dimensional free-stream values for computation ------------------ 
delta99   = 4.3798443993 
Pra   = 2.3947469379 
Ec   = 1.0000000000e-02 
Ma   = 2.0283025635e-01 
eos_dof   = 9.000000 
eos_ac   = 0.224000 
!------------ set dimensional free-stream values for computation ------------------ 
Tcrit   = 304.128200 
Pcrit   = 7377300.000000 
Vcrit   = 2.138580e-03 
eos_Rgas   = 188.924058 
! ------------ set wall BC for computation ------------------ 
wall_bc  = "isoth" 
Twall_bot   = 1.125000 
!------------ set reference free-stream values for computation ------------------ 
Tref   = 0.9207000000 
Pref   = 1.0840000000 
Rhoref   = 2.2168294138 
Cpref   = 14.9751171285 
SOSref   = 3.2816454038 
Rref   = 3.2133676093 
! ----------- set viscosity and conductivity ---------------- 
USE_VISC  = "Chung" 
Muref   = 9.6510017117e-05 
Kref   = 1.1401717046e-01 
