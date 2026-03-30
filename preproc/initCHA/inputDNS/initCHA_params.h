!------------ use equation of state ------------------ 
USE_EOS  = "PR" 
!------------ set non-dimensional reference values for computation ------------------ 
Pra   = 2.0787234544 
Ec   = 4.0000000000e-02 
Ma   = 1.7823787608e-01 
eos_dof   = 5.000000 
eos_ac   = 0.037200 
Ri_wall   = 0.000000 
Ri_unit   = 0.000000 
!------------ set dimensional values for computation ------------------ 
Tcrit   = 126.192000 
Pcrit   = 3395800.000000 
Vcrit   = 3.028575e-03 
eos_Rgas   = 296.803896 
! ------------ set wall BC for computation ------------------ 
Twall_bot   = 1.000000 
Twall_top   = 1.000000 
!------------ set non-dimensional reference values for computation ------------------ 
Tref   = 0.7900000000 
Pref   = 1.1400000000 
Rhoref   = 2.6191707535 
Cpref   = 7.4057624237 
SOSref   = 4.9403943702 
Rref   = 3.3133428316 
! ----------- set viscosity and conductivity ---------------- 
USE_VISC  = "Chung" 
Muref   = 9.0916915030e-05 
Kref   = 9.6136287658e-02 
