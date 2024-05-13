CASE  = "BoundaryLayer"
! ------------ use equation of state ------------------ 
USE_EOS  = "IG"
! ------------ set non-dimensional free-stream values for computation ------------------ 
delta99 =4.9141200166
Pra =0.7500000000
Ec =4.0000000000e-03
Ma =1.0000000000e-01
Tinf =300.000000
Pref =71.428571
ig_gam =1.400000
eos_Rgas =71.428571
eos_dof =9.000000
! ------------ set wall BC for computation ------------------ 
wall_bc  = "adiab"
Twall_bottom =-1.000000
! ------------ set viscosity and conductivity ------------------ 
USE_VISC  = "Sutherland"
Stref =273.0000000000
Muinf =1.8417735821e-05
Muref =1.7160000000e-05
Smuref =111.0000000000
Kinf =2.5866400541e-02
Kref =2.4100000000e-02
Skref =194.0000000000
