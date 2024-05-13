INIT  = "initBoundaryLayer"
! ------------ use equation of state ------------------ 
USE_EOS  = "PR"
! ------------ set non-dimensional free-stream values for computation ------------------ 
delta99 =4.9022850386
Pra =2.3945206885
Ec =2.4313796248e-03
Ma =1.0000000000e-01
eos_dof =9.000000
eos_ac =0.224000
! ------------ set dimensional free-stream values for computation ------------------ 
Tcrit =304.128200
Pcrit =7377300.000000
Vcrit =2.138580e-03
eos_Rgas =188.924058
! ------------ set wall BC for computation ------------------ 
wall_bc  = "adiab"
Twall =-1.000000
! ------------ set reference free-stream values for computation ------------------ 
Tref =0.9207000000
Pref =1.0844000000
Rhoref =2.2169228941
Cpref =14.9734319261
SOSref =3.2819075284
Rref =3.213368
! ------------ set viscosity and conductivity ------------------ 
USE_VISC  = "Chung"
Muref =9.6520026435e-05
Kref =1.1402693638e-01
