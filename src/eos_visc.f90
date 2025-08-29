! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
! viscosity and conductivity module

module mod_eos_visc
  use decomp_2d, only: mytype
  implicit none 

! Variables container for Constant Law
  type t_VISC_Constant
    real(mytype) :: mu_0, ka_0
  end type
  type(t_VISC_Constant) :: t_const
  !$acc declare create(t_const)
  !$acc declare create(t_const%mu_0,t_const%ka_0)

! Variables container for Power Law
  type t_VISC_PowerLaw
    real(mytype) :: mu_0, ka_0
  end type
  type(t_VISC_PowerLaw) :: t_pow
  !$acc declare create(t_pow)
  !$acc declare create(t_pow%mu_0,t_pow%ka_0)

! Variables container for Sutherland
  type t_VISC_Suth
    real(mytype) :: mu_0, ka_0, Sconst
  end type
  type(t_VISC_Suth) :: t_suth
  !$acc declare create(t_suth)
  !$acc declare create(t_suth%mu_0,t_suth%ka_0,t_suth%Sconst)

! Variables container for JST (Jossi, Stiel, Thodos)
  type t_VISC_JST
    real(mytype) :: mu_0, ka_0
    real(mytype) :: Kref,Muref,eos_Ru
    real(mytype) :: coeffrho_1,coeffrho_2,coeffrho_3,coeffrho_4,coeffrho_5,ka_factor,ZcPow5
  end type
  type(t_VISC_JST) :: t_jst
  !$acc declare create(t_jst)
  !$acc declare create(t_jst%mu_0,t_jst%ka_0,t_jst%Kref,t_jst%Muref,t_jst%eos_Ru, &
  !$acc                t_jst%coeffrho_1,t_jst%coeffrho_2,t_jst%coeffrho_3, &
  !$acc                t_jst%coeffrho_4,t_jst%coeffrho_5, &
  !$acc                t_jst%ka_factor,t_jst%ZcPow5)

! Variables container for Chung
  type t_VISC_Chung
    real(mytype) :: mu_0, ka_0
    real(mytype) :: mol,Tcrit,chung_Vcrit,ek,Zc,rp_Zc,Kref,Muref
    real(mytype) :: chung_A,chung_B,chung_C,chung_D,chung_E,chung_F,chung_G,chung_H,chung_S,chung_W
    real(mytype) :: chung_Fc,chung_alpha,chung_beta
    real(mytype), dimension(10) :: chung_a0,chung_a1,chung_Arho
    real(mytype), dimension(7)  :: chung_b0,chung_b1,chung_Brho
  end type
  type(t_VISC_Chung) :: t_chu
  !$acc declare create(t_chu)
  !$acc declare create(t_chu%mu_0,t_chu%ka_0,t_chu%mol,t_chu%Tcrit,t_chu%chung_Vcrit, &
  !$acc                t_chu%ek,t_chu%Zc,t_chu%rp_Zc,t_chu%Kref,t_chu%Muref, &
  !$acc                t_chu%chung_A,t_chu%chung_B,t_chu%chung_C,t_chu%chung_D, &
  !$acc                t_chu%chung_E,t_chu%chung_F,t_chu%chung_G,t_chu%chung_H, &
  !$acc                t_chu%chung_S,t_chu%chung_W, &
  !$acc                t_chu%chung_Fc,t_chu%chung_alpha,t_chu%chung_beta, &
  !$acc                t_chu%chung_a0,t_chu%chung_a1,t_chu%chung_Arho, &
  !$acc                t_chu%chung_b0,t_chu%chung_b1,t_chu%chung_Brho)
  contains 


! Include Constant law 
#ifdef Constant
! Initialize

  subroutine init_VISCModel()
    use mod_param
    implicit none
    integer                            :: ierr 
#if defined(BL) || defined(CHA)
    if (nrank==0) then
      if (USE_VISC=='Constant') then
        write(stdout,* ) 'Correct initialisation of USE_VISC'
        write(stdout,* )
      else
        write(stdout,* ) 'Mismatch USE_VISC between initBL and Makefile, check both again!'
        write(stdout,* )
        call decomp_2d_finalize
        call mpi_finalize(ierr)
        stop
      endif
    endif
#endif
    ! Declare variables
    t_const%mu_0 = 1.0_mytype/Re
    t_const%ka_0 = 1.0_mytype/(Re*Pra*Ec)
    !$acc update device(t_const)
    !$acc update device(t_const%mu_0,t_const%ka_0)
  end subroutine
! Calculate viscosity and conductivity according to Constant law
  subroutine calcVisc(tem, rho, mu, ka)
    !$acc routine seq
    implicit none
    real(mytype), intent(IN)  :: tem, rho
    real(mytype), intent(OUT) :: mu, ka
    mu = t_const%mu_0*tem
    ka = t_const%ka_0*tem
  end subroutine
#endif

! Include Power Law  
#ifdef PowerLaw
! Initialize
  subroutine init_VISCModel()
    use mod_param
    implicit none
    integer                            :: ierr 
#if defined(BL) || defined(CHA)
    if (nrank==0) then
      if (USE_VISC=='PowerLaw') then
        write(stdout,* ) 'Correct initialisation of USE_VISC'
        write(stdout,* )
      else
        write(stdout,* ) 'Mismatch USE_VISC between initBL and Makefile, check both again!'
        write(stdout,* )
        call decomp_2d_finalize
        call mpi_finalize(ierr)
        stop
      endif
    endif
#endif
    ! Declare variables
    t_pow%mu_0 = 1.0_mytype/Re
    t_pow%ka_0 = 1.0_mytype/(Re*Pra*Ec)
    !$acc update device(t_pow)
    !$acc update device(t_pow%mu_0,t_pow%ka_0)
  end subroutine
! Calculate viscosity and conductivity according to Power Law
  subroutine calcVisc(tem, rho, mu, ka)
    !$acc routine seq
    implicit none
    real(mytype), intent(IN)  :: tem, rho
    real(mytype), intent(OUT) :: mu, ka
    real(mytype) :: fact, exp_powerl
    exp_powerl = 3.0_mytype/4.0_mytype 
    fact = tem**exp_powerl
    mu = t_pow%mu_0*fact
    ka = t_pow%ka_0*fact
  end subroutine
#endif

! Include Sutherland
#ifdef Sutherland
! Initialize
  subroutine init_VISCModel()
    use mod_param
    implicit none
    integer                            :: ierr 
#if defined(BL) || defined(CHA)
    if (nrank==0) then
      if (USE_VISC=='Sutherland') then
        write(stdout,* ) 'Correct initialisation of USE_VISC'
        write(stdout,* )
      else
        write(stdout,* ) 'Mismatch USE_VISC between initBL and Makefile, check both again!'
        write(stdout,* )
        call decomp_2d_finalize
        call mpi_finalize(ierr)
        stop
      endif
    endif
#endif
    ! Declare variables
    t_suth%Sconst = Smuref/Tinf
    t_suth%mu_0 = 1.0_mytype/Re
    t_suth%ka_0 = 1.0_mytype/(Re*Pra*Ec)
    !$acc update device(t_suth)
    !$acc update device(t_suth%Sconst,t_suth%mu_0,t_suth%ka_0)
  end subroutine

! Calculate viscosity and conductivity according to Sutherland
  subroutine calcVisc(tem, rho, mu, ka)
    !$acc routine seq  
    implicit none
    real(mytype), intent(IN)  :: tem, rho
    real(mytype), intent(OUT) :: mu, ka
    real(mytype) :: suth
    suth = ((1.0_mytype + t_suth%Sconst)/(tem + t_suth%Sconst))*tem**(3.0_mytype/2.0_mytype)
    mu = t_suth%mu_0*suth
    ka = t_suth%ka_0*suth
  end subroutine
#endif


! Include JST
#ifdef JST
! Initialize

  subroutine init_VISCModel()
    use mod_param
    implicit none
    integer                            :: ierr                   
#if defined(BL) || defined(CHA)
    if (nrank==0) then
      if (USE_VISC=='JST') then
        write(stdout,* ) 'Correct initialisation of USE_VISC'
        write(stdout,* )
      else
        write(stdout,* ) 'Mismatch USE_VISC between initBL and Makefile, check both again!'
        write(stdout,* )
        call decomp_2d_finalize
        call mpi_finalize(ierr)
        stop
      endif
    endif
#endif
    ! Declare variables
    t_jst%mu_0 = 1.0_mytype/Re
    t_jst%ka_0 = 1.0_mytype/(Re*Pra*Ec)
    t_jst%Kref = Kref
    t_jst%Muref = Muref
    t_jst%eos_Ru = eos_Ru
    t_jst%coeffrho_1 = 0.10230_mytype 
    t_jst%coeffrho_2 = 0.023364_mytype 
    t_jst%coeffrho_3 = 0.058533_mytype 
    t_jst%coeffrho_4 = 0.040758_mytype 
    t_jst%coeffrho_5 = 0.0093324_mytype 
    t_jst%ka_factor  = 0.307*(eos_dof/2.0_mytype) + 0.539_mytype
#if defined(VdW)
    t_jst%ZcPow5=vdw_Zc**5
#elif defined(RK)
    t_jst%ZcPow5=rk_Zc**5
#elif defined(VdW)
    t_jst%ZcPow5=pr_Zc**5
#endif
    !$acc update device(t_jst)
    !$acc update device(t_jst%mu_0,t_jst%ka_0,t_jst%Kref,t_jst%Muref,t_jst%eos_Ru, &
    !$acc                t_jst%coeffrho_1,t_jst%coeffrho_2,t_jst%coeffrho_3, &
    !$acc                t_jst%coeffrho_4,t_jst%coeffrho_5, &
    !$acc                t_jst%ka_factor,t_jst%ZcPow5)
  end subroutine

! Calculate viscosity and conductivity according to JST
  subroutine calcVisc(tem_r, rho_r, mu, ka)
    !$acc routine seq
    implicit none
    real(mytype), intent(IN)  :: tem_r, rho_r
    real(mytype), intent(OUT) :: mu, ka
    real(mytype) :: mu1_r, mu_diff, mu_2, kappa1_r, f_rho_ka, ka_diff, ka_2, ZcPow5
    ! Viscosity
    if (tem_r .le. 1.50_mytype) then
        mu1_r = 34E-5_mytype*tem_r**0.94_mytype
    else
        mu1_r = 17.78E-5_mytype*(4.58_mytype*tem_r-1.67_mytype)**(5.0_mytype/8.0_mytype)
    endif 
    mu_diff =  (t_jst%coeffrho_1           &
                + t_jst%coeffrho_2*rho_r**1 &
                + t_jst%coeffrho_3*rho_r**2 &
                - t_jst%coeffrho_4*rho_r**3 & 
                + t_jst%coeffrho_5*rho_r**4)**4-1E-4_mytype
    mu_2     = mu_diff+mu1_r
    mu = mu_2/t_jst%Muref*t_jst%mu_0  
    ! Conductivity
    ZcPow5 = t_jst%ZcPow5
    kappa1_r = 15.0_mytype/4.0_mytype*t_jst%eos_Ru*mu1_r*t_jst%ka_factor
    if (rho_r < 0.50_mytype) then
      f_rho_ka = 14.0_mytype*(exp(0.535_mytype*rho_r)-1.0_mytype)
      ka_diff = f_rho_ka*1E-8_mytype/ZcPow5
      ka_2 = ka_diff*4.1868*1E2_mytype + kappa1_r
    else if ((rho_r .GE. 0.50_mytype) .and. (rho_r<2.0_mytype)) then
      f_rho_ka = 13.1_mytype*(exp(0.67_mytype*rho_r)-1.069_mytype)
      ka_diff = f_rho_ka*1E-8_mytype/ZcPow5
      ka_2 = ka_diff*4.1868*1E2_mytype + kappa1_r
    else
      f_rho_ka = 2.976_mytype*(exp(1.155_mytype*rho_r)+2.016_mytype)
      ka_diff = f_rho_ka*1E-8_mytype/ZcPow5
      ka_2 = ka_diff*4.1868*1E2_mytype + kappa1_r
    endif
    ka = ka_2/t_jst%Kref*t_jst%ka_0 
  end subroutine
#endif


! Include Chung
#ifdef Chung
! Initialize
  subroutine init_VISCModel()
    use mod_param
    implicit none
    integer                            :: ierr                        
#if defined(BL) || defined(CHA)
    if (nrank==0) then
      if (USE_VISC=='Chung') then
        write(stdout,* ) 'Correct initialisation of USE_VISC'
        write(stdout,* )
      else
        write(stdout,* ) 'Mismatch USE_VISC between initBL and Makefile, check both again!'
        write(stdout,* )
        call decomp_2d_finalize
        call mpi_finalize(ierr)
        stop
      endif
    endif
#endif
    ! Declare variables
    t_chu%mu_0 = 1.0_mytype/Re
    t_chu%ka_0 = 1.0_mytype/(Re*Pra*Ec)
    t_chu%mol     = eos_Ru/eos_Rgas
    t_chu%Tcrit   = Tcrit
    t_chu%chung_Vcrit = Vcrit*t_chu%mol*1E6_mytype
    t_chu%ek          = Tcrit/1.2593_mytype
    t_chu%rp_Zc = rp_Zc
#if defined(VdW)
      t_chu%Zc = vdw_Zc
#elif defined(RK)
      t_chu%Zc = rk_Zc      
#elif defined(PR)
      t_chu%Zc = pr_Zc
#endif
    t_chu%Kref  = Kref
    t_chu%Muref = Muref
    t_chu%chung_alpha = eos_dof/2 - 3.0_mytype/2.0_mytype
    t_chu%chung_beta  = 0.7862_mytype - 0.7109_mytype*eos_ac + 1.3168_mytype*eos_ac**2
    t_chu%chung_A =  1.16145_mytype
    t_chu%chung_B =  0.14874_mytype
    t_chu%chung_C =  0.52487_mytype
    t_chu%chung_D =  0.77320_mytype
    t_chu%chung_E =  2.16178_mytype
    t_chu%chung_F =  2.43787_mytype
    t_chu%chung_G =   -6.435_mytype*1E-4
    t_chu%chung_H =  7.27371_mytype
    t_chu%chung_S =  18.0323_mytype
    t_chu%chung_W = -0.76830_mytype
    t_chu%chung_Fc = 1 - 0.2756_mytype*eos_ac
    t_chu%chung_a0 = (/6.32402_mytype,0.0012102_mytype, 5.28346_mytype, 6.62263_mytype, 19.7454_mytype, &
                      -1.89992_mytype,  24.2745_mytype, 0.79716_mytype,-0.23816_mytype,0.068629_mytype/)
    t_chu%chung_a1 = (/50.4119_mytype,-0.0011536_mytype,254.209_mytype, 38.0957_mytype,7.63034_mytype, &
                     -12.5367_mytype,  3.44945_mytype, 1.11764_mytype,0.067695_mytype,0.34793_mytype/)
    t_chu%chung_Arho = t_chu%chung_a0 + t_chu%chung_a1*eos_ac
    t_chu%chung_b0 = (/2.41657_mytype,-0.50924_mytype,6.61069_mytype,14.54250_mytype,0.79274_mytype, &
                      -5.86340_mytype,81.171_mytype/)
    t_chu%chung_b1 = (/0.74824_mytype,-1.50936_mytype,5.62073_mytype,-8.91387_mytype,0.82019_mytype, &
                      12.80050_mytype,114.158_mytype/)
    t_chu%chung_Brho = t_chu%chung_b0 + t_chu%chung_b1*eos_ac
    !$acc update device(t_chu)
    !$acc update device(t_chu%mu_0,t_chu%ka_0,t_chu%mol,t_chu%Tcrit,t_chu%chung_Vcrit, &
    !$acc                t_chu%ek,t_chu%Zc,t_chu%rp_Zc,t_chu%Kref,t_chu%Muref, &
    !$acc                t_chu%chung_A,t_chu%chung_B,t_chu%chung_C,t_chu%chung_D, &
    !$acc                t_chu%chung_E,t_chu%chung_F,t_chu%chung_G,t_chu%chung_H, &
    !$acc                t_chu%chung_S,t_chu%chung_W, &
    !$acc                t_chu%chung_Fc,t_chu%chung_alpha,t_chu%chung_beta, &
    !$acc                t_chu%chung_a0,t_chu%chung_a1,t_chu%chung_Arho, &
    !$acc                t_chu%chung_b0,t_chu%chung_b1,t_chu%chung_Brho)
 end subroutine

! Calculate viscosity and conductivity according to Chung
  subroutine calcVisc(tem_r, rho_r, mu, ka)
    !$acc routine seq  
    implicit none
    real(mytype), intent(IN)  :: tem_r, rho_r
    real(mytype), intent(OUT) :: mu, ka
    real(mytype) :: tem_dimless,chung_psi_mu,mu_1,chung_Y,chung_G1,chung_G2,mu_k,mu_p
    real(mytype) :: chung_Z,chung_psi_ka,kappa_1,chung_H2,kappa_k,kappa_p,tem_dim
    ! Viscosity
    tem_dim = tem_r*t_chu%Tcrit
    tem_dimless  = tem_dim/t_chu%ek
    chung_psi_mu = t_chu%chung_A/tem_dimless**(t_chu%chung_B) &
                  +t_chu%chung_C/exp(t_chu%chung_D*tem_dimless) &
                  +t_chu%chung_E/exp(t_chu%chung_F*tem_dimless) &
                  +t_chu%chung_G*tem_dimless**(t_chu%chung_B)*sin(t_chu%chung_S*tem_dimless**t_chu%chung_W-t_chu%chung_H)
    mu_1 = 4.0785_mytype*1e-6_mytype*t_chu%chung_Fc*(t_chu%mol*1E3_mytype*tem_dim)**0.5_mytype &
         /(t_chu%chung_Vcrit**(2.0_mytype/3.0_mytype)*chung_psi_mu)
    chung_Y  = rho_r*t_chu%rp_Zc/t_chu%Zc/6.0_mytype
    chung_G1 = (1.0_mytype-0.5_mytype*chung_Y)/(1.0_mytype-chung_Y)**3
    chung_G2 = ( (t_chu%chung_Arho(1)*(1-exp(-t_chu%chung_Arho(4)*chung_Y))/chung_Y &
             + t_chu%chung_Arho(2)*chung_G1*exp(t_chu%chung_Arho(5)*chung_Y) + t_chu%chung_Arho(3)*chung_G1 ) &
             /( t_chu%chung_Arho(1)*t_chu%chung_Arho(4)+t_chu%chung_Arho(2)+t_chu%chung_Arho(3) ) )
    mu_k = mu_1*(1/chung_G2 + t_chu%chung_Arho(6)*chung_Y) 
    mu_p = ( 3.6344_mytype*1E-6_mytype*(t_chu%mol*1E3_mytype*t_chu%Tcrit)**0.5_mytype/(t_chu%chung_Vcrit**(2.0_mytype/3.0_mytype)) &
         *(t_chu%chung_Arho(7)*chung_Y**2*chung_G2 &
         *exp(t_chu%chung_Arho(8)+t_chu%chung_Arho(9)*tem_dimless**(-1)+t_chu%chung_Arho(10)*tem_dimless**(-2))) )
    mu = (mu_k + mu_p)/t_chu%Muref*t_chu%mu_0
    ! Conductivity
    chung_Z = 2.0_mytype+10.5_mytype*tem_r**2
    chung_psi_ka = ( 1.0_mytype+t_chu%chung_alpha*(0.215_mytype+0.28288_mytype*t_chu%chung_alpha &
                 - 1.061_mytype*t_chu%chung_beta+0.26665_mytype*chung_Z ) &
                 /(0.6366_mytype+t_chu%chung_beta*chung_Z+1.061_mytype*t_chu%chung_alpha*t_chu%chung_beta) )
    kappa_1 = 31.2_mytype*(mu_1/t_chu%mol)*chung_psi_ka
    chung_H2 = ( t_chu%chung_Brho(1)*(1.0_mytype-exp(-t_chu%chung_Brho(4)*chung_Y))/chung_Y + t_chu%chung_Brho(2)*chung_G1 &
            *exp(t_chu%chung_Brho(5)*chung_Y) &
            + t_chu%chung_Brho(3)*chung_G1)/(t_chu%chung_Brho(1)*t_chu%chung_Brho(4)+t_chu%chung_Brho(2)+t_chu%chung_Brho(3))
    kappa_k = kappa_1*(1.0_mytype/chung_H2+t_chu%chung_Brho(6)*chung_Y)
    kappa_p = ( 3.586_mytype*1E-3_mytype*(t_chu%Tcrit/t_chu%mol)**0.5_mytype/(t_chu%chung_Vcrit**(2.0_mytype/3.0_mytype)) &
            *( t_chu%chung_Brho(7)*chung_Y**2*chung_H2*tem_r**0.5_mytype ) )
    ka = (kappa_k+kappa_p)/t_chu%Kref*t_chu%ka_0
  end subroutine
#endif

  
end module mod_eos_visc
