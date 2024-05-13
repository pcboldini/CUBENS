! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_eos_var
  use decomp_2d, only: mytype
  implicit none 
! Variables container for ideal gas EoS   
  type t_EOS_IG
    real(mytype) :: cp, cv, cvInv, Rgas, gam, prefac_r
  end type 
  type(t_EOS_IG) :: t_ig
  !$acc declare create(t_ig)
  !$acc declare create(t_ig%cp,t_ig%cv,t_ig%cvInv,t_ig%Rgas,t_ig%gam,t_ig%prefac_r)
! Variables container for Van der Waals EoS   
  type t_EOS_VDW
    real(mytype) :: Efac_r, cvOverR, prefac_r, Cpfac_r, Zc, Rref, a
  end type 
  type(t_EOS_VDW) :: t_vdw
  !$acc declare create(t_vdw)
  !$acc declare create(t_vdw%Efac_r,t_vdw%cvOverR, & 
  !$acc                t_vdw%prefac_r,t_vdw%Cpfac_r,t_vdw%Zc,t_vdw%Rref,t_vdw%a)
! Variables container for Peng-Robinson EoS  
  type t_EOS_PR
    real(mytype) :: Rref, Efac_r, cvOverR, prefac_r, Cpfac_r, K, Zc, Zc1, Zc2, a, b
  end type 
  type(t_EOS_PR) :: t_pr
  !$acc declare create(t_pr)
  !$acc declare create(t_pr%Rref,t_pr%Efac_r,t_pr%cvOverR,t_pr%prefac_r, &
  !$acc                 t_pr%Cpfac_r,t_pr%K,t_pr%Zc,t_pr%Zc1,t_pr%Zc2,t_pr%a,t_pr%b)
  contains 
! Include ideal gas EoS 
#ifdef IG 
! Initialize 
  subroutine init_EOSModel()
    use mod_param
    implicit none 
    real(mytype) :: prefac_r
    integer      :: ierr                         ! ierr
    if (nrank==0) then
      if (USE_EOS=='IG') then
        write(stdout,* ) 'Correct initialisation of USE_EOS'   
      else
        write(stdout,* ) 'Mismatch USE_EOS between initBL and Makefile, check both again!'
        call decomp_2d_finalize
        call mpi_finalize(ierr) 
        stop
      endif
    endif
    ! Calculate variables
    prefac_r = 1.0_mytype
    Rhoref = 1.0_mytype
    Tref = 1.0_mytype
    ! Declare variables 
    t_ig%gam   = ig_gam
    t_ig%Rgas  = eos_Rgas
    t_ig%cv    = eos_Rgas/(ig_gam - 1.0_mytype)
    t_ig%cp    = t_ig%cv*ig_gam
    t_ig%cvInv = 1.0_mytype/t_ig%cv
    t_ig%prefac_r = prefac_r
    !$acc update device(t_ig)
    !$acc update device(t_ig%gam,t_ig%Rgas,t_ig%cp,t_ig%cv,t_ig%cvInv,t_ig%prefac_r)
  end subroutine
! Calculate EoS from density and internal energy
  subroutine calcEOS_re(rho,ien,pre,tem) 
    !$acc routine seq
    implicit none
    real(mytype), intent(IN)  :: rho,ien
    real(mytype), intent(OUT) :: pre,tem
    tem = t_ig%cvInv*ien
    pre = rho*t_ig%Rgas*tem
  end subroutine
! Calculate EoS from density and pressure
  subroutine calcEOS_rP(rho,pre,ien,tem)
    !$acc routine seq
    implicit none
    real(mytype), intent(IN)  :: rho,pre
    real(mytype), intent(OUT) :: ien,tem
    tem = pre/t_ig%Rgas/rho 
    ien = t_ig%cv*tem 
  end subroutine
! Calculate EoS from density and temperature
  subroutine calcEOS_rT(rho,tem,ien,pre)
    !$acc routine seq
    implicit none
    real(mytype), intent(IN)  :: rho,tem
    real(mytype), intent(OUT) :: ien,pre
    ien = t_ig%cv*tem
    pre = rho*t_ig%Rgas*tem
  end subroutine
! Calculate EoS from pressure and temperature
  subroutine calcEOS_PT(pre,tem,rho,ien)
    !$acc routine seq
    implicit none
    real(mytype), intent(IN)  :: pre,tem
    real(mytype), intent(OUT) :: rho,ien
    rho = pre/(tem*t_ig%Rgas)
    ien = t_ig%cv * tem
  end subroutine
! Calculate speed of sound from density and internal energy
  subroutine calcSOS_re(rho,ien,sos)
    !$acc routine seq
    implicit none
    real(mytype), intent(IN)  :: rho,ien
    real(mytype), intent(OUT) :: sos
    sos = sqrt(t_ig%gam*(t_ig%gam-1.0_mytype)*ien) 
  end subroutine
! Calculate c_p/alpha_v from density and internal energy
  subroutine calcFac_re(rho,ien,fac,tref,rhoref)
    !$acc routine seq
    implicit none
    real(mytype), intent(IN)  :: rho,ien,tref,rhoref
    real(mytype), intent(OUT) :: fac
    fac = t_ig%gam*ien
  end subroutine
#endif
! Include Van der Waals EoS 
#ifdef VdW
! Initialize 
  subroutine init_EOSModel()
    use mod_param
    implicit none 
    real(mytype) :: vdw_cvR, prefac_r, Cpfac_r
    integer      :: ierr                         ! ierr
    if (nrank==0) then
      if (USE_EOS=='VdW') then
        write(stdout,* ) 'Correct initialisation of USE_EOS'   
      else
        write(stdout,* ) 'Mismatch USE_EOS between initBL and Makefile, check both again!'
        call decomp_2d_finalize
        call mpi_finalize(ierr) 
        stop
      endif
    endif
    ! Calculate variables
    vdw_cvR=eos_dof/2
    prefac_r = vdw_Zc/Rhoref/Tref/Ec/Cpref
    Cpfac_r = Tref/Ma**2/SOSref**2/vdw_Zc
    ! Declare variables
    t_vdw%Rref= Rref
    t_vdw%a = vdw_a
    t_vdw%prefac_r = prefac_r 
    t_vdw%Efac_r = vdw_Zc/Tref/Ec/Cpref
    t_vdw%Zc = vdw_Zc
    t_vdw%cvOverR = vdw_cvR
    t_vdw%Cpfac_r = Cpfac_r
    !$acc update device(t_vdw) 
    !$acc update device(t_vdw%Rref,t_vdw%a,t_vdw%prefac_r,t_vdw%Efac_r,t_vdw%Zc,t_vdw%cvOverR,t_vdw%Cpfac_r)
  end subroutine 
! Calculate EoS from density and internal energy
  subroutine calcEOS_re(rho_r,ien,pre,tem_r)
    !$acc routine seq
    implicit none
    real(mytype), intent(IN)  :: rho_r,ien
    real(mytype), intent(OUT) :: pre,tem_r
    real(mytype) :: ien_r, pre_r
    real(mytype) :: Efac_r, prefac_r, vdw_a, Rref, cvOverR
    ! Reduced quantities
    Efac_r = t_vdw%Efac_r
    prefac_r = t_vdw%prefac_r
    vdw_a = t_vdw%a
    Rref=t_vdw%Rref
    cvOverR=t_vdw%cvOverR
    ien_r = ien/Efac_r
    tem_r = 1.0_mytype/(Rref*cvOverR)*(ien_r + vdw_a*rho_r)
    pre_r = 8.0_mytype*tem_r/(3.0_mytype/rho_r-1.0_mytype) - vdw_a*rho_r**2
    ! Non-dimensional quantities
    pre = pre_r*prefac_r
  end subroutine
! Calculate EoS from density and pressure
  subroutine calcEOS_rP(rho_r,pre,ien,tem_r)
    !$acc routine seq  
    use decomp_2d, only: mytype
    implicit none
    real(mytype), intent(IN)  :: rho_r,pre
    real(mytype), intent(OUT) :: ien,tem_r
    real(mytype) :: pre_r, ien_r
    real(mytype) :: Efac_r, prefac_r, vdw_a, Rref, cvOverR
    ! Reduced quantities
    Efac_r = t_vdw%Efac_r
    prefac_r = t_vdw%prefac_r
    vdw_a = t_vdw%a
    Rref=t_vdw%Rref
    cvOverR=t_vdw%cvOverR
    pre_r = pre/prefac_r
    tem_r = 1.0_mytype/8.0_mytype*(3.0_mytype/rho_r-1.0_mytype)*(pre_r+vdw_a*rho_r**2)
    ien_r = Rref*cvOverR*tem_r - vdw_a*rho_r
    ! Non-dimensional quantities
    ien = Efac_r*ien_r
  end subroutine
! Calculate EoS from density and temperature
  subroutine calcEOS_rT(rho_r,tem_r,ien,pre)
    !$acc routine seq  
    use decomp_2d, only: mytype
    implicit none
    real(mytype), intent(IN)  :: rho_r,tem_r
    real(mytype), intent(OUT) :: ien,pre
    real(mytype) :: pre_r, ien_r
    real(mytype) :: Efac_r, prefac_r, vdw_a, Rref, cvOverR
    ! Reduced quantities
    Efac_r = t_vdw%Efac_r
    prefac_r = t_vdw%prefac_r
    vdw_a = t_vdw%a
    Rref=t_vdw%Rref
    cvOverR=t_vdw%cvOverR
    ien_r = Rref*cvOverR*tem_r - vdw_a*rho_r
    pre_r = 8.0_mytype*tem_r/(3.0_mytype/rho_r-1.0_mytype) - vdw_a*rho_r**2
    ! Non-dimensional quantities
    pre = pre_r*prefac_r
    ien = Efac_r*ien_r
  end subroutine
! Calculate EoS from pressure and temperature
  subroutine calcEOS_PT(pre,tem_r,rho_r,ien)
    !$acc routine seq  
    use decomp_2d, only: mytype
    use mod_math
    implicit none
    real(mytype), intent(IN)  :: pre,tem_r
    real(mytype), intent(OUT) :: rho_r,ien
    real(mytype) :: pre_r,v_r,ien_r,A,B,C,D
    real(mytype) :: Efac_r, prefac_r, vdw_a, Rref, cvOverR
    ! Reduced quantities
    Efac_r = t_vdw%Efac_r
    prefac_r = t_vdw%prefac_r
    vdw_a = t_vdw%a
    Rref=t_vdw%Rref
    cvOverR=t_vdw%cvOverR
    pre_r = pre/prefac_r
    A =  3.0_mytype*pre_r
    B = -(pre_r+8*tem_r)
    C =  9.0_mytype
    D = -3.0_mytype
    call cubic_root(A,B,C,D,v_r)
    rho_r = 1/v_r
    ien_r = Rref*cvOverR*tem_r - vdw_a*rho_r
    ! Non-dimensional quantities
    ien = Efac_r*ien_r
  end subroutine
! Calculate speed of sound from density and internal energy
  subroutine calcSOS_re(rho_r,ien,sos)
    !$acc routine seq
    use decomp_2d, only: mytype
    implicit none
    real(mytype), intent(IN)  :: rho_r,ien
    real(mytype), intent(OUT) :: sos
    real(mytype) :: sosfac_r,ien_r,sos_r,cvOverR
    ! Reduced quantities
    sosfac_r = t_vdw%Efac_r ! sosfac_r is equal to Efac_r
    cvOverR = t_vdw%cvOverR
    ien_r = ien/t_vdw%Efac_r
    sos_r = (1.0_mytype+1.0_mytype/cvOverR)/cvOverR*(ien_r+3.0_mytype*rho_r)*(3.0_mytype/(3.0_mytype-rho_r))**2 - 6.0_mytype*rho_r
    ! Non-dimensional quantities
    sos   = sqrt(sosfac_r*sos_r)
  end subroutine
! Calculate c_p/alpha_v from density and internal energy
  subroutine calcFac_re(rho_r,ien,fac,tref,rhoref)
    !$acc routine seq
    use decomp_2d, only: mytype
    implicit none
    real(mytype), intent(IN)  :: rho_r,ien,tref,rhoref
    real(mytype), intent(OUT) :: fac
    real(mytype) :: ien_r,tem_r,cp_r,cp
    real(mytype) :: dPdT_rho_r,dPdrho_T_r,dPdT_rho,dPdrho_T,alpha_v
    real(mytype) :: prefac_r, vdw_a, cvOverR
    ! Reduced quantities
    vdw_a   = t_vdw%a
    cvOverR = t_vdw%cvOverR
    ien_r   = ien/t_vdw%Efac_r
    tem_r   = t_vdw%Zc/cvOverR*(ien_r + vdw_a*rho_r)
    cp_r    = cvOverR+1.0_mytype/(1.0_mytype-((3.0_mytype/rho_r-1.0_mytype)**2/(4.0_mytype*tem_r*(1.0_mytype/rho_r)**3)))
    dPdT_rho_r = 8.0_mytype*rho_r/(3.0_mytype-rho_r)
    dPdrho_T_r = 24.0_mytype*tem_r/(3.0_mytype-rho_r)**2-2.0_mytype*vdw_a*rho_r
    ! Non-dimensional quantities
    cp      = cp_r*t_vdw%Cpfac_r
    dPdT_rho = dPdT_rho_r*t_vdw%prefac_r*tref
    dPdrho_T = dPdrho_T_r*t_vdw%prefac_r*rhoref
    alpha_v = 1.0_mytype*rhoref/rho_r*dPdT_rho/dPdrho_T
    fac     = cp/alpha_v
  end subroutine
#endif
! Include Peng-Robinson EoS 
#ifdef PR
! Initialize 
  subroutine init_EOSModel()
    use mod_param
    implicit none
    real(mytype) :: pr_cvR, Cpfac_r
    integer      :: ierr                         ! ierr
    if (nrank==0) then
      if (USE_EOS=='PR') then
        write(stdout,* ) 'Correct initialisation of USE_EOS'   
      else
        write(stdout,* ) 'Mismatch USE_EOS between initBL and Makefile, check both again!'
        call decomp_2d_finalize
        call mpi_finalize(ierr) 
        stop
      endif
    endif
    ! Calculate variables
    pr_cvR=eos_dof/2
    Cpfac_r = Tref/Ma**2/SOSref**2/pr_Zc
    ! Declare variables
    t_pr%K   = 0.37464_mytype+1.54226_mytype*eos_ac-0.26992_mytype*eos_ac**2
    t_pr%Zc  = pr_Zc
    t_pr%Zc1 = pr_Zc**(-1)
    t_pr%Zc2 = pr_Zc**(-2)
    t_pr%a = pr_a
    t_pr%b = pr_b
    t_pr%prefac_r = pr_Zc/Rhoref/Tref/Ec/Cpref
    t_pr%Efac_r   = pr_Zc/Tref/Ec/Cpref
    t_pr%cvOverR  = pr_cvR
    t_pr%Cpfac_r  = Cpfac_r
    t_pr%Rref= Rref
    !$acc update device(t_pr) 
    !$acc update device(t_pr%Rref,t_pr%K,t_pr%Zc,t_pr%Zc1,t_pr%Zc2,t_pr%a,t_pr%b, &
    !$acc               t_pr%prefac_r,t_pr%Efac_r,t_pr%cvOverR,t_pr%Cpfac_r)
  end subroutine
! Calculate EoS from density and internal energy
  subroutine calcEOS_re(rho_r,ien,pre,tem_r)
    !$acc routine seq
    use decomp_2d, only: mytype
    implicit none
    real(mytype), intent(IN) :: rho_r,ien
    real(mytype), intent(OUT) :: pre,tem_r
    real(mytype) :: T1, F1, F2, F3, alpha
    real(mytype) :: ien_r, pre_r
    real(mytype) :: Efac_r, prefac_r, cvOverR, t_pr_a, t_pr_b, t_pr_Zc1, t_pr_Zc2, t_pr_K, sqrt2
    ! Reduced quantities
    Efac_r = t_pr%Efac_r
    prefac_r = t_pr%prefac_r
    cvOverR = t_pr%cvOverR
    t_pr_a = t_pr%a
    t_pr_b = t_pr%b
    t_pr_Zc1 = t_pr%Zc1
    t_pr_Zc2 = t_pr%Zc2
    t_pr_K = t_pr%K
    sqrt2=2**(0.5_mytype)
    ien_r=ien/Efac_r
    T1 = log( (1+t_pr_b*t_pr_Zc1*rho_r*(1-sqrt2))/(1+t_pr_b*t_pr_Zc1*rho_r*(1+sqrt2)) ) 
    F1 = ( t_pr_a*t_pr_Zc1*(t_pr_K+1)**2*T1/(2*sqrt2*t_pr_b) - ien_r )
    F2 = ( t_pr_a*t_pr_Zc1*t_pr_K*(t_pr_K+1)*T1/(2*sqrt2*t_pr_b) )
    F3 = ( cvOverR*t_pr_Zc1 )     
    tem_r = ( (F2 + sqrt( F2**2-4*F1*F3 ))/2/F3 )**2
    alpha = ( 1 + t_pr_K*(1-sqrt(tem_r)) )**2
    pre_r = tem_r*t_pr_Zc1/( 1.0_mytype/rho_r-t_pr_b*t_pr_Zc1 )- &
            alpha*t_pr_a*t_pr_Zc2/(1.0_mytype/rho_r*(1.0_mytype/rho_r+2*t_pr_b*t_pr_Zc1)-t_pr_b**2*t_pr_Zc2 )
    ! Non-dimensional quantities
    pre = pre_r*prefac_r
  end subroutine
! Calculate EoS from density and pressure
  subroutine calcEOS_rP(rho_r,pre,ien,tem_r)
    !$acc routine seq
    use decomp_2d, only: mytype
    implicit none
    real(mytype), intent(IN) :: rho_r,pre
    real(mytype), intent(OUT) :: ien,tem_r
    real(mytype) :: N1, N2, F1, F2, F3
    real(mytype) :: pre_r, ien_r
    real(mytype) :: Efac_r, prefac_r, cvOverR, t_pr_a, t_pr_b, t_pr_Zc1, t_pr_Zc2, t_pr_K, sqrt2
    ! Reduced quantities
    Efac_r = t_pr%Efac_r
    prefac_r = t_pr%prefac_r
    cvOverR = t_pr%cvOverR
    t_pr_a = t_pr%a
    t_pr_b = t_pr%b
    t_pr_Zc1 = t_pr%Zc1
    t_pr_Zc2 = t_pr%Zc2
    t_pr_K = t_pr%K
    sqrt2=2**(0.5_mytype)
    pre_r = pre/prefac_r
    N1 = ( 1.0_mytype/rho_r-t_pr_b*t_pr_Zc1 )
    N2 = ( (1.0_mytype/rho_r)*(1.0_mytype/rho_r+2*t_pr_b*t_pr_Zc1)-t_pr_b**2*t_pr%Zc2 )
    F1 = ( (1.0_mytype + 2.0_mytype*t_pr_K + t_pr_K**2)*t_pr_a*t_pr_Zc2/N2 + pre_r )
    F2 = ( t_pr_a*t_pr_Zc2*(2.0_mytype*t_pr_K+2.0_mytype*t_pr_K**2)/N2 )
    F3 = ( t_pr_K**2*t_pr_a*t_pr_Zc2/N2 - t_pr_Zc1/N1 )
    tem_r = ( (F2 - sqrt( F2**2-4*F1*F3 ))/2/F3 )**2
    ien_r = ( cvOverR*tem_r*t_pr_Zc1 + t_pr_a*t_pr_Zc1/(2.0_mytype*sqrt2*t_pr_b) &
            *( (t_pr_K+1.0_mytype)**2-t_pr_K*(t_pr_K+1.0_mytype)*sqrt(tem_r) ) &
            *log( (1.0_mytype+t_pr_b*t_pr_Zc1*rho_r*(1.0_mytype-sqrt2))/(1.0_mytype+t_pr_b*t_pr_Zc1*rho_r*(1.0_mytype+sqrt2)) ) )
    ! Non-dimensional quantities
    ien = Efac_r*ien_r
  end subroutine
! Calculate EoS from density and temperature
  subroutine calcEOS_rT(rho_r,tem_r,ien,pre)
    !$acc routine seq
    use decomp_2d, only: mytype
    implicit none
    real(mytype), intent(IN) :: rho_r,tem_r
    real(mytype), intent(OUT) :: ien,pre
    real(mytype) :: alpha
    real(mytype) :: pre_r, ien_r
    real(mytype) :: Efac_r, prefac_r, cvOverR, t_pr_a, t_pr_b, t_pr_Zc1, t_pr_Zc2, t_pr_K, sqrt2
    ! Reduced quantities
    Efac_r = t_pr%Efac_r
    prefac_r = t_pr%prefac_r
    cvOverR = t_pr%cvOverR
    t_pr_a = t_pr%a
    t_pr_b = t_pr%b
    t_pr_Zc1 = t_pr%Zc1
    t_pr_Zc2 = t_pr%Zc2
    t_pr_K = t_pr%K
    sqrt2=2**(0.5_mytype)
    alpha = ( 1.0_mytype + t_pr_K*(1-sqrt(tem_r)) )**2
    ien_r = ( cvOverR*tem_r*t_pr_Zc1 + t_pr_a*t_pr_Zc1/(2.0_mytype*sqrt2*t_pr%b) &
            *( (t_pr%K+1.0_mytype)**2-t_pr_K*(t_pr_K+1.0_mytype)*sqrt(tem_r) ) &
            *log( (1.0_mytype+t_pr_b*t_pr_Zc1*rho_r*(1.0_mytype-sqrt2))/(1.0_mytype+t_pr_b*t_pr_Zc1*rho_r*(1.0_mytype+sqrt2)) ) )
    pre_r = tem_r*t_pr_Zc1/( 1.0_mytype/rho_r-t_pr%b*t_pr_Zc1 )&
            - alpha*t_pr_a*t_pr_Zc2/( (1.0_mytype/rho_r)*(1.0_mytype/rho_r+2.0_mytype*t_pr_b*t_pr_Zc1)-t_pr_b**2*t_pr_Zc2 )
    ! Non-dimensional quantities
    pre = pre_r*prefac_r
    ien = Efac_r*ien_r
  end subroutine
! Calculate EoS from pressure and temperature
  subroutine calcEOS_PT(pre,tem_r,rho_r,ien)
    !$acc routine seq
    use decomp_2d
    implicit none
    integer :: ierr
    real(mytype), intent(IN) :: pre,tem_r
    real(mytype), intent(OUT) :: rho_r,ien
    write(*,*) "Not included yet! Select different boundary conditions"
    call decomp_2d_finalize
    call mpi_finalize(ierr)
    stop
  end subroutine
! Calculate speed of sound from density and internal energy
  subroutine calcSOS_re(rho_r,ien,sos)
    !$acc routine seq
    use decomp_2d, only: mytype
    implicit none
    real(mytype), intent(IN)  :: rho_r,ien
    real(mytype), intent(OUT) :: sos
    real(mytype) :: T1, F1, F2, F3, tem_r, cv_r, dPdrho_T_r, dPdT_rho_r, alpha
    real(mytype) :: ien_r,sos_r
    real(mytype) :: sosfac_r, cvOverR, t_pr_a, t_pr_b, t_pr_Zc1, t_pr_Zc2, t_pr_K, sqrt2
    ! Reduced quantities
    sosfac_r = t_pr%Efac_r  ! sosfac_r is equal to Efac_r
    ien_r = ien/t_pr%Efac_r
    cvOverR = t_pr%cvOverR
    t_pr_a = t_pr%a
    t_pr_b = t_pr%b
    t_pr_Zc1 = t_pr%Zc1
    t_pr_Zc2 = t_pr%Zc2
    t_pr_K = t_pr%K
    sqrt2=2**(0.5_mytype)
    T1 = log( (1.0_mytype+t_pr_b*t_pr_Zc1*rho_r*(1.0_mytype-sqrt2))/(1.0_mytype+t_pr_b*t_pr_Zc1*rho_r*(1.0_mytype+sqrt2)) ) 
    F1 = ( t_pr_a*t_pr_Zc1*(t_pr_K+1.0_mytype)**2*T1/(2.0_mytype*sqrt2*t_pr_b) - ien_r )
    F2 = ( t_pr_a*t_pr_Zc1*t_pr_K*(t_pr_K+1.0_mytype)*T1/(2.0_mytype*sqrt2*t_pr_b) )
    F3 = ( cvOverR*t_pr_Zc1 )
    tem_r = ( (F2 + sqrt( F2**2-4.0_mytype*F1*F3 ))/2/F3 )**2
    alpha = ( 1.0_mytype + t_pr_K*(1.0_mytype-sqrt(tem_r)) )**2
    cv_r = ( cvOverR-(t_pr_a*t_pr_K*(t_pr_K+1.0_mytype))/(4.0_mytype*t_pr_b*sqrt(2.0_mytype*tem_r)) &
      *log( (1.0_mytype+t_pr_b*t_pr_Zc1*rho_r*(1.0-sqrt2))/(1.0_mytype+t_pr_b*t_pr_Zc1*rho_r*(1.0_mytype+sqrt2)) ) )
    dPdrho_T_r = ( t_pr_Zc1*tem_r/(t_pr_Zc1*t_pr_b*rho_r-1.0_mytype)**2 &
               - t_pr_a*t_pr_Zc2*alpha*(2.0_mytype*rho_r+2.0_mytype*t_pr_b*t_pr_Zc1*rho_r**2) &
               /(1.0_mytype+2.0_mytype*t_pr_b*t_pr_Zc1*rho_r-rho_r**2*t_pr_b**2*t_pr_Zc2)**2 )
    dPdT_rho_r = ( t_pr_K*sqrt(alpha/tem_r)*(rho_r**2*t_pr_a*t_pr_Zc2) & 
               /(1.0_mytype+2.0_mytype*t_pr_b*t_pr_Zc1*rho_r-t_pr_b**2*t_pr_Zc2*rho_r**2) & 
               - (rho_r*t_pr_Zc1)/(rho_r*t_pr_b*t_pr_Zc1-1.0_mytype) )
    sos_r = dPdrho_T_r+t_pr%Zc*tem_r/(rho_r**2*cv_r)*dPdT_rho_r**2
    ! Non-dimensional quantities
    sos = sqrt(sosfac_r*sos_r)
  end subroutine
! Calculate c_p/alpha_v from density and internal energy
  subroutine calcFac_re(rho_r,ien,fac,tref,rhoref)
    !$acc routine seq
    use decomp_2d, only: mytype
    implicit none
    real(mytype), intent(IN)  :: rho_r,ien,tref,rhoref
    real(mytype), intent(OUT) :: fac
    real(mytype) :: ien_r,tem_r,cp_r,dPdT_rho_r,dPdrho_T_r
    real(mytype) :: T1, F1, F2, F3, alpha, alpha_v, cv_r
    real(mytype) :: cp,dPdT_rho,dPdrho_T
    real(mytype) :: cvOverR, t_pr_a, t_pr_b, t_pr_Zc1, t_pr_Zc2, t_pr_K, sqrt2
    ! Reduced quantities
    cvOverR = t_pr%cvOverR
    t_pr_a = t_pr%a
    t_pr_b = t_pr%b
    t_pr_Zc1 = t_pr%Zc1
    t_pr_Zc2 = t_pr%Zc2
    t_pr_K = t_pr%K
    sqrt2=2**(0.5_mytype)
    ien_r = ien/t_pr%Efac_r
    T1 = log( (1.0_mytype+t_pr_b*t_pr_Zc1*rho_r*(1.0_mytype-sqrt2))/(1.0_mytype+t_pr_b*t_pr_Zc1*rho_r*(1.0_mytype+sqrt2)) ) 
    F1 = ( t_pr_a*t_pr_Zc1*(t_pr_K+1.0_mytype)**2*T1/(2.0_mytype*sqrt2*t_pr_b) - ien_r )
    F2 = ( t_pr_a*t_pr_Zc1*t_pr_K*(t_pr_K+1.0_mytype)*T1/(2.0_mytype*sqrt2*t_pr_b) )
    F3 = ( cvOverR*t_pr_Zc1 )
    tem_r = ( (F2 + sqrt( F2**2-4.0_mytype*F1*F3 ))/2/F3 )**2
    alpha = ( 1.0_mytype + t_pr_K*(1-sqrt(tem_r)) )**2
    cv_r  = ( cvOverR-(t_pr_a*t_pr_K*(t_pr_K+1.0_mytype))/(4.0_mytype*t_pr_b*sqrt(2.0_mytype*tem_r))&
            *log( (1.0_mytype+t_pr_b*t_pr_Zc1*rho_r*(1.0_mytype-sqrt2))/(1.0_mytype+t_pr_b*t_pr_Zc1*rho_r*(1.0_mytype+sqrt2)) ) )
    dPdrho_T_r = ( t_pr_Zc1*tem_r/(t_pr_Zc1*t_pr_b*rho_r-1.0_mytype)**2 &
               - t_pr_a*t_pr_Zc2*alpha*(2.0_mytype*rho_r+2.0_mytype*t_pr_b*t_pr_Zc1*rho_r**2) &
               /(1.0_mytype+2.0_mytype*t_pr_b*t_pr_Zc1*rho_r-rho_r**2*t_pr_b**2*t_pr_Zc2)**2 )
    dPdT_rho_r = ( t_pr_K*sqrt(alpha/tem_r)*(rho_r**2*t_pr_a*t_pr_Zc2) &
               /(1.0_mytype+2.0_mytype*t_pr_b*t_pr_Zc1*rho_r-t_pr_b**2*t_pr_Zc2*rho_r**2) & 
               - (rho_r*t_pr_Zc1)/(rho_r*t_pr_b*t_pr_Zc1-1.0_mytype) )
    cp_r = ( cv_r+tem_r/rho_r**2*t_pr%Zc*dPdT_rho_r**2/dPdrho_T_r )
    ! Non-dimensional quantities
    cp   = cp_r*t_pr%Cpfac_r
    dPdT_rho = dPdT_rho_r*t_pr%prefac_r*tref
    dPdrho_T = dPdrho_T_r*t_pr%prefac_r*rhoref
    alpha_v = 1.0_mytype*rhoref/rho_r*dPdT_rho/dPdrho_T
    fac     = cp/alpha_v
  end subroutine
#endif  
end module mod_eos_var
