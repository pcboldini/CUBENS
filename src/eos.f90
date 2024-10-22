! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -

module mod_eos
  use decomp_2d, only: mytype
  !@acc use openacc
  implicit none

! Variables container for global variables  
! EoS: Equation of State, TP: transport properties 
  type t_INIT_PARAM
    real(mytype) :: Rhoref, Tref
  end type 
  type(t_INIT_PARAM) :: t_param
  !$acc declare create(t_param)
  !$acc declare create(t_param%Rhoref,t_param%Tref)
  contains 

! Initialize global variables
  subroutine init_PARAM_EOS()
    use mod_param
    implicit none 
    t_param%Rhoref = Rhoref
    t_param%Tref = Tref
    !$acc update device(t_param) 
    !$acc update device(t_param%Rhoref,t_param%Tref)
  end subroutine init_PARAM_EOS


! Call secondary EoS and TP variables from density and internal energy
  subroutine calcState_RE(rho,ien,pre,tem,mu,ka,i1,i2,j1,j2,k1,k2)
    use mod_param, only: nHalo
    use mod_eos_visc
    use mod_eos_var
    implicit none
    integer :: i,j,k
    integer, intent(IN) :: i1,i2,j1,j2,k1,k2
    real(mytype), intent(IN),  dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,ien
    real(mytype), intent(OUT), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: pre,tem,mu,ka
    !$acc parallel loop collapse(3) default(present)
    do k=k1,k2
      do j=j1,j2
        do i=i1,i2
          call calcEOS_re(rho(i,j,k)*t_param%Rhoref, ien(i,j,k), pre(i,j,k), tem(i,j,k))
          call calcVisc(tem(i,j,k), rho(i,j,k)*t_param%Rhoref, mu(i,j,k), ka(i,j,k))
          tem(i,j,k) = tem(i,j,k) / t_param%Tref
          enddo 
      enddo 
    enddo 
  end subroutine


! Call secondary EoS and TP variables from density and pressure
  subroutine calcState_rP(rho,pre,ien,tem,mu,ka,i1,i2,j1,j2,k1,k2)
    use mod_param, only: nHalo
    use mod_eos_visc
    use mod_eos_var
    implicit none
    integer :: i,j,k
    integer, intent(IN) :: i1,i2,j1,j2,k1,k2
    real(mytype), intent(IN),  dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,pre
    real(mytype), intent(OUT), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: ien,tem,mu,ka
    !$acc parallel loop collapse(3) default(present)
    do k=k1,k2
      do j=j1,j2
        do i=i1,i2
          call calcEOS_rP(rho(i,j,k)*t_param%Rhoref, pre(i,j,k), ien(i,j,k), tem(i,j,k))
          call calcVisc(tem(i,j,k), rho(i,j,k)*t_param%Rhoref, mu(i,j,k), ka(i,j,k))    
          tem(i,j,k) = tem(i,j,k) / t_param%Tref   
          enddo 
      enddo 
    enddo 
  end subroutine


! Call secondary EoS and TP variables from density and temperature
  subroutine calcState_rT(rho,tem,ien,pre,mu,ka, i1,i2,j1,j2,k1,k2)
    use mod_param, only: nHalo
    use mod_eos_visc
    use mod_eos_var
    implicit none
    integer :: i,j,k
    integer, intent(IN) :: i1,i2,j1,j2,k1,k2
    real(mytype), intent(IN),  dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,tem
    real(mytype), intent(OUT), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: ien,pre,mu,ka
    !$acc parallel loop collapse(3) default(present)
    do k=k1,k2
      do j=j1,j2
        do i=i1,i2
          call calcEOS_rT(rho(i,j,k)*t_param%Rhoref, tem(i,j,k)*t_param%Tref, ien(i,j,k), pre(i,j,k)) 
          call calcVisc(tem(i,j,k)*t_param%Tref, rho(i,j,k)*t_param%Rhoref, mu(i,j,k), ka(i,j,k)) 
          enddo 
      enddo 
    enddo 
  end subroutine


! Call secondary EoS and TP variables from pressure and temperature
  subroutine calcState_PT(pre,tem,rho,ien,mu,ka, i1,i2,j1,j2,k1,k2)
    use mod_param, only: nHalo
    use mod_eos_visc
    use mod_eos_var
    implicit none
    integer :: i,j,k
    integer, intent(IN) :: i1,i2,j1,j2,k1,k2
    real(mytype), intent(IN),  dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: pre,tem
    real(mytype), intent(OUT), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,ien,mu,ka
    !$acc parallel loop collapse(3) default(present)
    do k=k1,k2
      do j=j1,j2
        do i=i1,i2
          call calcEOS_PT(pre(i,j,k), tem(i,j,k)*t_param%Tref, rho(i,j,k), ien(i,j,k))
          call calcVisc(tem(i,j,k)*t_param%Tref, rho(i,j,k), mu(i,j,k), ka(i,j,k)) 
          rho(i,j,k) = rho(i,j,k) / t_param%Rhoref 
          enddo 
      enddo 
    enddo 
  end subroutine

  
! Call speed of sound from density and internal energy
  subroutine calcSOS(rho,ien,sos) 
    !$acc routine seq
    use mod_param, only: nHalo
    use mod_eos_var
    implicit none
    real(mytype), intent(IN)  :: rho,ien
    real(mytype), intent(OUT) :: sos
    call calcSOS_re(rho*t_param%Rhoref, ien, sos)
  end subroutine
  
! Call c_p/alpha_v from density and internal energy (for boundary conditions)
  subroutine calcFac(rho,ien,fac) 
    !$acc routine seq
    use mod_param, only: nHalo
    use mod_eos_var
    implicit none
    real(mytype), intent(IN)  :: rho,ien
    real(mytype), intent(OUT) :: fac
    call calcFac_re(rho*t_param%Rhoref, ien, fac, t_param%Tref, t_param%Rhoref)
  end subroutine
end module mod_eos
