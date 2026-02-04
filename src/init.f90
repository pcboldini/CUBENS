! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
! initial conditions module

module mod_init
use io_std_units
implicit none
integer, parameter :: typeReadBlasius = 8
contains


! boundary layer: Blasius solution loaded from external files
subroutine initField_BL(part,xcoord,zcoord,rho,u,v,w,ien,pre,tem,mu,ka)
  use decomp_2d
  use mod_param
  use mod_eos
  use mod_eos_var
  use mod_grid
  use mod_halo
  use mod_math
  implicit none
  integer :: i,j,k, npts, fs
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka
  real(8),      allocatable, dimension(:) :: xRead, wRead, uRead, rRead
  real(mytype), allocatable, dimension(:) ::        wIntp, uIntp, rIntp
  real(mytype), dimension(:) :: xcoord, zcoord
  real(mytype) :: scaling_delta,rho_integ
  TYPE (DECOMP_INFO), intent(IN) :: part
  ! wall-normal mesh
  inquire(FILE='preproc/initBL/inputDNS/prof_x.bin', SIZE=fs); npts = fs/8
  allocate(xRead(npts))
  allocate(rRead(npts)); allocate(rIntp(npts))
  allocate(wRead(npts)); allocate(wIntp(npts))
  allocate(uRead(npts)); allocate(uIntp(npts))
  ! load mesh, density, streamwise velocity, and wall-normal velocity from files
  open(9,file='preproc/initBL/inputDNS/prof_x.bin',access='direct',recl=fs); read(9,rec=1) (xRead(i),i=1,npts);close(9)
  open(9,file='preproc/initBL/inputDNS/prof_r.bin',access='direct',recl=fs); read(9,rec=1) (rRead(i),i=1,npts);close(9)
  open(9,file='preproc/initBL/inputDNS/prof_w.bin',access='direct',recl=fs); read(9,rec=1) (wRead(i),i=1,npts);close(9)
  open(9,file='preproc/initBL/inputDNS/prof_u.bin',access='direct',recl=fs); read(9,rec=1) (uRead(i),i=1,npts);close(9)
  call spline(xRead, rRead, npts, rIntp)
  call spline(xRead, wRead, npts, wIntp)
  call spline(xRead, uRead, npts, uIntp)
  ! interpolating and rescaling over flat-plate length
  do k=1,part%xsz(3)
    do j=1,part%xsz(2)
      do i=1,part%xsz(1)
        scaling_delta = sqrt((zcoord(k)+zStartDNS)/zStartDNS) ! scaling factor for Blasius thickness
        call splint(xRead, rRead, rIntp, npts, xcoord(i)/scaling_delta, rho(i,j,k))
        call splint(xRead, wRead, wIntp, npts, xcoord(i)/scaling_delta,   w(i,j,k))
        call splint(xRead, uRead, uIntp, npts, xcoord(i)/scaling_delta,   u(i,j,k))
        u(i,j,k) = u(i,j,k)/scaling_delta/Redelta_start ! wall-normal velocity has to be rescaled again
      enddo
    enddo
  enddo
  ! spanwise velocity
  v = 0.0_mytype
  ! pressure
#if defined(IG)
  pre = Pref*t_ig%prefac_r  
#elif defined(VdW)
  pre = Pref*t_vdw%prefac_r 
#elif defined(RK)
  pre = Pref*t_rk%prefac_r 
#elif defined(PR)
  pre = Pref*t_pr%prefac_r 
#endif
  ! pressure correction if buoyancy is active (self-similarity is lost)
  do k=1,part%xsz(3)
     do j=1,part%xsz(2)
        rho_integ = 0.0_mytype
        do i=part%xsz(1)-1,1,-1
           rho_integ = rho_integ + (rho(i,j,k)+rho(i+1,j,k)-2.0_mytype)*0.5_mytype*(xcoord(i)-xcoord(i+1))
           pre(i,j,k) = pre(i,j,k) - Ri_unit*rho_integ
        enddo
     enddo
  enddo
  !$acc update device(rho,pre,u,v,w)
  ! calculation of the secondary variables
  call calcState_rP(rho,pre,ien,tem,mu,ka,1,part%xsz(1),1,part%xsz(2),1,part%xsz(3)) 
  !$acc update host(rho,u,v,w,ien,pre,tem,mu,ka)
  deallocate(xRead)
  deallocate(rRead); deallocate(rIntp)
  deallocate(wRead); deallocate(wIntp)
  deallocate(uRead); deallocate(uIntp)
end subroutine


! 1-D wave for advection test case
subroutine initField_1D(part,rho,u,v,w,ien,pre,tem,mu,ka)
  use decomp_2d
  use mod_param
  use mod_eos
  use mod_eos_var
  use mod_grid
  use mod_math
  implicit none
  integer :: i,j,k,jj,kk
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka
  real(mytype) :: fz, rho_max, rho_min
  TYPE (DECOMP_INFO), intent(IN) :: part
  ! 1-D velocity and max/min density
  rho_max = 8.2984_mytype
  rho_min = 0.5967_mytype
  u   = 0.0_mytype
  v   = 0.0_mytype
  w   = 1.0_mytype
  do k=1,part%xsz(3)
    do j=1,part%xsz(2)
      do i=1,part%xsz(1)
        kk = k + xstart(3) - 1
        fz = len_z*(kk-0.5_mytype)/(nz_global)
        rho(i,j,k) = 0.5_mytype*(rho_max+rho_min) + 0.5_mytype*(rho_max-rho_min)*(sin(2.0_mytype*pi_const*fz))
      enddo
    enddo
  enddo
  ! pressure
#if defined(IG)
  pre = Pref*t_ig%prefac_r  
#elif defined(VdW)
  pre = Pref*t_vdw%prefac_r 
#elif defined(RK)
  pre = Pref*t_rk%prefac_r 
#elif defined(PR)
  pre = Pref*t_pr%prefac_r 
#endif
  !$acc update device(rho,pre,u,v,w)
  ! calculation of the secondary variables
  call calcState_rP(rho,pre,ien,tem,mu,ka, 1,xsize(1),1,xsize(2),1,xsize(3))
  !$acc update host(rho,u,v,w,ien,pre,tem,mu,ka)
end subroutine


! turbulent channel: Poiseuille laminar profile plus vortex pair according to Henningson & Kim, JFM 228, 1991
subroutine initField_CHA(part,rho,u,v,w,ien,pre,tem,mu,ka)
  use decomp_2d
  use mod_param
  use mod_eos
  use mod_eos_var
  use mod_grid
  use mod_finitediff
  implicit none
  integer :: i,j,k,jj,kk
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka
  real(mytype) :: c_0, c_1, c_2, xa, ya, za, psi, psi_y, psi_x
  TYPE (DECOMP_INFO), intent(IN) :: part
  ! integration constants for laminar velocity profile
  c_0 = 0.5_mytype*dpdz
  c_1 = dpdz 
  c_2 = 0.0_mytype
  ! density is initialized with one
  rho = 1.0_mytype
  u   = 0.0_mytype
  v   = 0.0_mytype
  w   = 0.0_mytype
  ! streamfunction psi, origin of the vortex pair is located in the domain center 
  do k=1,part%xsz(3)
    do j=1,part%xsz(2)
      do i=1,part%xsz(1)
        xa = x(i)
        w(i,j,k) = - c_0*xa**2 + c_1*xa + c_2
        xa = 2*x(i)/len_x - 1
        ya = (y(j)-0.5*len_y)*2/len_x
        za = (z(k)-0.5*len_z)*2/len_x
        psi      = (1-xa**2)**2 * ya*exp(-16*za**2-4*ya**2)
        psi_y    = (1-xa**2)**2 *exp(-16*za**2-4*ya**2)*(1 - 8*ya**2)
        psi_x    = -4*xa*(1-xa**2)*ya*exp(-16*za**2-4*ya**2) 
        u(i,j,k) = psi_y
        v(i,j,k) = -psi_x
      enddo
    enddo
  enddo
  ! linear temperature profile for the initialization of the temperature in case Twall_top not equal to Twall_bot
  do i=1,xsize(1)
    tem(i,:,:) = Twall_bot + 0.5_mytype*x(i)*(Twall_top - Twall_bot)
  enddo
  !$acc update device(rho,tem,u,v,w)
  ! calculation of the secondary variables
  call calcState_rT(rho,tem,ien,pre,mu,ka,1,part%xsz(1),1,part%xsz(2),1,part%xsz(3))
  !$acc update host(rho,u,v,w,ien,pre,tem,mu,ka)
end subroutine


! Taylor-Green Vortex: laminar-turbulent transition of a decaying vortex
subroutine initField_TGV(part,rho,u,v,w,ien,pre,tem,mu,ka)
  use decomp_2d
  use mod_param
  use mod_eos
  use mod_grid
  implicit none
  integer :: i,j,k,jj,kk
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka
  real(mytype) :: fx,fy,fz, V0,P0,R0,T0
  TYPE (DECOMP_INFO), intent(IN) :: part
  ! initialization parameters, an ideal pressure is imposed
  V0 = 1.0_mytype
  T0 = 1.0_mytype
  R0 = 1.0_mytype
  P0 = R0*eos_Rgas*T0
  do i=1,xsize(1)
    fx = (i - 0.5_mytype)/nx_global*len_x
    do j=1,xsize(2)
      jj = j + xstart(2) - 1
      fy = (jj - 0.5_mytype)/ny_global*len_y
      do k=1,xsize(3)
        kk = k + xstart(3) - 1
        fz = (kk - 0.5_mytype)/nz_global*len_z
        u(i,j,k) =  V0*sin(fx)*cos(fy)*cos(fz)
        v(i,j,k) = -V0*cos(fx)*sin(fy)*cos(fz)
        w(i,j,k) =  0.0
        pre(i,j,k) = P0 + 1.0/16.0*R0*V0**2 * (cos(2.0*fx) + cos(2.0*fy)) * (cos(2.0*fz) + 2.0)
        rho(i,j,k) = R0
      enddo
    enddo
  enddo
!$acc update device(rho,pre,u,v,w)
! calculation of the secondary variables, EoS and TP are here considered
call calcState_rP(rho,pre,ien,tem,mu,ka,1,part%xsz(1),1,part%xsz(2),1,part%xsz(3)) 
!$acc update host(rho,u,v,w,ien,pre,tem,mu,ka)
end subroutine


end module mod_init
