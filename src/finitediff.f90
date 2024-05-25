! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_finitediff
  use decomp_2d
  use mod_param
  use mod_grid
  implicit none
! definition of the first and second derivative coefficients
  real(mytype), target, dimension(1)    :: c1_CD2 = (/ 1.0_mytype/2.0_mytype                     /)
  real(mytype), target, dimension(1:2)  :: c1_CD4 = (/ 2.0_mytype/3.0_mytype, -1.0_mytype/12.0_mytype          /)
  real(mytype), target, dimension(1:3)  :: c1_CD6 = (/ 0.75_mytype,    -0.15_mytype,    1.0_mytype/60.0_mytype /)

  real(mytype), target, dimension( 0:1) :: c2_CD2 = (/-2.0_mytype,       1.0_mytype /)
  real(mytype), target, dimension( 0:2) :: c2_CD4 = (/- (5.0_mytype/2.0_mytype)  , &
                                                        (4.0_mytype/3.0_mytype)  , &
                                                      - (1.0_mytype/12.0_mytype)   /)
  real(mytype), target, dimension( 0:3) :: c2_CD6 = (/- (49.0_mytype/18.0_mytype), & 
                                                        (3.0_mytype/2.0_mytype)  , &
                                                        (-3.0_mytype/20.0_mytype), & 
                                                        (1.0_mytype/90.0_mytype) /)
  real(mytype), target, dimension( 0:2) :: c1_FD2 = (/-1.5_mytype, 2.0_mytype,-0.5_mytype/)
  real(mytype), target, dimension(-2:0) :: c1_BD2 = (/ 0.5_mytype,-2.0_mytype, 1.5_mytype/)
  !$acc declare create(c1_CD2,c1_CD4,c1_CD6)
  !$acc declare create(c2_CD2,c2_CD4,c2_CD6)
  !$acc declare create(c1_FD2,c1_BD2)
  ! define allocation coefficients for convective and diffusive derivatives
  real(mytype), allocatable, dimension(:  ) :: conv_ddx, conv_ddy, conv_ddz
  real(mytype), allocatable, dimension(:  ) :: visc_ddx,   visc_ddy,   visc_ddz
  real(mytype), allocatable, dimension(:,:) :: visc_d2dx2, visc_d2dz2
  real(mytype), allocatable, dimension(: )  :: visc_d2dy2
contains
! initialization finite difference coefficients
  subroutine init_derivCoeffs()
    integer :: i, k, c
    real(mytype), pointer, dimension(:) :: c1_conv, c1_visc, c2_visc
    allocate(conv_ddx(  nStencilConv))
    allocate(conv_ddy(  nStencilConv))
    allocate(conv_ddz(  nStencilConv))
    allocate(visc_ddx(  nStencilVisc))
    allocate(visc_ddy(  nStencilVisc))
    allocate(visc_ddz(  nStencilVisc))
    allocate(visc_d2dx2(-nStencilVisc:nStencilVisc, xsize(1)))
    allocate(visc_d2dy2(0:nStencilVisc))
    allocate(visc_d2dz2(-nStencilVisc:nStencilVisc, xsize(3)))
    ! selection of the convective order (1st derviative)
    if     (nStencilConv == 1) then;  c1_conv => c1_CD2; 
    elseif (nStencilConv == 2) then;  c1_conv => c1_CD4; 
    elseif (nStencilConv == 3) then;  c1_conv => c1_CD6;
    endif 
    ! selection of the diffusive order (1st derivative)
    if     (nStencilVisc == 1) then;  c1_visc => c1_CD2; 
    elseif (nStencilVisc == 2) then;  c1_visc => c1_CD4; 
    elseif (nStencilVisc == 3) then;  c1_visc => c1_CD6;
    endif 
    ! selection of the diffusive order (2nd derivative)
    if     (nStencilVisc == 1) then;  c2_visc => c2_CD2; 
    elseif (nStencilVisc == 2) then;  c2_visc => c2_CD4; 
    elseif (nStencilVisc == 3) then;  c2_visc => c2_CD6; 
    endif 
    ! calculation of 1st order
    conv_ddx = c1_conv/dx
    conv_ddy = c1_conv/dy
    conv_ddz = c1_conv/dz
    visc_ddx = c1_visc/dx
    visc_ddy = c1_visc/dy
    visc_ddz = c1_visc/dz
    ! calculation of 2nd order (with metrics): x- and z-direction stretched
    ! x-direction (wall-normal)
    visc_d2dx2 = 0.0_mytype
    do i=1,xsize(1)
      visc_d2dx2(0,i) = c2_visc(0)*xp(i)**2/dx**2  !! c1_visc(0)=0
      do c = 1, nStencilVisc
        visc_d2dx2( c,i) = c2_visc(c)*xp(i)**2/dx**2 - c1_visc(c)/dx*xpp(i)*xp(i)**3 
        visc_d2dx2(-c,i) = c2_visc(c)*xp(i)**2/dx**2 + c1_visc(c)/dx*xpp(i)*xp(i)**3
      enddo
    enddo
    ! y-direction (spanwise)
    visc_d2dy2 = c2_visc/dy**2
    ! z-direction (streamwise)
    visc_d2dz2 = 0.0_mytype
    do k=1,xsize(3)
      visc_d2dz2(0,k) = c2_visc(0)*zp(k)**2/dz**2  !! c1_visc(0)=0
      do c = 1, nStencilVisc
        visc_d2dz2( c,k) = c2_visc(c)*zp(k)**2/dz**2 - c1_visc(c)/dz*zpp(k)*zp(k)**3 
        visc_d2dz2(-c,k) = c2_visc(c)*zp(k)**2/dz**2 + c1_visc(c)/dz*zpp(k)*zp(k)**3
      enddo
    enddo
    !$acc enter data copyin(conv_ddx,conv_ddy,conv_ddz) async
    !$acc enter data copyin(visc_ddx,visc_ddy,visc_ddz) async
    !$acc enter data copyin(visc_d2dx2,visc_d2dy2,visc_d2dz2) async
  end subroutine
! forward convective difference at x-boundary
  subroutine calc_conv_FD_ddx(d1,arr,i,j,k,nH,xp,dx)
    !$acc routine seq
    implicit none
    real(mytype), dimension(1-nH:,1-nH:,1-nH:) :: arr
    real(mytype) :: d1, xp, dx
    integer :: i,j,k,c,nH
    d1 = 0.0_mytype
    do c = 0,2
      d1 = d1 + c1_FD2(c)*arr(i+c,j,k)
    enddo
    d1 = d1*xp/dx
  end subroutine
! backward convective difference at x-boundary
  subroutine calc_conv_BD_ddx(d1,arr,i,j,k,nH,xp,dx)
    !$acc routine seq
    implicit none
    real(mytype), dimension(1-nH:,1-nH:,1-nH:)  :: arr
    real(mytype) :: d1, xp, dx
    integer :: i,j,k,c,nH
    d1 = 0.0_mytype
    do c = -2, 0
      d1 = d1 + c1_BD2(c)*arr(i+c,j,k)
    enddo
    d1 = d1*xp/dx
  end subroutine
! Forward convective difference at z-boundary
  subroutine calc_conv_FD_ddz(d1,arr,i,j,k,nH,zp,dz)
    !$acc routine seq
    implicit none
    real(mytype), dimension(1-nH:,1-nH:,1-nH:)  :: arr
    real(mytype) :: d1, zp, dz
    integer :: i,j,k,c,nH
    d1 = 0.0_mytype
    do c = 0, 2
      d1 = d1 + c1_FD2(c)*arr(i,j,k+c)
    enddo
    d1 = d1*zp/dz
  end subroutine
! backward convective difference at z-boundary
  subroutine calc_conv_BD_ddz(d1,arr,i,j,k,nH,zp,dz)
    !$acc routine seq
    implicit none
    real(mytype), dimension(1-nH:,1-nH:,1-nH:) :: arr
    real(mytype) :: d1, zp, dz
    integer :: i,j,k,c,nH
    d1 = 0.0_mytype
    do c = -2, 0
      d1 = d1 + c1_BD2(c)*arr(i,j,k+c)
    enddo
    d1 = d1*zp/dz
  end subroutine
! central diffusive difference at x-boundary
  subroutine calc_visc_ddx(d1,arr,i,j,k)
    use decomp_2d
    use mod_param
    implicit none
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: arr
    real(mytype) :: d1
    integer :: i,j,k,c
    d1 = 0.0_mytype
    do c = 1, nStencilVisc
      d1 = d1 + visc_ddx(c)*(arr(i+c,j,k) - arr(i-c,j,k))
    enddo
    d1 = d1*xp(i)
  end subroutine
! central diffusive difference at y-boundary
  subroutine calc_visc_ddy(d1,arr,i,j,k)
    use decomp_2d
    use mod_param
    implicit none
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: arr
    real(mytype) :: d1
    integer :: i,j,k,c
    d1 = 0.0_mytype
    if (xsize(2) == 1) return 
    do c = 1, nStencilVisc
      d1 = d1 + visc_ddy(c)*(arr(i,j+c,k) - arr(i,j-c,k))
    enddo
  end subroutine
! central diffusive difference at z-boundary
  subroutine calc_visc_ddz(d1,arr,i,j,k)
    use decomp_2d
    use mod_param
    implicit none
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: arr
    real(mytype) :: d1
    integer :: i,j,k,c
    d1 = 0.0_mytype
    do c = 1, nStencilVisc
      d1 = d1 + visc_ddz(c)*(arr(i,j,k+c) - arr(i,j,k-c))
    enddo
    d1 = d1*zp(k)
  end subroutine
! central diffusive differences with loop (1st order)
  subroutine calc_visc_ddxyz(ddx,ddy,ddz,arr,i,j,k)
    use decomp_2d
    use mod_param
    implicit none
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: arr
    real(mytype) :: ddx,ddy,ddz
    integer :: i,j,k,c
    ddx = 0.0_mytype
    ddy = 0.0_mytype
    ddz = 0.0_mytype
    do c = 1, nStencilVisc
      ddx = ddx + visc_ddx(c)*(arr(i+c,j,k) - arr(i-c,j,k))
      ddy = ddy + visc_ddy(c)*(arr(i,j+c,k) - arr(i,j-c,k))
      ddz = ddz + visc_ddz(c)*(arr(i,j,k+c) - arr(i,j,k-c))
    enddo
    ddx = ddx*xp(i)
    ddz = ddz*zp(k)
  end subroutine
! central diffusive differences with loop (2nd order)
  subroutine calc_d2dxyz2(d2dx2,d2dy2,d2dz2,arr,i,j,k)
    use decomp_2d
    use mod_param
    implicit none
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: arr
    real(mytype) :: d2dx2,d2dy2,d2dz2
    integer :: i,j,k,c
    d2dx2 = visc_d2dx2(0,i)*arr(i,j,k)
    d2dy2 = visc_d2dy2(0)  *arr(i,j,k)
    d2dz2 = visc_d2dz2(0,k)*arr(i,j,k)
    do c = 1, nStencilVisc
      d2dx2 = d2dx2 + visc_d2dx2(c,i)*  arr(i+c,j,k) + visc_d2dx2(-c,i)*arr(i-c,j,k)
      d2dy2 = d2dy2 + visc_d2dy2(c)  * (arr(i,j+c,k) +                  arr(i,j-c,k) )
      d2dz2 = d2dz2 + visc_d2dz2(c,k)*  arr(i,j,k+c) + visc_d2dz2(-c,k)*arr(i,j,k-c) 
    enddo
  end subroutine
! print order of finite differences 
  subroutine print_init_derivCoeffs()
  use decomp_2d
  use mod_param
  implicit none
  if (nrank == 0) then
    write(stdout,* ) 'Finite differences'
    write(stdout,'(A)') 'o--------------------------------------------------o'
    write(stdout,'(A, I10)') 'Order convective fluxes:              ',2*nStencilConv 
    write(stdout,'(A, I10)') 'Order diffusive fluxes:               ',2*nStencilVisc
    write(stdout,* )
    else
  endif        
  end subroutine
end module
