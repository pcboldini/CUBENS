!----------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                                  !
!      mod_init                                                                                                                    !
!                                                                                                                                  !
!      DESCRIPTION:                                                                                                                !
!      ------------                                                                                                                !
!      This module is for the initialisation of the flow variables inside the domain                                               !
!                                                                                                                                  !
!                                                                                                                                  !
!      CHANGELOG:                                                                                                                  !
!      ----------                                                                                                                  !
!         xx.xx.2021: module created (Rene)                                                                                        !
!         03.12.2021: namelist initBlasius deleted (Pietro)                                                                        !
!         07.12.2021: removed initChannelFromRestartFile                                                                           !
!                                                                                                                                  !
!----------------------------------------------------------------------------------------------------------------------------------!

module mod_init

use io_std_units

implicit none

integer, parameter :: typeReadBlasius = 8

contains

subroutine initBlasius1D(part,xcoord,zcoord,rho,u,v,w,ien,pre,tem,mu,ka)

  use decomp_2d
  use mod_param
  use mod_eos
  use mod_eos_var
  ! use mod_table
  use mod_grid
  use mod_halo
  use mod_math

  implicit none

  integer :: i,j,k, npts, fs, ierr, kk, lenr
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka
  
  real(8),      allocatable, dimension(:) :: xRead, wRead, uRead, rRead
  real(mytype), allocatable, dimension(:) ::        wIntp, uIntp, rIntp
  real(mytype), dimension(:) :: xcoord, zcoord
  real(mytype) :: dummy, scaling_delta, x0, fact, integral !!!
  TYPE (DECOMP_INFO), intent(IN) :: part
  real(mytype), dimension(1:xsize(1)) :: vector, siny
  real(mytype) :: factz, ztem_1, ztem_2

  inquire(FILE='initBL/inputDNS/prof_x.bin', SIZE=fs); npts = fs/8

  if (nrank == 0) write(stdout,*) "init profiles with npts = ", npts

  allocate(xRead(npts))
  allocate(rRead(npts)); allocate(rIntp(npts))
  allocate(wRead(npts)); allocate(wIntp(npts))
  allocate(uRead(npts)); allocate(uIntp(npts))


  open(9,file='initBL/inputDNS/prof_x.bin',access='direct',recl=fs); read(9,rec=1) (xRead(i),i=1,npts);close(9)
  open(9,file='initBL/inputDNS/prof_r.bin',access='direct',recl=fs); read(9,rec=1) (rRead(i),i=1,npts);close(9)
  open(9,file='initBL/inputDNS/prof_w.bin',access='direct',recl=fs); read(9,rec=1) (wRead(i),i=1,npts);close(9)
  open(9,file='initBL/inputDNS/prof_u.bin',access='direct',recl=fs); read(9,rec=1) (uRead(i),i=1,npts);close(9)

  call spline(xRead, rRead, npts, rIntp)
  call spline(xRead, wRead, npts, wIntp)
  call spline(xRead, uRead, npts, uIntp)


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

  v = 0.0_mytype

  select case (t_param%USE_EOS)
  case("IG")
    pre = Pref*t_ig%prefac_r  
  case("VdW")
    pre = Pref*t_vdw%prefac_r 
  case("PR")
    pre = Pref*t_pr%prefac_r 
  end select

  !$acc update device(rho,pre,u,v,w)

  call calcState_rP(rho,pre,ien,tem,mu,ka, 1,part%xsz(1),1,part%xsz(2),1,part%xsz(3)) !


  ! if (BC_bot == "temp_std") then

  !   ztem_1 = ReTem1**2/Re/zEndDNS ! zStart_tem
  !   ztem_2 = ReTem2**2/Re/zEndDNS ! zEnd_tem

  !   do k=1,part%xsz(3)
  !     kk = k + xstart(3) - 1
  !     factz = (kk - 1.0_mytype)/(nz_global-1.0_mytype)
  !     do j=1,part%xsz(2)
  !       tem(1,j,k) = Twall_new+0.5_mytype*(Twall-Twall_new)*&
  !                      (2.0_mytype-tanh((factz-ztem_1)/(deltatem1))+tanh((factz-ztem_2)/(deltatem2)))
  !     enddo
  !   enddo

  !   call calcState_pT(pre,tem,rho,ien,mu,ka, 1,part%xsz(1),1,part%xsz(2),1,part%xsz(3)) !

  ! endif

  deallocate(xRead)
  deallocate(rRead); deallocate(rIntp)
  deallocate(wRead); deallocate(wIntp)
  deallocate(uRead); deallocate(uIntp)


end subroutine initBlasius1D


subroutine initField_1D(rho,u,v,w,ien)

  use decomp_2d
  use mod_param
  use mod_eos
  use mod_grid
  use mod_math

  implicit none

  integer :: i,j,k,jj,kk
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien
  real(mytype), dimension(1-nHalo:xsize(1)+nHalo,1-nHalo:xsize(2)+nHalo,1-nHalo:xsize(3)+nHalo) :: pres
  real(mytype) :: fz,P0,T0

  u   = 0.0_mytype
  v   = 0.0_mytype
  w   = 1.0_mytype
  T0  = 1.0_mytype

  do i=1,xsize(1)
    do j=1,xsize(2)
      do k=1,xsize(3)

        kk = k + xstart(3) - 1
        fz = len_z*(kk-1.0_mytype)/(nz_global-1)
        rho(i,j,k) = 2.0_mytype + 1.0_mytype*(sin(2.0_mytype*pi_const*fz))
        P0 = rho(i,j,k)*eos_Rgas*T0
        ien(i,j,k) = P0/(ig_gam-1.0)

      enddo
    enddo
  enddo

end subroutine



subroutine initField(rho,u,v,w,ien)

  use decomp_2d
  use mod_param
  use mod_eos
  use mod_grid

  implicit none

  integer :: i,j,k,jj,kk
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien
  real(mytype), dimension(1-nHalo:xsize(1)+nHalo,1-nHalo:xsize(2)+nHalo,1-nHalo:xsize(3)+nHalo) :: pres
  real(mytype) :: pInit,temp,rhoInit,tempInit,wb ! deltax
  real(mytype) :: rho_0,visc_0,cp_0,sos_0,inte_0,U_0
  real(mytype) :: inte,dummy, fx,fy,fz, tmp, V0
  real(mytype) :: rnum

  rho = 1.0_mytype
  u   = 0.0_mytype
  v   = 0.0_mytype
  ! w   = 1.0

! set BC for velocity 
  ! w(1,:,:) = 0.0
  ! w(xsize(1),:,:) = 0.0

  call random_number(rnum)

  do i=1,xsize(1)
    do j=1-nHalo,xsize(2)+nHalo

      jj = j + xstart(2) - 1
      fy = 1.0_mytype*(jj)/(ny_global)

      do k=1,xsize(3)

        kk = k + xstart(3) - 1
        fz = len_z*(kk-1.0_mytype)/(nz_global-1)

        w(i,j,k) = 2.0_mytype/2.0_mytype*x(i)*(2.0_mytype-x(i))

!        if (x(i) < 1.0) then
!          yplus = x(i)*450.0
!          if  (yplus .lt. 11.6)  w(i,j,k) = yplus
!          if  (yplus .gt. 11.6)  w(i,j,k) = (2.5*log(yplus)+5.5)
!        else
!          yplus = (2.0-x(i))*450.0
!          if  (yplus .lt. 11.6)  w(i,j,k) = yplus
!          if  (yplus .gt. 11.6)  w(i,j,k) = (2.5*log(yplus)+5.5)
!        endif
      enddo
    enddo
  enddo



  pres = Pref !!!! only IG
  ien = pres/rho/(ig_gam-1.0_mytype)


end subroutine initField


subroutine initField_TGV(part,rho,u,v,w,ien,pre,tem,mu,ka)

  use decomp_2d
  use mod_param
  use mod_eos
  use mod_grid

  implicit none

  integer :: i,j,k,jj,kk
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka
  real(mytype) :: fx,fy,fz, V0,P0,R0,T0,press
  TYPE (DECOMP_INFO), intent(IN) :: part


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

call calcState_rP(rho,pre,ien,tem,mu,ka, 1,part%xsz(1),1,part%xsz(2),1,part%xsz(3)) !

end subroutine

end module mod_init




