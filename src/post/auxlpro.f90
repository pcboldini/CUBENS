! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -

module mod_auxlpro
  use decomp_2d
  use mod_param
  use mod_grid
  use mod_finitediff
  use mod_eos
  use mod_eos_var
  implicit none
contains


! calculate Q-criterion
  subroutine calcQ(lambda,u,v,w,part)
    use mpi
    use decomp_2d
    use decomp_2d_io
    implicit none
    integer i,j,k,c
    character*7 cha
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    TYPE(DECOMP_INFO), intent(IN) :: part
    real(mytype), dimension(:,:,:) :: lambda
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: u,v,w
    real(mytype) :: dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz
    real(mytype) :: sxx, sxy, sxz, syy, syz, szz
    real(mytype) ::      oxx, oxy, oxz, oyy, oyz, ozz

    do k=1,part%xsz(3)
      do j=1,part%xsz(2)
        do i=1,part%xsz(1)
          dux   = 0.0_mytype
          duy   = 0.0_mytype
          duz   = 0.0_mytype
          dvx   = 0.0_mytype
          dvy   = 0.0_mytype
          dvz   = 0.0_mytype
          dwx   = 0.0_mytype
          dwy   = 0.0_mytype
          dwz   = 0.0_mytype
          do c = 1, nStencilVisc
            dux   = dux   + visc_ddx(c)*( u(i+c,j,k) -  u(i-c,j,k))*xp(i)
            duy   = duy   + visc_ddy(c)*( u(i,j+c,k) -  u(i,j-c,k))
            duz   = duz   + visc_ddz(c)*( u(i,j,k+c) -  u(i,j,k-c))*zp(k)
            dvx   = dvx   + visc_ddx(c)*( v(i+c,j,k) -  v(i-c,j,k))*xp(i)
            dvy   = dvy   + visc_ddy(c)*( v(i,j+c,k) -  v(i,j-c,k))
            dvz   = dvz   + visc_ddz(c)*( v(i,j,k+c) -  v(i,j,k-c))*zp(k)
            dwx   = dwx   + visc_ddx(c)*( w(i+c,j,k) -  w(i-c,j,k))*xp(i)
            dwy   = dwy   + visc_ddy(c)*( w(i,j+c,k) -  w(i,j-c,k))
            dwz   = dwz   + visc_ddz(c)*( w(i,j,k+c) -  w(i,j,k-c))*zp(k)
          enddo
          sxx =  dux
          sxy = (duy + dvx)/2.0_mytype
          sxz = (duz + dwx)/2.0_mytype
          syy =  dvy
          syz = (dvz + dwy)/2.0_mytype
          szz =  dwz
          oxx =  0.0_mytype
          oxy = (duy - dvx)/2.0_mytype
          oxz = (duz - dwx)/2.0_mytype
          oyy =  0.0_mytype
          oyz = (dvz - dwy)/2.0_mytype
          ozz =  0.0_mytype
          lambda(i,j,k)= 0.5_mytype*(  oxx**2 + oxy**2 + oxz**2 &
                              + oxy**2 + oyy**2 + oyz**2 &
                              + oxz**2 + oyz**2 + ozz**2 & 
                              - sxx**2 - sxy**2 - sxz**2 &
                              - sxy**2 - syy**2 - syz**2 &
                              - sxz**2 - syz**2 - szz**2 )
        enddo
      enddo
    enddo
  end subroutine


! Calculate gradients
  subroutine calcGrad(part,istep,dir,loc,tmp01,name011,name012,tmp02,tmp03,tmp04,tmp05,tmp06,tmp07,tmp08,tmp09,tmp10,tmp11,tmp12)
    use mpi
    use decomp_2d
    use decomp_2d_io
    implicit none
    integer i,j,k,c,loc,dir,istep,ierr,comm
    character*7 cha
    character(len=1024) :: cha2
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    TYPE(DECOMP_INFO), intent(IN) :: part
    real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:), optional :: &
                        tmp01,tmp02,tmp03,tmp04,tmp05,tmp06,tmp07,tmp08,tmp09,tmp10,tmp11,tmp12
    character(len=*), optional :: name011 !,name021,name031,name041,name051,name061,name071,name081,name091,name101,name111,name121
    character(len=*), optional :: name012 !,name022,name032,name042,name052,name062,name072,name082,name092,name102,name112,name122
    real(mytype), dimension(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3)) :: &
                        Gtmp01,Gtmp02,Gtmp03,Gtmp04,Gtmp05,Gtmp06,Gtmp07,Gtmp08,Gtmp09,Gtmp10,Gtmp11,Gtmp12
    real(mytype) :: d01x, d01y, d01z     
    real(mytype), allocatable, dimension(:,:,:) :: tmp

    write(cha,'(I0.7)') istep  
    write(cha2,'(I0)') loc      
    do k=1,part%xsz(3)
      do j=1,part%xsz(2)
        do i=1,part%xsz(1)
          d01x   = 0.0_mytype
          d01y   = 0.0_mytype
          d01z   = 0.0_mytype
          do c = 1, nStencilVisc
            d01x   = d01x   + visc_ddx(c)*( tmp01(i+c,j,k) -  tmp01(i-c,j,k))*xp(i)
            d01y   = d01y   + visc_ddy(c)*( tmp01(i,j+c,k) -  tmp01(i,j-c,k))
            d01z   = d01z   + visc_ddz(c)*( tmp01(i,j,k+c) -  tmp01(i,j,k-c))*zp(k)
          enddo
          Gtmp01(i,j,k) = sqrt(  d01x**2 + d01y**2 + d01z**2 )
        enddo
      enddo
    enddo
    if(present(tmp01)) then 
      tmp = Gtmp01(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','output/planes/'//trim(name011)//trim(cha2)//'.'//trim(name012)//cha//'.bin', &
                                 'dummy',part)
    endif
  end subroutine


! Calculate stress tensor
  subroutine calcStrain(dilla2,sxx,sxy,sxz,syy,syz,szz,u,v,w) 
    use decomp_2d
    implicit none
    integer i,j,k
    real(mytype), dimension(:,:,:)   :: dilla2,sxx,sxy,sxz,syy,syz,szz
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: u,v,w
    real(mytype) :: dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz

    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          call calc_visc_ddxyz(dux, duy, duz, u, i,j,k)
          call calc_visc_ddxyz(dvx, dvy, dvz, v, i,j,k)
          call calc_visc_ddxyz(dwx, dwy, dwz, w, i,j,k)
          dilla2(i,j,k) = dux + dvy + dwz ! divergence
          sxx(i,j,k) = 2.0_mytype*dux - 2.0_mytype/3.0_mytype*dilla2(i,j,k) ! stress-tensor (2,2)
          sxy(i,j,k) =     duy + dvx ! stress-tensor (2,3)=(3,2)
          sxz(i,j,k) =     duz + dwx ! stress-tensor (2,1)=(1,2)
          syy(i,j,k) = 2.0_mytype*dvy - 2.0_mytype/3.0_mytype*dilla2(i,j,k) ! stress-tensor (3,3)
          syz(i,j,k) =     dvz + dwy ! stress-tensor (3,1)=(1,3)
          szz(i,j,k) = 2.0_mytype*dwz - 2.0_mytype/3.0_mytype*dilla2(i,j,k) ! stress-tensor (1,1)
        enddo
      enddo
    enddo
  end subroutine


! Calculate temperature gradient
  subroutine calcTemp(tmp_x_arr,tmp_y_arr,tmp_z_arr,tem) 
    use decomp_2d
    implicit none
    integer i,j,k
    real(mytype), dimension(:,:,:)   :: tmp_x_arr,tmp_y_arr,tmp_z_arr
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: tem
    real(mytype) :: dTx, dTy, dTz

    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          call calc_visc_ddxyz(dTx, dTy, dTz, tem, i,j,k)
          tmp_x_arr(i,j,k) = dTx
          tmp_y_arr(i,j,k) = dTy
          tmp_z_arr(i,j,k) = dTz
        enddo
      enddo
    enddo
  end subroutine


! Calculate specific heat at constant pressure
! Ideal gas
#ifdef IG
  subroutine calcCp(Cp,rho,ien)
    use decomp_2d
    implicit none
    integer i,j,k
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:), intent(IN) :: rho,ien
    real(mytype), dimension(:,:,:), intent(OUT) :: Cp

    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          Cp(i,j,k) = t_ig%cp
        enddo
      enddo
    enddo
  end subroutine
#endif

  
! Van der Waals
#ifdef VdW
  subroutine calcCp(Cp,rho,ien)
    use decomp_2d
    implicit none
    integer i,j,k
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:), intent(IN) :: rho,ien
    real(mytype), dimension(:,:,:), intent(OUT) :: Cp
    real(mytype) :: rho_r, v_r, tem_r, ien_r, cp_r

    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          rho_r   = t_param%Rhoref*rho(i,j,k)  
          v_r     = 1.0_mytype/rho_r
          ien_r   = ien(i,j,k)/t_vdw%Efac_r
          tem_r   = t_vdw%Zc/t_vdw%cvOverR*(ien_r + t_vdw%a*rho_r)
          cp_r    = t_vdw%cvOverR+1.0_mytype/(1.0_mytype-((3.0_mytype*v_r-1.0_mytype)**2/(4.0_mytype*tem_r*v_r**3)))
          Cp(i,j,k) = cp_r*t_vdw%Cpfac_r
        enddo
      enddo
    enddo
  end subroutine
#endif


! Peng-Robinson
#ifdef PR
  subroutine calcCp(Cp,rho,ien)
    use decomp_2d
    implicit none
    integer i,j,k
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:), intent(IN) :: rho,ien
    real(mytype), dimension(:,:,:), intent(OUT) :: Cp
    real(mytype) :: rho_r, v_r, tem_r, ien_r, cv_r, alpha, dPdrho_T_r, dPdT_rho_r, cp_r
    real(mytype) :: T1, F1, F2, F3
    real(mytype) :: cvOverR, t_pr_a, t_pr_b, t_pr_Zc1, t_pr_Zc2, t_pr_K, sqrt2
    ! Reduced quantities
    cvOverR = t_pr%cvOverR
    t_pr_a = t_pr%a
    t_pr_b = t_pr%b
    t_pr_Zc1 = t_pr%Zc1
    t_pr_Zc2 = t_pr%Zc2
    t_pr_K = t_pr%K
    sqrt2=2**(0.5_mytype)
    ! Non-dimensional quantities
    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          ien_r(i,j,k) = ien(i,j,k)/t_pr%Efac_r
          rho_r(i,j,k) = t_pr%Rhoref*rho(i,j,k) 
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
          Cp(i,j,k) = cp_r*t_pr%Cpfac_r
        enddo
      enddo
    enddo
  end subroutine
#endif


! FFT in spanwsie direction
  subroutine spectray(pen,speca,nfiles,part1,partf,a,b)
  use decomp_2d
  use decomp_2d_fft
  use decomp_2d_io
  use decomp_2d_constants
  use mod_param
  implicit none
  integer :: nfiles,i,j,k,pen
  complex(mytype), dimension(:,:,:) :: speca
  real(mytype), dimension(:,:,:), optional :: a,b 
  TYPE(DECOMP_INFO), intent(IN) :: part1,partf
  real(mytype), allocatable, dimension(:,:,:)    :: ayt,byt 
  complex(mytype), allocatable, dimension(:,:,:) :: aspec, bspec,specrxt,specryt,specrzt
  allocate(  ayt(1:part1%ysz(1), 1:part1%ysz(2), 1:part1%ysz(3))   )
  allocate(  byt(1:part1%ysz(1), 1:part1%ysz(2), 1:part1%ysz(3))   )
  allocate(specrxt(1:partf%xsz(1), 1:partf%xsz(2), 1:partf%xsz(3))   )
  allocate(specryt(1:partf%ysz(1), 1:partf%ysz(2), 1:partf%ysz(3))   )
  allocate(specrzt(1:partf%zsz(1), 1:partf%zsz(2), 1:partf%zsz(3))   )
  allocate( aspec(1:partf%ysz(1), 1:partf%ysz(2), 1:partf%ysz(3))   )
  allocate( bspec(1:partf%ysz(1), 1:partf%ysz(2), 1:partf%ysz(3))   )

  call transpose_x_to_y(a, ayt) 
  call decomp_2d_fft_1d(ayt,aspec,2)
   if(present(b)) then 
       call transpose_x_to_y(b, byt) 
       call decomp_2d_fft_1d(byt,bspec,2)
       bspec = conjg(bspec) ! compute complex conjugate of bspec
       specryt = aspec*bspec
   else
       specryt = aspec
   endif
  call transpose_y_to_z(specryt,specrzt,partf)
  call transpose_y_to_x(specryt,specrxt,partf)   !output in the form of z pencil
  speca = specryt/(nfiles-1)
  deallocate(ayt)
  deallocate(byt)
  deallocate(specrxt)
  deallocate(specryt)
  deallocate(specrzt)
  deallocate( aspec)
  deallocate( bspec)
  end subroutine
end module
