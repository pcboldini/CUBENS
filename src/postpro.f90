
module mod_postpro

  use decomp_2d
  use mod_param
  use mod_grid
  use mod_finitediff
  use mod_eos
  use mod_eos_var

  implicit none

!   interface
!    subroutine assemble_globalz(local1,global1,nx,ny,nz)
!      use decomp_2d
!      integer, intent(IN) :: nx,ny,nz
!      real(mytype), dimension(:,:,:), intent(IN) :: local1 !!!!
!      real(mytype), dimension(nx,ny,nz_global), intent(OUT) :: global1 !!!!
!    end subroutine assemble_globalz
! end interface


contains


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
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','postproc/planes/'//trim(name011)//trim(cha2)//'.'//trim(name012)//cha//'.bin', &
                                 'dummy',part)
    endif

  end subroutine


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

          cp_r    = t_vdw%cvOverR+1.0/(1.0-((3.0*v_r-1.0)**2/(4.0*tem_r*v_r**3.0)))
          Cp(i,j,k) = cp_r*t_vdw%Cpfac_r

        enddo
      enddo
    enddo

  end subroutine

#endif


 subroutine spectraz(pen,speca,nfiles,part1,partf,a,b)                                                                                                                                                                

  use decomp_2d
  use decomp_2d_fft
  use decomp_2d_io
  use mod_param

  implicit none

  integer :: nfiles,i,j,k,pen
  real(mytype), dimension(:,:,:) :: speca 
  real(mytype), dimension(:,:,:), optional :: a,b 
  TYPE(DECOMP_INFO), intent(IN) :: part1,partf

  real(mytype), allocatable, dimension(:,:,:)    :: ayt,azt,byt,bzt,specrxt,specryt,specrzt 
  complex(mytype), allocatable, dimension(:,:,:) :: aspec, bspec, spec
!  TYPE(DECOMP_INFO) :: partf

!  call decomp_2d_fft_init(PHYSICAL_IN_Z) !PHYSICAL IN Z IS IMPORTANT ONLY WHEN WE ARE DOING 3D FFT
  
!  call decomp_info_init(imax,jmax,kmax/2+1,partf)
  if (nrank==0) write(*,*) 'xpencil', partf%xsz(1), partf%xsz(2), partf%xsz(3)
  if (nrank==0) write(*,*) 'ypencil', partf%ysz(1), partf%ysz(2), partf%ysz(3) 
  if (nrank==0) write(*,*) 'zpencil', partf%zsz(1), partf%zsz(2), partf%zsz(3) 
  
  allocate(  ayt(1:part1%ysz(1), 1:part1%ysz(2), 1:part1%ysz(3))   )
  allocate(  byt(1:part1%ysz(1), 1:part1%ysz(2), 1:part1%ysz(3))   )
  allocate(  azt(1:part1%zsz(1), 1:part1%zsz(2), 1:part1%zsz(3))   )
  allocate(  bzt(1:part1%zsz(1), 1:part1%zsz(2), 1:part1%zsz(3))   )

  allocate(specrxt(1:partf%xsz(1), 1:partf%xsz(2), 1:partf%xsz(3))   )
  allocate(specryt(1:partf%ysz(1), 1:partf%ysz(2), 1:partf%ysz(3))   )
  allocate(specrzt(1:partf%zsz(1), 1:partf%zsz(2), 1:partf%zsz(3))   )

  allocate( aspec(1:partf%zsz(1), 1:partf%zsz(2), 1:partf%zsz(3))   )
  allocate( bspec(1:partf%zsz(1), 1:partf%zsz(2), 1:partf%zsz(3))   )
  allocate( spec(1:partf%zsz(1), 1:partf%zsz(2), 1:partf%zsz(3))   )

 ! write(*,*) 'a=', a(:,1,1)

  !call transpose_x_to_y(a,ayt) 
  !call transpose_y_to_z(ayt,azt) 

  

!  write(*,*) 'a=', a(1,1,2)
!  write(*,*) 'ayt=', ayt(1,1,2)

!  write(*,*) 'azt=', azt(,1,1)
 
  call decomp_2d_fft_1d(a,aspec,2)

 ! write(*,*) 'aspec=', aspec(2,:,100)

  stop


  !if (nrank==0) write(*,*) 'Imag part of spec is', aimag(aspec(1,1,2))

   if(present(b)) then 
      write(*,*) 'b present'
       call transpose_x_to_y(b, byt) 
       call transpose_y_to_z(byt, bzt) 

       call decomp_2d_fft_1d(bzt,bspec,2)

       bspec = conjg(bspec) ! compute complex conjugate of bspec
       spec = aspec*bspec
   else
       spec = aspec
   endif



   !if (nrank==0) write(*,*) 'Imag part of spec is', aimag(aspec(1,1,1)), aimag(bspec(200,1,100))

   !specrzt = abs(spec)

  

   call transpose_z_to_y(specrzt,specryt,partf)

!      do k=1,partf%ysz(3)
!        do i=1,partf%ysz(1)
!            specryt(i,1,k) = sum(specryt(i,:,k))/ny_global
!        enddo
!      enddo
        do j=1,partf%ysz(2)
            specryt(:,j,:) = specryt(:,1,:) 
        enddo

   call transpose_y_to_x(specryt,specrxt,partf)   !output in the form of x pencil
   call transpose_y_to_z(specryt,specrzt,partf)   !output in the form of z pencil

   ! if (pen == 1) then
   !       speca(:,:) = speca(:,:) + specrxt(:,1,:)/nfiles
   ! elseif (pen == 2) then
   !       speca(:,:) = speca(:,:) + specryt(:,1,:)/nfiles
   ! elseif (pen == 3) then
   !       speca(:,:) = speca(:,:) + specrzt(:,1,:)/nfiles
   ! endif

   deallocate(ayt)
   deallocate(byt)
   deallocate(azt)
   deallocate(bzt)

  deallocate(specrxt)
  deallocate(specryt)
  deallocate(specrzt)

  deallocate( aspec)
  deallocate( bspec)
  deallocate( spec)

  !call decomp_2d_fft_finalize

  end subroutine spectraz 


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

  real(mytype), allocatable, dimension(:,:,:)    :: ayt,azt,byt,bzt 
  complex(mytype), allocatable, dimension(:,:,:) :: aspec, bspec, spec,specrxt,specryt,specrzt
!  TYPE(DECOMP_INFO) :: partf

  !call decomp_2d_fft_init() !PHYSICAL IN Z IS IMPORTANT ONLY WHEN WE ARE DOING 3D FFT

!  if (nrank==0) write(*,*) 'xpencil', partf%xsz(1), partf%xsz(2), partf%xsz(3)
!  if (nrank==0) write(*,*) 'ypencil', partf%ysz(1), partf%ysz(2), partf%ysz(3) 
!  if (nrank==0) write(*,*) 'zpencil', partf%zsz(1), partf%zsz(2), partf%zsz(3) 
  
  allocate(  ayt(1:part1%ysz(1), 1:part1%ysz(2), 1:part1%ysz(3))   )
  allocate(  byt(1:part1%ysz(1), 1:part1%ysz(2), 1:part1%ysz(3))   )
  allocate(  azt(1:part1%zsz(1), 1:part1%zsz(2), 1:part1%zsz(3))   )
  allocate(  bzt(1:part1%zsz(1), 1:part1%zsz(2), 1:part1%zsz(3))   )

  allocate(specrxt(1:partf%xsz(1), 1:partf%xsz(2), 1:partf%xsz(3))   )
  allocate(specryt(1:partf%ysz(1), 1:partf%ysz(2), 1:partf%ysz(3))   )
  allocate(specrzt(1:partf%zsz(1), 1:partf%zsz(2), 1:partf%zsz(3))   )

  allocate( aspec(1:partf%ysz(1), 1:partf%ysz(2), 1:partf%ysz(3))   )
  allocate( bspec(1:partf%ysz(1), 1:partf%ysz(2), 1:partf%ysz(3))   )
  allocate(  spec(1:partf%ysz(1), 1:partf%ysz(2), 1:partf%ysz(3))   )

  call transpose_x_to_y(a, ayt) 

  call decomp_2d_fft_1d(ayt,aspec,2)

  !write(*,*) 'aspec(134,:,167)=', aspec(134,:,167)

   if(present(b)) then 
       call transpose_x_to_y(b, byt) 

       call decomp_2d_fft_1d(byt,bspec,2)

       bspec = conjg(bspec) ! compute complex conjugate of bspec
       specryt = aspec*bspec
   else
       specryt = aspec
   endif

!  specryt = abs(spec)

  call transpose_y_to_z(specryt,specrzt,partf)

      ! do j=1,partf%zsz(2)
      !   do i=1,partf%zsz(1)
      !       specrzt(i,j,1) = sum(specrzt(i,j,:))/nz_global
      !   enddo
      ! enddo
      !   do k=1,partf%zsz(3)
      !       specrzt(:,:,k) = specrzt(:,:,1) 
      !   enddo

   !call transpose_z_to_y(specrzt,specryt,partf)   !output in the form of x pencil
   call transpose_y_to_x(specryt,specrxt,partf)   !output in the form of z pencil

   ! if (pen == 1) then
   !       speca(:,:,:) = speca(:,:,:) + specrxt(:,:,1)/nfiles
   ! elseif (pen == 2) then
   !       speca(:,:,:) = speca(:,:,:) + specryt(:,:,1)/nfiles
   ! elseif (pen == 3) then
   !       speca(:,:,:) = speca(:,:,:) + specrzt(:,:,1)/nfiles
   ! endif


   speca = specryt/(nfiles-1)
  ! write(*,*) 'speca(134,:,167)=', speca(134,:,167)

   deallocate(ayt)
   deallocate(byt)
   deallocate(azt)
   deallocate(bzt)

  deallocate(specrxt)
  deallocate(specryt)
  deallocate(specrzt)

  deallocate( aspec)
  deallocate( bspec)
  deallocate( spec)

  end subroutine spectray

  ! subroutine avgField(r,u,v,w)

  !   use decomp_2d
  !   implicit none
  !   integer i,j,k
  !   real(mytype), dimension(:,:,:) :: lambda
  !   real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: u,v,w
  !   real(mytype) :: dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz
  !   real(mytype) :: dil, sxx, sxy, sxz, syy, syz, szz
  !   real(mytype) ::      oxx, oxy, oxz, oyy, oyz, ozz

  !   allocate(tau(xsize(1),xsize(2),xsize(3)))



  !   do k=1,xsize(3)
  !     do j=1,xsize(2)
  !       do i=1,xsize(1)
  !         call ddxyz(dux, duy, duz, u, i,j,k)
  !         call ddxyz(dvx, dvy, dvz, v, i,j,k)
  !         call ddxyz(dwx, dwy, dwz, w, i,j,k)

  !         ! dil = dux + dvy + dwz
  !         sxx =  dux
  !         sxy = (duy + dvx)/2.0
  !         sxz = (duz + dwx)/2.0
  !         syy =  dvy
  !         syz = (dvz + dwy)/2.0
  !         szz =  dwz

  !         oxx =  0.0
  !         oxy = (duy - dvx)/2.0
  !         oxz = (duz - dwx)/2.0
  !         oyy =  0.0
  !         oyz = (dvz - dwy)/2.0
  !         ozz =  0.0
  
  !         lambda(i,j,k)= 0.5*(  oxx**2.0 + oxy**2.0 + oxz**2.0 &
  !                             + oxy**2.0 + oyy**2.0 + oyz**2.0 &
  !                             + oxz**2.0 + oyz**2.0 + ozz**2.0 & 
  !                             - sxx**2.0 - sxy**2.0 - sxz**2.0 &
  !                             - sxy**2.0 - syy**2.0 - syz**2.0 &
  !                             - sxz**2.0 - syz**2.0 - szz**2.0 )
  !       enddo
  !     enddo
  !   enddo

  ! end subroutine

end module

! subroutine assemble_globalz(local1,global1,nx,ny,nz)
  
!   use decomp_2d
!   use mpi
  
!   implicit none
  
!   integer, intent(IN) :: nx,ny,nz
!   real(mytype), dimension(:,:,:), intent(IN) :: local1 !!!!
!   real(mytype), dimension(nx,ny,nz_global), intent(OUT) :: global1 !!!!
  
!   real(mytype), allocatable, dimension(:,:,:) :: rbuf_1 !!!!!
!   integer, dimension(9) :: sbuf1, rbuf1
  
!   integer :: ierror, i,j,k,m, i1,i2,j1,j2,k1,k2, count
!   integer, dimension(MPI_STATUS_SIZE) :: status
  
!   if (nrank==0) then
!      ! master writes its own data to a global array
!         i1 = 1
!         i2 = nx_global
!         j1 = 1
!         j2 = ny
!         k1 = 1
!         k2 = nz
!      do k=k1,k2
!         do j=j1,j2
!            do i=i1,i2
!               ! 'local' is assumbed shape array
!               ! but it is OK as starting index for rank 0 always 1
!               global1(i,j,k)=local1(i,j,k)
!            end do
!         end do
!      end do
!      ! then loop through all other ranks to collect data
!      do m=1,nproc-1
!         CALL MPI_RECV(rbuf1,9,MPI_INTEGER,m,m,MPI_COMM_WORLD, &
!              status,ierror)
!         allocate(rbuf_1(rbuf1(1):rbuf1(2),rbuf1(4):rbuf1(5), &
!              rbuf1(7):rbuf1(8)))
!         CALL MPI_RECV(rbuf_1,rbuf1(3)*rbuf1(6)*rbuf1(9),real_type,m, &
!              m+nproc,MPI_COMM_WORLD,status,ierror) !!!!!
!         do k=rbuf1(7),rbuf1(8)
!            do j=rbuf1(4),rbuf1(5)
!               do i=rbuf1(1),rbuf1(2)
!                  global1(i,j,k)=rbuf_1(i,j,k)
!               end do
!            end do
!         end do
!         deallocate(rbuf_1)
!      end do
!   else
!      ! slaves send data to mater
!         sbuf1(1) = 1
!         sbuf1(2) = nx_global
!         sbuf1(3) = nx_global
!         sbuf1(4) = 1
!         sbuf1(5) = ny
!         sbuf1(6) = ny
!         sbuf1(7) = nrank*nz+1
!         sbuf1(8) = (nrank+1)*nz
!         sbuf1(9) = (nrank+1)*nz
!         count = nx*ny*nz
!      ! send partition information
!      CALL MPI_SEND(sbuf1,9,MPI_INTEGER,0,nrank,MPI_COMM_WORLD,ierror)
!      ! send data array
!      CALL MPI_SEND(local1,count,real_type,0, &
!           nrank+nproc,MPI_COMM_WORLD,ierror) !!!!!
!   end if
  
!   return
! end subroutine assemble_globalz






