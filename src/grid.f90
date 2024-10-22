! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
! physical and numerical grid module

module mod_grid
  use decomp_2d
  use mod_param
  implicit none
! definition of the grid sizes, coordinate stencils, and their derivatives
  integer :: Kmult
  real(mytype) :: dx,dy,dz,dt_FFT
  real(mytype), allocatable, dimension(:) :: x,y,z,z_global
  real(mytype), allocatable, dimension(:) :: xp,xpp,zp,zpp,zp_global,zpp_global
  contains


! initialization of the grid variables and subroutines
  subroutine init_grid()
    use decomp_2d
    use mod_param
    use mod_halo
    implicit none
    ! definition of the variables
    integer      :: i,j,k,kk,c,lenr,ierr
    character*4 cha
    real(mytype), allocatable, dimension(:) :: tmp, tmp1, tmp2
    ! Computational mesh for non_equid: dx, dz; physical mesh elsewhere
    ! x-direction (wall-normal)
    if (xmesh_type == "equid") then
      ! x goes from 0 to len_x
      dx = len_x/(nx_global-1)  
      if (perBC(1) .eqv. .true.) dx = len_x/(nx_global)
    else if ( xmesh_type == "non_equid" ) then
      ! eta goes from 0 to 1
      dx = 1.0_mytype/(nx_global-1)
    endif
    ! y-direction (spanwise)
    dy = len_y/(ny_global)
    ! z-direction (streamwise)
    if (zmesh_type == "equid") then
      ! z goes from 0 to len_z
      dz = len_z/(nz_global-1) 
      if (perBC(3) .eqv. .true.) dz = len_z/(nz_global)
    else if ( zmesh_type == "non_equid" ) then
      ! xi goes from 0 to 1
      dz = 1.0_mytype/(nz_global-1) 
    endif
    ! allocation
    allocate(  x(xsize(1)))
    allocate( xp(xsize(1)))
    allocate(xpp(xsize(1)))
    allocate(  y(ny_global))
    allocate(  z(xsize(3)))
    allocate( zp(xsize(3)))
    allocate(zpp(xsize(3)))
    allocate(z_global(nz_global))
    allocate(zp_global(nz_global))
    allocate(zpp_global(nz_global))
    allocate(tmp(nz_global))
    allocate(tmp1(nz_global))
    allocate(tmp2(nz_global))
#if !defined(CHA)
      ! for boundary layer and TGV
      call xDistribution(xmesh_type,ReTau, gridStretchX, len_x, xsize(1), x, xp, xpp)       
#elif defined(CHA)
      ! for channel
      call xDistribution_CHA(gridStretchX, len_x, xsize(1), x, xp, xpp)
#endif   
    call yDistribution(ny_global, y, dy)
    call zDistribution(zmesh_type, xstart(3), xsize(3) ,nz_global, len_z &
                       ,z_1, z_2, bumpz1, bumpz2, zplus_min, zplus_max &
                       ,zpert_1, zpert_2, bumpzpert1, bumpzpert2, zpluspert_min &
                       ,z, zp, zpp, z_global, zp_global, zpp_global)

    call mpi_allreduce(z_global, tmp, nz_global,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    z_global = tmp
    call mpi_allreduce(zp_global, tmp, nz_global,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    zp_global = tmp
    call mpi_allreduce(zpp_global, tmp1, nz_global,real_type,MPI_MIN,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(zpp_global, tmp2, nz_global,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    zpp_global = 2*(tmp1+tmp2)
    ! calculate inverse to avoid divisions later
    xp = 1.0_mytype/xp    
    zp = 1.0_mytype/zp 
    ! write grid for visualization
    inquire (iolength=lenr) x(1)
    ! in planes folder
    if (nrank == 0) then
      open(11,file='output/planes/x.bin', status='REPLACE', access='direct', recl=nx_global*lenr)
      write(11,rec=1) x(1:nx_global)
      close(11)
      open(11,file='output/planes/y.bin', status='REPLACE', access='direct', recl=ny_global*lenr)
      write(11,rec=1) y(1:ny_global)
      close(11)
      open(11,file='output/planes/z.bin', status='REPLACE', access='direct', recl=nz_global*lenr)
      write(11,rec=1) z_global(1:nz_global) 
      close(11)
    endif
    ! in restart folder
    if ((nrank == 0) .AND. (ny_global > 1)) then
      open(11,file='output/restart/x.bin', status='REPLACE', access='direct', recl=nx_global*lenr)
      write(11,rec=1) x(1:nx_global)
      close(11)
      open(11,file='output/restart/y.bin', status='REPLACE', access='direct', recl=ny_global*lenr)
      write(11,rec=1) y(1:ny_global)
      close(11)
      open(11,file='output/restart/z.bin', status='REPLACE', access='direct', recl=nz_global*lenr)
      write(11,rec=1) z_global(1:nz_global) 
      close(11)
    endif
    ! as .txt file
    if (nrank.eq.0) then
      open(11,file = 'postproc/xgrid.txt')
        write(11,'(F12.6)') ReTau
        write(11,'(F12.6)') gridStretchX
        do i=1,xsize(1)
          write(11,'(i5,4F12.6)') i,x(i),xp(i),xpp(i)
        enddo
      close(11)
      open(11,file = 'postproc/ygrid.txt')
        do j=1,ny_global
          write(11,'(i5,1F12.6)') j,y(j)
        enddo
      close(11)
      open(11,file = 'postproc/zgrid.txt')
      write(11,'(F12.6)') z_1
      write(11,'(F12.6)') z_2
      write(11,'(F12.6)') bumpz1
      write(11,'(F12.6)') bumpz2
      write(11,'(F12.6)') zplus_min
      write(11,'(F12.6)') zplus_max
        do k=1,nz_global
          write(11,'(i5,4F15.8)') k,z_global(k),zp_global(k),zpp_global(k)
        enddo
      close(11)
    endif
    !$acc enter data copyin(x,xp,xpp,y,z,zp,zpp)
    !$acc enter data copyin(dx,dy,dz)
    !$acc enter data copyin(xstart,ystart,zstart,xend,yend,zend)
  end subroutine init_grid


! calculation of the x-stencil (wall-normal) for flat-plate
  subroutine xDistribution(xmesh_type, ReTau, stretchx, len_x, npts, x, xp, xpp)
    implicit none
    integer :: npts, i, ierr
    real(mytype) :: stretchx, len_x, yplFactx, ReTau, factx
    real(mytype), dimension(:) :: x, xp, xpp
    character(len=30) :: xmesh_type
    ! equidistant
    if ( xmesh_type == "equid" ) then
      do i=1,npts
        if (perBC(1) .eqv. .true.)  factx  = (i-0.5_mytype)/(npts)
        if (perBC(1) .eqv. .false.) factx  = (1.0_mytype*i-1.0_mytype)/(npts-1.0_mytype)
        x(i)   = factx*len_x
        xp(i)  = 1.0_mytype
        xpp(i) = 0.0_mytype
      enddo
    ! stretched towards wall
    else if ( xmesh_type == "non_equid" ) then 
      yplFactx = 0.6_mytype/ReTau*(npts-1)/(len_x)
      if (yplFactx .gt. 2.0_mytype) then
        if  (nrank .eq. 0) then
          write(*,*) "yplFactx:", yplFactx
          write(*,*) "ERROR! Wall-normal mesh distribution too fine for this ReTau!"
        endif
        call decomp_2d_finalize
        call mpi_finalize(ierr) 
        stop
      endif
      do i=1,npts
        factx  =  (1.0_mytype*i-1.0_mytype)/(npts-1.0_mytype)
        x(i)   = (factx*yplFactx + (1.0_mytype + tanh(stretchx*(factx-1)/2.0_mytype) / & 
                      (tanh(stretchx*0.5_mytype)))*(1.0_mytype-yplFactx))*len_x
        xp(i)  = (yplFactx + (1-yplFactx)*0.5_mytype*stretchx/tanh(0.5_mytype*stretchx)/cosh(stretchx*(factx-1)/2)**2)*len_x
        xpp(i) = (-0.5_mytype*(1-yplFactx)*(stretchx**2)*tanh(stretchx*(factx-1)/2) &
                 /tanh(0.5_mytype*stretchx)/(cosh(stretchx*(factx-1)/2))**2)*len_x
      enddo
    endif
  end subroutine


! calculation of the x-stencil (wall-normal) for channel
  subroutine xDistribution_CHA(stretchx, len_x, npts, x, xp, xpp)
    implicit none
    integer :: npts, i, ierr
    real(mytype) :: stretchx, len_x, factx
    real(mytype), dimension(:) :: x, xp, xpp
    ! stretched towards both walls
    do i=1,npts
      factx    =  (i-0.5_mytype)/(npts) - 0.5_mytype
      x(i)    =  0.5_mytype*(1.0+tanh(stretchx*factx)/tanh(stretchx*0.5_mytype))*len_x
      xp(i)   =  0.5_mytype*stretchx/tanh(stretchx*0.5_mytype)/cosh(stretchx*factx)**2*len_x
      xpp(i)  =  -stretchx**2*tanh(stretchx*factx)/tanh(stretchx*0.5_mytype)/cosh(stretchx*factx)**2*len_x
    enddo
    end subroutine


! calculation of the y-stencil (spanwise)
  subroutine yDistribution(npts, y, dy)
    implicit none
    integer :: npts, j
    real(mytype) :: dy, facty
    real(mytype), dimension(:) :: y
    ! Equidistant
    do j=1,npts
      facty= 1.0_mytype*j
      y(j)= facty*dy
    enddo
  end subroutine


! calculation of the z-stencil (streamwise)
  subroutine zDistribution(zmesh_type, xst, npts, kmax, len_z, z_1, z_2, bumpz1, bumpz2, zplus_min, zplus_max &
                          ,z_pert1, z_pert2, bumpzpert1, bumpzpert2, zpluspert_min &
                          ,z, zp, zpp, z_global, zp_global, zpp_global)
    implicit none
    integer :: npts, k, kk, kmax, xst
    real(mytype) :: len_z, z_1, z_2, bumpz1, bumpz2, zplus_min, zplus_max
    real(mytype) :: z_pert1, z_pert2, bumpzpert1, bumpzpert2, zpluspert_min
    real(mytype) :: factz, delta1, delta2, deltapert1, deltapert2, z_int1, z_max
    real(mytype), dimension(:) :: z, zp, zpp, z_global, zp_global, zpp_global
    character(len=30) :: zmesh_type
    z_global   = 0.0_mytype
    zp_global  = 0.0_mytype
    zpp_global = 0.0_mytype
    ! equidistant
    if ( zmesh_type == "equid" ) then
      do k=1,npts
        kk=(1.0_mytype*k+xst-1.0_mytype)
        if (perBC(3) .eqv. .true.)    factz = (kk-0.5_mytype)/(kmax)
        if (perBC(3) .eqv. .false.)   factz = (kk-1.0_mytype)/(kmax-1.0_mytype)
        z(k) = factz*len_z
        zp(k) = 1.0_mytype
        zpp(k) = 0.0_mytype
        z_global(kk)=z(k)
        zp_global(kk)=zp(k)
        zpp_global(kk)=zpp(k)
      enddo
    ! stretched
    else if ( zmesh_type == "non_equid" ) then
      delta1=z_1*bumpz1
      delta2=z_2*bumpz2
      deltapert1=z_pert1*bumpzpert1
      deltapert2=z_pert2*bumpzpert2
      z_int1=deltapert1*0.5_mytype*(zpluspert_min-zplus_max)*log(cosh((z_pert1)/deltapert1)) &
            +deltapert2*0.5_mytype*(zplus_max-zpluspert_min)*log(cosh((z_pert2)/deltapert2)) &
            +delta1*0.5_mytype*(zplus_min-zplus_max)*log(cosh((z_1)/delta1)) &
            +delta2*0.5_mytype*(zplus_max-zplus_min)*log(cosh((z_2)/delta2))
      do k=1,npts
        kk = k + xst - 1
        factz = (kk - 1.0_mytype)/(kmax-1.0_mytype)
        z(k)  = ( delta1*0.5_mytype*(zplus_min-zplus_max)*log(cosh((z_1-factz)/delta1)) & 
                + delta2*0.5_mytype*(zplus_max-zplus_min)*log(cosh((z_2-factz)/delta2)) & 
                + deltapert1*0.5_mytype*(zpluspert_min-zplus_max)*log(cosh((z_pert1-factz)/deltapert1)) & 
                + deltapert2*0.5_mytype*(zplus_max-zpluspert_min)*log(cosh((z_pert2-factz)/deltapert2)) & 
                + factz*zplus_max-z_int1 )
        zp(k) = ( zplus_min+0.5_mytype*(zplus_max-zplus_min)*( 2.0_mytype-tanh((factz-z_1)/(delta1))+tanh((factz-z_2)/(delta2)) ) &
                + 0.5_mytype*(zplus_max-zpluspert_min)*( -tanh((factz-z_pert1)/(deltapert1))+tanh((factz-z_pert2)/(deltapert2)) ) )
        zpp(k)= 0.5_mytype*(zplus_max-zplus_min) &
                 *( 1.0_mytype/delta2/(cosh((z_2-factz)/delta2))**2 - 1.0_mytype/delta1/(cosh((z_1-factz)/delta1))**2 ) &
                + 0.5_mytype*(zplus_max-zpluspert_min) &
                 *( 1.0_mytype/deltapert2/(cosh((z_pert2-factz)/deltapert2))**2 &
                    - 1.0_mytype/deltapert1/(cosh((z_pert1-factz)/deltapert1))**2 )
        z_global(kk)   = z(k)
        zp_global(kk)  = zp(k)
        zpp_global(kk) = zpp(k)
      enddo 
      z_max   =  ( delta1*0.5_mytype*(zplus_min-zplus_max)*log(cosh((z_1-1.0_mytype)/delta1)) & 
                 + delta2*0.5_mytype*(zplus_max-zplus_min)*log(cosh((z_2-1.0_mytype)/delta2)) & 
                 + deltapert1*0.5_mytype*(zpluspert_min-zplus_max)*log(cosh((z_pert1-1.0_mytype)/deltapert1)) & 
                 + deltapert2*0.5_mytype*(zplus_max-zpluspert_min)*log(cosh((z_pert2-1.0_mytype)/deltapert2)) & 
                 + 1.0_mytype*zplus_max-z_int1 )
      z   = z/z_max*len_z
      zp  = zp/z_max*len_z
      zpp = zpp/z_max*len_z
      z_global   = z_global/z_max*len_z
      zp_global  = zp_global/z_max*len_z
      zpp_global = zpp_global/z_max*len_z
    endif
  end subroutine


! display the grid settings and for pert_calc==1, calculate the fixed timestep
  subroutine print_init_grid()
    use decomp_2d
    use mod_param
    use mod_math
    implicit none
    integer :: k, ierr
    integer :: p_total, ngrid_total, index_dz_min
    real(mytype) :: cphase, alpha_wl, wl, npoints_x_wl1, npoints_x_wl2, dz1, dz2, dz_min, yplus, freq1
    real(mytype), allocatable, dimension(:) :: dz_real
    allocate(dz_real(nz_global))
    p_total=p_row*p_col
    ngrid_total=nx_global*ny_global*nz_global
    ! calculation of streamwise grid spacing
    dz1=z_global(2)  
    dz_real=1.0e30_mytype
    do k=2,nz_global
      dz_real(k)=z_global(k)-z_global(k-1)
    enddo
    index_dz_min=minloc(dz_real,1)
    dz2=dz_real(index_dz_min)
    ! calculation of a new time step when perturbation is active
    if ((pert_calc==1) .OR. (BC_inl == "inlet_lst")) then
      freq1=pert_F(1)
      omega1=Re*freq1
      ! multiple of the amount of samples during each forcing period
      Kmult = fft_samples*fft_step
      dt_FFT = 2.0_mytype*pi_const/(omega1*Kmult)
      ! assumption: classical incompressible phase-speed 
      cphase=0.3_mytype
      alpha_wl=omega1/cphase
      wl=2*pi_const/alpha_wl
      npoints_x_wl1=wl/dz1
      npoints_x_wl2=wl/dz2
    endif
    ! calculation of yplus at the domain inlet
    yplus = x(2)*ReTau
    if (xmesh_type == "non_equid") then
      if (yplus .gt. 2.0_mytype) then
        if  (nrank .eq. 0) then
          write(stdout,*) "yplus:", yplus
          write(stdout,*) "ERROR! yplus too large! Reduce ReTau"
        endif
        call decomp_2d_finalize
        call mpi_finalize(ierr) 
        stop
      endif
    endif

    ! display settings
    if (nrank == 0) then
      write(stdout,* ) 'Mesh'
      write(stdout,'(A)') 'o--------------------------------------------------o'
      if (nrank == 0) then
        write(stdout,'(A, F10.4)') 'Grid initialisation                        done!'
      endif
      if (ny_global == 1) then
        write(stdout,'(A)') 'Simulation grid dimensions                    2D'
      else
        write(stdout,'(A)') 'Simulation grid dimensions                    3D'
      end if
      write(stdout,'(A, I10)') 'Total number of domains:              ',p_total
      write(stdout,'(A, I10)') 'Domains in z-direction (x-y-plane):   ',p_col
      write(stdout,'(A, I10)') 'Domains in y-direction (x-z-plane):   ',p_row
      write(stdout,* ) 
      write(stdout,'(A, F10.4)') 'Length domain in z-direction:         ',len_z
      write(stdout,'(A, F10.4)') 'Length domain in y-direction:         ',len_y
      write(stdout,'(A, F10.4)') 'Length domain in x-direction:         ',len_x
      write(stdout,* ) 
      write(stdout,'(A, A10)') 'Wall-normal mesh:                          ',xmesh_type
      write(stdout,'(A, I10)') 'Grid points nz (streamwise):          ',nz_global
      write(stdout,'(A, A10)') 'Streamwise mesh:                           ',zmesh_type
      if (zmesh_type == "equid") then
        write(stdout,'(A, F10.4)') 'Grid size:                            ',dz
#if defined(BL)
        write(stdout,'(A, F10.4)') 'z_plus in laminar region:             ',z(2)*ReTau
#endif
        if (pert_calc==1) then
          write(stdout,'(A, F10.4)') 'Characteristic wavelength:            ',wl
          write(stdout,'(A, F10.4)') 'Number of points per wavelength:      ',npoints_x_wl1
        endif
      else if (zmesh_type == "non_equid") then
        write(stdout,'(A, F10.4)') 'Laminar region length:                ',z_1*len_z
        write(stdout,'(A, F10.4)') 'Grid size dz1 in laminar region:      ',dz1
        write(stdout,'(A, F10.4)') 'z_plus in laminar region:             ',z(2)*ReTau
        write(stdout,'(A, F10.4)') 'Turbulent region length:              ',len_z-z_1*len_z-(1-z_2)*len_z
        write(stdout,'(A, F10.4)') 'Grid size dz2 in turbulent region:    ',dz2
        if (pert_calc==1) then
          write(stdout,'(A, F10.4)') 'Characteristic wavelength:            ',wl
          write(stdout,'(A, F10.4)') 'Points per wavelength (laminar):      ',npoints_x_wl1
          write(stdout,'(A, F10.4)') 'Points per wavelength (turbulent):    ',npoints_x_wl2
        endif
      endif
      write(stdout,'(A, I10)') 'Grid points ny (spanwise):            ',ny_global  
      write(stdout,'(A, F10.4)') 'Grid size dy                          ',dy    
      write(stdout,'(A, I10)') 'Grid points nx (wall-normal):         ',nx_global
      if (xmesh_type == "non_equid") then
        write(stdout,'(A, F10.4)') 'Wall distance to first point:         ',x(2)
        write(stdout,'(A, F10.4)') 'y_plus at domain inlet                ',yplus
        write(stdout,'(A, F10.4)') 'Stretching factor                     ',gridStretchX
      else
        write(stdout,'(A, F10.4)') 'Grid size dx                          ',dx
      endif
      write(stdout,'(A, I10)') 'Total number of grid points           ',ngrid_total
      write(stdout,* ) 
      write(stdout,* ) 'Grid transformations:                             '
      if (xmesh_type == "equid") then
        write(stdout,* ) 'no transformation in x-direction (equidistant) '
      else if (xmesh_type == "non_equid") then
        write(stdout,* ) 'x-direction (0,len_x) -> eta-direction (0,1) '
      endif 
      if (zmesh_type == "equid") then
        write(stdout,* ) 'no transformation in z-direction (equidistant) '  
      else if (zmesh_type == "non_equid") then
        write(stdout,* ) 'z-direction (0,len_z) -> xi-direction (0,1)  '
      endif 
      write(stdout,'(A)') 'o--------------------------------------------------o'
      write(stdout,* )
    endif
    deallocate(dz_real)
  end subroutine


end module mod_grid
