! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
! auxiliary subroutines for parameters calculation and I/O
module mod_auxl
  use decomp_2d
  use mod_param
  use mod_grid
  use mod_finitediff
  implicit none
  contains
  ! calculation of the vorticity
  subroutine calcVort(vortx,vorty,vortz,u,v,w) 
    use decomp_2d
    implicit none
    integer i,j,k,c
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: vortx,vorty,vortz
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: u,v,w
    real(mytype) :: dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz, xa, za
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dwz,dwy,dwx,dvz,dvy,dvx,duz,duy,dux)
    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          dux   = 0.0_mytype
          duy   = 0.0_mytype
          duz   = 0.0_mytype
          dvx   = 0.0_mytype
          dvy   = 0.0_mytype
          dvz   = 0.0_mytype
          dwx   = 0.0_mytype
          dwy   = 0.0_mytype
          dwz   = 0.0_mytype
          xa = xp(i)
          za = zp(k)
          !$acc loop seq
          do c = 1, nStencilVisc
            dux = dux + visc_ddx(c)*(u(i+c,j,k) - u(i-c,j,k))*xa
            duy = duy + visc_ddy(c)*(u(i,j+c,k) - u(i,j-c,k))
            duz = duz + visc_ddz(c)*(u(i,j,k+c) - u(i,j,k-c))*za
            dvx = dvx + visc_ddx(c)*(v(i+c,j,k) - v(i-c,j,k))*xa
            dvy = dvy + visc_ddy(c)*(v(i,j+c,k) - v(i,j-c,k))
            dvz = dvz + visc_ddz(c)*(v(i,j,k+c) - v(i,j,k-c))*za
            dwx = dwx + visc_ddx(c)*(w(i+c,j,k) - w(i-c,j,k))*xa
            dwy = dwy + visc_ddy(c)*(w(i,j+c,k) - w(i,j-c,k))
            dwz = dwz + visc_ddz(c)*(w(i,j,k+c) - w(i,j,k-c))*za
          enddo
          ! x: wall-normal, y: spanwise, z: streamwise
          vortx(i,j,k) = dwy - dvz ! du/dz - dw/dx -> vorticity wall-normal
          vorty(i,j,k) = duz - dwx ! dv/dx - du/dy -> vorticity spanwise
          vortz(i,j,k) = dvx - duy ! dw/dy - dv/dz -> vorticity streamwise 
        enddo
      enddo
    enddo
  end subroutine
! ????????
  subroutine cmpbulkvel(r,w,wb)
    use decomp_2d
    use mod_param
    use mpi
    use mod_grid
    implicit none 
    integer ierr,i,j,k, nxnynz
    real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: w,r
    real(mytype) wb, da
    wb  = 0.0_mytype
    !$acc parallel loop collapse(3) default(present) private(da) reduction(+:wb)
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=2,xsize(1)
             da = x(i)-x(i-1)
             wb = wb  + 0.5_mytype*(r(i,j,k)*w(i,j,k)+r(i-1,j,k)*w(i-1,j,k))*da
          enddo
       enddo
    enddo
    nxnynz = jmax*kmax*len_x
    call mpi_allreduce(MPI_IN_PLACE, wb, 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr);   wb = wb/nxnynz
  end subroutine cmpbulkvel
! calculation of bulk properties for quick monitoring of the simulation
  subroutine cmpbulk(istep,wt1,time,dt,CFL_new,rho,u,v,w,ien,pre,tem,mu,ka,vortx,vorty,vortz,wb)
    use decomp_2d
    use mod_param
    use mpi
    use mod_grid
    use mod_eos
    implicit none 
    integer ierr,i,j,k,istep,ioutput,nxnynz
    real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka,vortx,vorty,vortz
    real(mytype), dimension(:,:,:), allocatable :: Mach
    real(mytype) ddx,ddy,ddz,time,dt,isNan,isNanGlobal,CFL_new
    real(mytype) rb,wb,preb,kib,eb,str,tmp,enst,sos,max_Mach
    real(8) :: wt1
    allocate(Mach(xsize(1), xsize(2), xsize(3)))
    ! initialization
    isNan = 0.0_mytype
    isNanGlobal = 0.0_mytype
    rb   = 0.0_mytype
    wb   = 0.0_mytype
    preb = 0.0_mytype
    kib  = 0.0_mytype
    eb   = 0.0_mytype
    str  = 0.0_mytype
    enst = 0.0_mytype
    ! domain volume
    nxnynz = len_x*len_z*len_y
    ! grid spacing in y- and z-direction
    ddy=y(2)-y(1)
    ddz=z(2)-z(1)
    ! calculation of bulk properties
    !$acc parallel loop collapse(3) default(present) private(ddx,ddy,ddz) reduction(+:mub,rb,wb,preb,rwb,kib,eb,enst) async(1)
    do k=1,xsize(3)
      do j=1,xsize(2)
          do i=2,xsize(1)
            ddx  = x(i)-x(i-1)
            ! density
            rb  = rb  + 0.5*(rho(i,j,k)+rho(i-1,j,k))*ddx
            ! streamwise velocity
            wb  = wb  + 0.5*(w(i,j,k)+w(i-1,j,k))*ddx
            ! pressure
            preb = preb + 0.5*(pre(i,j,k)+pre(i-1,j,k))*ddx
            ! kinetic energy
            kib = kib + 0.5*rho(i,j,k)*(w(i,j,k)**2+v(i,j,k)**2+u(i,j,k)**2)*ddx
            ! internal energy
            eb  = eb  + 0.5*(ien(i,j,k)+ien(i-1,j,k))**2*ddx
            ! enstrophy
            enst = enst + mu(i,j,k)*(vortx(i,j,k)**2 + vorty(i,j,k)**2 + vortz(i,j,k)**2)*ddx
            call calcSOS(rho(i,j,k),ien(i,j,k),sos)
            ! Mach number
            Mach(i,j,k)=sqrt(w(i,j,k)**2+v(i,j,k)**2+u(i,j,k)**2)/sos
          enddo
       enddo
    enddo
    ! Wall shear stress
    !$acc parallel loop collapse(2) default(present) reduction(+:str) async(1)
    do k=1,xsize(3)
       do j=1,xsize(2)
          ! approximation
          str = str + (mu(2,j,k)*w(2,j,k)/x(2)) 
       enddo
    enddo
    kib=kib*ddy*ddz
    enst=enst*ddy*ddz
    !$acc wait
    ! combine all MPI processes
    call mpi_allreduce(MPI_IN_PLACE, rb  , 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr);   rb   = rb/nxnynz
    call mpi_allreduce(MPI_IN_PLACE, wb  , 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr);   wb   = wb/nxnynz
    call mpi_allreduce(MPI_IN_PLACE, preb  , 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr); preb = preb/nxnynz
    call mpi_allreduce(MPI_IN_PLACE, kib , 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr);   kib  = kib/nxnynz
    call mpi_allreduce(MPI_IN_PLACE, eb  , 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr);   eb   = eb/nxnynz
    call mpi_allreduce(MPI_IN_PLACE, str , 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr);   str  = str/nxnynz
    call mpi_allreduce(MPI_IN_PLACE, enst, 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr);   enst = enst/nxnynz
    call mpi_allreduce(MPI_IN_PLACE, Mach, 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr);   
    ! calculate maximal Mach number
    max_Mach=maxval(Mach)
    ! print the step, timestep, CFL, and the previously calculated bulk quantities
    if  (nrank .eq. 0) then
       if ((mod(istep,20*intvPrint).eq.0)) then 
          write(*,'(A7,13A16)') 'step', 'wall-time', 'dt',  'CFL', 'time', 'wall-stress', &
                                'bulk w-vel', 'bulk rho', 'bulk P', 'ke-bulk', 'eint-bulk' , 'enstrophy', 'max_Mach'
       endif
       if (mod(istep,intvPrint).eq.0) write(*,'(I7,13E16.7)') &
                     istep, MPI_WTIME()-wt1, dt, CFL_new, time, str, wb, rb, preb, kib, eb, enst, max_Mach
    endif
    ! check if there is any NaN, in case stop the simulation
    if (wb .ne. wb) isNan = 1
    call mpi_allreduce(isNan, isNanGlobal,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    if (isNanGlobal>0) then
      call decomp_2d_finalize
      call mpi_finalize(ierr)
      stop
    endif
    ! deallocation
    deallocate(Mach)
  end subroutine
! I/O planes
  subroutine output2dPlane(part,nHaloIn,istep,dir,loc,tmp01,name011,name012,tmp02,name021,name022,tmp03,name031,name032,tmp04,   &
                                                      name041,name042,tmp05,name051,name052,tmp06,name061,name062,tmp07,name071, &
                                                      name072,tmp08,name081,name082,tmp09,name091,name092,tmp10,name101,name102, &
                                                      tmp11,name111,name112,tmp12,name121,name122)
    ! -- dir specifies the direction of the plane:
    ! - 1 --> x-normal plane
    ! - 2 --> y-normal plane
    ! - 3 --> z-normal plane
    ! loc is the location of the plane in global coordinates
    use decomp_2d
    use decomp_2d_io
    use mod_param
    implicit none
    integer :: istep,dir,loc,nHaloIn
    character*7 cha
    character(len=1024) :: cha2
    character(len=*), optional :: name011,name021,name031,name041,name051,name061,name071,name081,name091,name101,name111,name121
    character(len=*), optional :: name012,name022,name032,name042,name052,name062,name072,name082,name092,name102,name112,name122
    real(mytype), dimension(1-nHaloIn:,1-nHaloIn:,1-nHaloIn:), optional :: &
                        tmp01,tmp02,tmp03,tmp04,tmp05,tmp06,tmp07,tmp08,tmp09,tmp10,tmp11,tmp12
    TYPE (DECOMP_INFO) :: part                    
    real(mytype), allocatable, dimension(:,:,:) :: tmp
    allocate(tmp(part%xsz(1),part%xsz(2),part%xsz(3)))
    write(cha,'(I0.7)') istep
    write(cha2,'(I0)') loc
    ! all the planes are written in ./postproc/planes/
    if(present(tmp01)) then 
      tmp = tmp01(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','postproc/planes/'//trim(name011)//trim(cha2)//'.'//trim(name012)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp02)) then 
      tmp = tmp02(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))  
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','postproc/planes/'//trim(name021)//trim(cha2)//'.'//trim(name022)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp03)) then 
      tmp = tmp03(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','postproc/planes/'//trim(name031)//trim(cha2)//'.'//trim(name032)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp04)) then 
      tmp = tmp04(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3)) 
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','postproc/planes/'//trim(name041)//trim(cha2)//'.'//trim(name042)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp05)) then 
      tmp = tmp05(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','postproc/planes/'//trim(name051)//trim(cha2)//'.'//trim(name052)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp06)) then 
      tmp =  tmp06(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3)) 
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','postproc/planes/'//trim(name061)//trim(cha2)//'.'//trim(name062)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp07)) then 
      tmp =  tmp07(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3)) 
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','postproc/planes/'//trim(name071)//trim(cha2)//'.'//trim(name072)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp08)) then 
      tmp =  tmp08(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','postproc/planes/'//trim(name081)//trim(cha2)//'.'//trim(name082)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp09)) then 
      tmp =  tmp09(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','postproc/planes/'//trim(name091)//trim(cha2)//'.'//trim(name092)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp10)) then 
      tmp =  tmp10(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','postproc/planes/'//trim(name101)//trim(cha2)//'.'//trim(name102)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp11)) then 
      tmp =  tmp11(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','postproc/planes/'//trim(name111)//trim(cha2)//'.'//trim(name112)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp12)) then 
      tmp =  tmp12(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','postproc/planes/'//trim(name121)//trim(cha2)//'.'//trim(name122)//cha//'.bin', &
                                 'dummy',part)
    endif
    if (nrank == 0) then
      if (istep == 0) then
        write(stdout,* ) 'Initial planes'
        write(stdout,* ) '-------------------------'
        write(stdout,* ) 'Calculation and writing                              completed'
        write(stdout,* ) 
      else if ( (istep > 0) .AND. (dir==1)) then 
        write(stdout,* ) 'Writing x-plane at step                             ', istep
      else if ( (istep > 0) .AND. (dir==2)) then 
        write(stdout,* ) 'Writing y-plane at step                             ', istep
      else if ( (istep > 0) .AND. (dir==3)) then
        write(stdout,* ) 'Writing z-plane at step                             ', istep
      endif
    endif
    deallocate(tmp)
  end subroutine
! I/O load restart
  subroutine loadRestart(istep,time,rho,u,v,w,ien,nHaloIn,part)
    use mpi
    use decomp_2d
    use decomp_2d_io
    use mod_param
    implicit none
    integer :: istep,ierr,nHaloIn
    character*7 cha
    real(mytype), dimension(1-nHaloIn:,1-nHaloIn:,1-nHaloIn:) :: rho,u,v,w,ien
    real(mytype), allocatable, dimension(:,:,:) :: tmp
    real(mytype), dimension(5) :: infoRestart
    real(mytype) :: time
    LOGICAL :: file_exists
    integer :: fh
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    TYPE(DECOMP_INFO), intent(IN) :: part
    file_exists = .FALSE.
    allocate(tmp(part%xsz(1),part%xsz(2),part%xsz(3)))
    write(cha,'(I0.7)') istep
    ! all the restart are located in ./restart
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'restart/ruvwe.'//cha//'.bin', &
                       MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
    call MPI_FILE_GET_SIZE(fh,filesize,ierr)
    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_read_scalar(fh,disp,5,infoRestart)
    INQUIRE (FILE='restart/ruvwe.'//cha//'.bin', EXIST=file_exists)
    ! check if restart file exists
    if (file_exists) then
      if(nrank==0) write(stdout,*) 'restart FILE exists!'
    else
      if(nrank==0) write(stdout,*) '!restart FILE does not exist!'
      call MPI_FILE_CLOSE(fh,ierr)
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      stop
    end if
    ! check if the restart file precision is correct
    if (infoRestart(1).ne.1.0*mytype) then 
      if (nrank.eq.0) write(stdout,*) 'Simulation stopped!'
      if (nrank.eq.0) write(stdout,*) 'Restart file has different precision.'
      call MPI_FILE_CLOSE(fh,ierr)
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      stop
    endif
    ! check if the restart file has the correct size
    if ((infoRestart(2).ne.1.0*nx_global) .or. &
        (infoRestart(3).ne.1.0*ny_global) .or. &
        (infoRestart(4).ne.1.0*nz_global)) then 
      if (nrank.eq.0) write(stdout,*) 'Simulation stopped!'
      if (nrank.eq.0) write(stdout,*) 'Restart has imax/jmax/kmax : ', nint(infoRestart(2)), &
                                                                       nint(infoRestart(3)), &
                                                                       nint(infoRestart(4))
      if (nrank.eq.0) write(stdout,*) 'This sim has imax/jmax/kmax : ', nx_global,ny_global,nz_global
      call MPI_FILE_CLOSE(fh,ierr)
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      stop
    endif 
    ! open the five primary non-conservative variables: density, velocities, and internal energy
    call decomp_2d_read_var(fh,disp,1,tmp,part); rho(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3)) = tmp
    call decomp_2d_read_var(fh,disp,1,tmp,part);   u(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3)) = tmp
    call decomp_2d_read_var(fh,disp,1,tmp,part);   v(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3)) = tmp
    call decomp_2d_read_var(fh,disp,1,tmp,part);   w(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3)) = tmp
    call decomp_2d_read_var(fh,disp,1,tmp,part); ien(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3)) = tmp
    call MPI_FILE_CLOSE(fh,ierr)
    deallocate(tmp)
    time = infoRestart(5)
    if(nrank==0) write(stdout,'(A31,E12.5,A12,I10)') 'The simulation restarts at t =', time, 'and at step', istep
    if(nrank==0) write(stdout,*) ''
  end subroutine
! I/O save restart
  subroutine saveRestart(istep,time,rho,u,v,w,ien,nHaloIn,part,interpol)
    use mpi
    use decomp_2d
    use decomp_2d_io
    use mod_param
    implicit none
    integer :: istep,ierr,nHaloIn
    character*7 cha
    real(mytype), dimension(1-nHaloIn:,1-nHaloIn:,1-nHaloIn:) :: rho,u,v,w,ien
    real(mytype), allocatable, dimension(:,:,:) :: tmp
    real(mytype), dimension(5) :: infoRestart
    real(mytype) :: time
    integer :: fh
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    TYPE(DECOMP_INFO), intent(IN) :: part
    character(len=*), optional :: interpol
    ! read the five infoRestart arguments
    infoRestart(1) = 1.0*mytype
    ! distinguish interpolated file
    if(interpol == 'aI') then 
      infoRestart(2) = 1.0*inew
      infoRestart(3) = 1.0*jnew
      infoRestart(4) = 1.0*knew
    else
      infoRestart(2) = 1.0*imax
      infoRestart(3) = 1.0*jmax
      infoRestart(4) = 1.0*kmax
    endif
    infoRestart(5) = time 
    allocate(tmp(part%xsz(1),part%xsz(2),part%xsz(3)))
    write(cha,'(I0.7)') istep
    if (nrank == 0) then
      if (istep == 0) then
        write(stdout,* ) 
        write(stdout,* ) 'Initial restart condition'
        write(stdout,* ) '-------------------------'
        write(stdout,* ) 'Calculation and writing                              completed'
        write(stdout,* ) 
      else if ( istep > 0 ) then 
        write(stdout,* ) 'Writing restart at step                             ', istep    
      endif
    endif
    ! all the restart files are written in ./restart
    if(interpol == 'bI') then 
      if (nrank == 0) write(stdout,* ) 'before interpolation'
      call MPI_FILE_OPEN(MPI_COMM_WORLD, 'restart/ruvwe.'//trim(interpol)//'.'//cha//'.bin', &
                        MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
    else
      call MPI_FILE_OPEN(MPI_COMM_WORLD, 'restart/ruvwe.'//cha//'.bin', &
                        MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
    endif
    if (nrank == 0) then
      write(stdout,* ) '-------------------------'
      write(stdout,* ) 
    endif
    ! writing of the restart file
    filesize = 0_MPI_OFFSET_KIND
    ! guarantee overwriting
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)  
    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_write_scalar(fh,disp,5,infoRestart)
    ! write the five primary non-conservative variables: density, velocities, and internal energy
    tmp = rho(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3));  call decomp_2d_write_var(fh,disp,1,tmp,part)
    tmp =   u(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3));  call decomp_2d_write_var(fh,disp,1,tmp,part)
    tmp =   v(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3));  call decomp_2d_write_var(fh,disp,1,tmp,part)
    tmp =   w(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3));  call decomp_2d_write_var(fh,disp,1,tmp,part)
    tmp = ien(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3));  call decomp_2d_write_var(fh,disp,1,tmp,part)
    call MPI_FILE_CLOSE(fh,ierr)
    deallocate(tmp)
  end subroutine
end module
