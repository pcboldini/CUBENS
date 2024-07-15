! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
program interpolate
  use mod_interpolate
  use mod_eos
  use mod_eos_var
  use mod_eos_visc
  use mod_param
  use iso_fortran_env
  use decomp_2d
  use decomp_2d_io
  use decomp_2d_constants
  use factor
  implicit none


  !===============================================================================================!
  !                                       CUBENS interpol data
  !===============================================================================================!
  integer :: ierr,i,j,k,lenr,h
  real(mytype) :: time
  real(8) :: wt_start
  character*2 :: cha
  CHARACTER(100) :: commLinePar1
  CHARACTER(100) :: commLinePar2
  real(mytype), allocatable, dimension(:) :: tmp, tmp_old
  character(len=30) :: date
  ! CUBENS Version number
   real(mytype), parameter                    :: version = 1.0       


!===============================================================================================!
!                                       INITIALIZATION
!===============================================================================================!
! generate today's date
  call fdate(date) 
!---------------------------------------!
!  Read the config.h and init_BL files  !
!---------------------------------------!
! Reading parameters
  call read_config()
#if defined(BL) || defined(CHA)
  call read_init_params()
#endif
!-----------------------------!
!         MPI library         !
!-----------------------------!
! users can specify the number of cores (p_row and p_col) also from the command line
  if (COMMAND_ARGUMENT_COUNT().eq.2) then
    call GET_COMMAND_ARGUMENT(1, commLinePar1)
    call GET_COMMAND_ARGUMENT(2, commLinePar2)
    read(commLinePar1,*) p_row_intp
    read(commLinePar2,*) p_col_intp
  endif
  call mpi_init(ierr)
! initialize the first partition
  call decomp_2d_init(imax,jmax,kmax,p_row_intp,p_col_intp)
! copy the main partion 
  call get_decomp_info(part1)
! get partition for the new grid
  call decomp_info_init(inew,jnew,knew,part2)
! print a welcome message
  if (nrank == 0) then
    write (stdout, *)
    write (stdout, *) "o--------------------------------------------------o"
    write (stdout, *) "|         ________  ______  _______   ________     |"
    write (stdout, *) "|        / ____/ / / / __ )/ ____/ | /  / ___/     |"
    write (stdout, *) "|       / /   / / / / /_/ / __/ /  |/  //____      |"
    write (stdout, *) "|      / /___/ /_/ / /_/ / /___/ __   /___/ /      |"
    write (stdout, *) "|     /_____/_____/_____/_____/_/  |_//____/       |"
    write (stdout, *) "|                                                  |"
    write (stdout, *) "|                                                  |"
    write (stdout, *) "|       CUBic Equation of state Navier-Stokes      |"
    write (stdout, *) "|                                                  |"
    write (stdout, '(A,F4.2,A,A,A)') " |      Version ",version,": ",date,"|"                                                                         
    write (stdout, *) "|                                                  |" 
    write (stdout, *) "o--------------------------------------------------o"
    write (stdout, *)
    write (stdout, '(A,A)') " compiled with ", compiler_version()
    write (stdout, *)
    write (stdout, '(A,I0,A,I0,A,I0,A)') " Using ", nproc, " MPI processes, ", p_col_intp, &
                                         " number of procs in z and ", p_row_intp, " in y " 
    write (stdout, *)
    write (stdout, *) "o--------------------------------------------------o"
    write (stdout, *) "|                      INTERPOL                    |"
    write (stdout, *) "o--------------------------------------------------o"
    write (stdout, *)
  endif
! Print the initial simulation parameters
  call print_init_params()
! initialize EoS 
  if (nrank==0) then
    write(stdout,*) "Initializing EoS"
  endif
  call init_EOSModel()
! initialize transport properties
  if (nrank==0) then
    write(stdout,*) "Initializing viscosity and conductivity"
  endif
  call init_VISCModel()
! initialize EoS and TP global parameters
  if (nrank==0) then
    write(stdout,*) 'Initializing global variables for EOS and transport properties'
    write(stdout,*)
  endif
  call init_PARAM_EOS()


!===============================================================================================!
!                                       ALLOCATION
!===============================================================================================!
  allocate(tmp_old(kmax))
  allocate(tmp(knew))


!===============================================================================================!
!                                       CALCULATION
!===============================================================================================!
  wt_start = MPI_WTIME()  
! init the domain of the original partition 
  call initDomain(xmesh_type,part1,imax,jmax,kmax,len_x,len_y,len_z,  &
                  ReTau,gridStretchX, &
                  zmesh_type, z_1, z_2, bumpz1, bumpz2, zplus_min, zplus_max, &
                  zpert_1, zpert_2, bumpzpert1, bumpzpert2, zpluspert_min, &
                  xold,yold,zold,zold_global,rbo,ubo,vbo,wbo,ebo)
  call mpi_allreduce(zold_global, tmp_old, kmax,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
  zold_global = tmp_old
  if (nrank.eq.0) then
    write(stdout,'(A)') 'o--------------------------------------------------o'
    write(stdout,* ) 'OLD mesh:                                        '
    write(stdout,* ) 
    write(stdout,'(A, F10.4)') 'Length domain in z-direction:         ',len_z
    write(stdout,'(A, F10.4)') 'Length domain in y-direction:         ',len_y
    write(stdout,'(A, F10.4)') 'Length domain in x-direction:         ',len_x
    write(stdout,* )
    write(stdout,'(A, I10)') 'Grid points nx (wall-normal):         ',nx_global
    write(stdout,'(A, I10)') 'Grid points ny (spanwise):            ',ny_global     
    write(stdout,'(A, I10)') 'Grid points nz (streamwise):          ',nz_global
    write(stdout,* )
    write(stdout,'(A, A10)') 'Wall-normal mesh:                          ',xmesh_type
    if (xmesh_type == "non_equid") then
      write(stdout,'(A, F10.4)') 'Wall distance to first point:         ',xold(2)
      write(stdout,'(A, F10.4)') 'y_plus at domain inlet                ',xold(2)*ReTau
      write(stdout,'(A, F10.4)') 'Stretching factor                     ',gridStretchX
    else
      write(stdout,'(A, F10.4)') 'Grid size dx                          ',dx
    endif
    write(stdout,'(A, A10)') 'Streamwise mesh:                           ',zmesh_type
    write(stdout,* ) 
  endif
! init the domain of the new partition
  call initDomain(xmesh_type_new, part2, inew,jnew,knew,len_x_new, len_y_new, len_z_new, &
                  ReTau_new, gridStretchX_new, &
                  zmesh_type_new, z_1_new, z_2_new, bumpz1_new, bumpz2_new, zplus_min_new, zplus_max_new, &
                  zpert_1, zpert_2, bumpzpert1, bumpzpert2, zpluspert_min, &
                  xnew,ynew,znew,znew_global,rbn,ubn,vbn,wbn,ebn)
  call mpi_allreduce(znew_global, tmp, knew,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
  znew_global = tmp
  if (nrank.eq.0) then
    write(stdout,'(A)') 'o--------------------------------------------------o'
    write(stdout,* ) 'NEW mesh:                                        '
    write(stdout,* ) 
    write(stdout,'(A, F10.4)') 'Length domain in z-direction:         ',len_z_new
    write(stdout,'(A, F10.4)') 'Length domain in y-direction:         ',len_y_new
    write(stdout,'(A, F10.4)') 'Length domain in x-direction:         ',len_x_new
    write(stdout,* )
    write(stdout,'(A, I10)') 'Grid points nx (wall-normal):         ',inew
    write(stdout,'(A, I10)') 'Grid points ny (spanwise):            ',jnew     
    write(stdout,'(A, I10)') 'Grid points nz (streamwise):          ',knew
    write(stdout,* )
    write(stdout,'(A, A10)') 'Wall-normal mesh:                          ',xmesh_type_new
    if (xmesh_type == "non_equid") then
      write(stdout,'(A, F10.4)') 'Wall distance to first point:         ',xnew(2)
      write(stdout,'(A, F10.4)') 'y_plus at domain inlet                ',xnew(2)*ReTau_new
      write(stdout,'(A, F10.4)') 'Stretching factor                     ',gridStretchX_new
    else
      write(stdout,'(A, F10.4)') 'Grid size dx                          ',dx
    endif
    write(stdout,'(A, A10)') 'Streamwise mesh:                           ',zmesh_type_new
    write(stdout,* ) 
  endif
! initialize solution on the new partition (mesh)
  call initSolution(part2,xnew,znew,rbn,ubn,vbn,wbn,ebn)
  if (nrank.eq.0) then
    write(stdout,'(A)',ADVANCE='NO') 'Reading old restart files:'
    write(stdout,* ) 
  endif
! read and save the original solution, this is named with _bI
  call loadRestart(timeStepRead,time,rbo,ubo,vbo,wbo,ebo,nHaloInterpol,part1)
  call saveRestart(timeStepRead,time,rbo,ubo,vbo,wbo,ebo,nHaloInterpol,part1,'bI')
! write the old planes
  if (yi_plane(1).gt.0)  then
    do i=1,size(yi_plane)
      call output2dPlane(part1,nHaloInterpol,timeStepRead, 2, yi_plane(i), &
                         rbo,'ypl_bI.','r.', ubo,'ypl_bI.','u.', vbo,'ypl_bI.','v.', wbo,'ypl_bI.','w.', ebo,'ypl_bI.','e.')
      enddo
  endif
  if (xi_plane(1).gt.0)  then
    do i=1,size(xi_plane)  
      call output2dPlane(part1,nHaloInterpol,timeStepRead, 1, xi_plane(i), &
                         rbo,'xpl_bI.','r.', ubo,'xpl_bI.','u.', vbo,'xpl_bI.','v.', wbo,'xpl_bI.','w.', ebo,'xpl_bI.','e.')
    enddo
  endif
  if (zi_plane(1).gt.0)  then
    do i=1,size(zi_plane)  
      call output2dPlane(part1,nHaloInterpol,timeStepRead, 3, zi_plane(i), &
                         rbo,'zpl_bI.','r.', ubo,'zpl_bI.','u.', vbo,'zpl_bI.','v.', wbo,'zpl_bI.','w.', ebo,'zpl_bI.','e.')
    enddo
  endif
! write the new mesh
  inquire (iolength=lenr) xnew(1)
  if (nrank == 0) then
    open(11,file='output/planes/x_I.bin', status='REPLACE', access='direct', recl=inew*lenr)
    write(11,rec=1) xnew(1:inew)
    close(11)

    open(11,file='output/planes/y_I.bin', status='REPLACE', access='direct', recl=jnew*lenr)
    write(11,rec=1) ynew(1:jnew)
    close(11)

    open(11,file='output/planes/z_I.bin', status='REPLACE', access='direct', recl=knew*lenr)
    write(11,rec=1) znew_global(1:knew)
    close(11)
  endif
  if ((nrank == 0) .AND. (jnew > 1)) then
    open(11,file='output/restart/x_I.bin', status='REPLACE', access='direct', recl=inew*lenr)
    write(11,rec=1) xnew(1:inew)
    close(11)
    open(11,file='output/restart/y_I.bin', status='REPLACE', access='direct', recl=jnew*lenr)
    write(11,rec=1) ynew(1:jnew)
    close(11)
    open(11,file='output/restart/z_I.bin', status='REPLACE', access='direct', recl=knew*lenr)
    write(11,rec=1) znew_global(1:knew) 
    close(11)
  endif
! start the interpolation from source to destination
  if(nrank==0) write(stdout,*)
  if(nrank==0) write(stdout,*) "interpolating density"
  call interp3D(rbo,rbn, imax,jmax,kmax, inew,jnew,knew)
  if(nrank==0) write(stdout,*) "interpolating u-vel"
  call interp3D(ubo,ubn, imax,jmax,kmax, inew,jnew,knew)
  if(nrank==0) write(stdout,*) "interpolating v-vel"
  call interp3D(vbo,vbn, imax,jmax,kmax, inew,jnew,knew)
  if(nrank==0) write(stdout,*) "interpolating w-vel"
  call interp3D(wbo,wbn, imax,jmax,kmax, inew,jnew,knew)
  if(nrank==0) write(stdout,*) "interpolating internal energy"
  call interp3D(ebo,ebn, imax,jmax,kmax, inew,jnew,knew)
  if(nrank==0) write(stdout,*) 
! write the new interpolated restart file and planes
  if(nrank.eq.0) write(stdout,*) "writing the new restart files" 
  call saveRestart(timeStepSave,time,rbn,ubn,vbn,wbn,ebn,nHaloInterpol,part2,'aI')
! write the new interpolated planes
  if (yi_plane_new(1).gt.0)  then
    do i=1,size(yi_plane_new)
      call output2dPlane(part2,nHaloInterpol,timeStepSave, 2, yi_plane_new(i), &
                         rbn,'ypl_aI.','r.', ubn,'ypl_aI.','u.', vbn,'ypl_aI.','v.', wbn,'ypl_aI.','w.', ebn,'ypl_aI.','e.')
    enddo
  endif
  if (xi_plane_new(1).gt.0)  then
    do i=1,size(xi_plane_new)
      call output2dPlane(part2,nHaloInterpol,timeStepSave, 1, xi_plane_new(i), &
                         rbn,'xpl_aI.','r.', ubn,'xpl_aI.','u.', vbn,'xpl_aI.','v.', wbn,'xpl_aI.','w.', ebn,'xpl_aI.','e.')
    enddo
  endif
  if (zi_plane_new(1).gt.0)  then
    do i=1,size(zi_plane_new)
      call output2dPlane(part2,nHaloInterpol,timeStepSave, 3, zi_plane_new(i), &
                         rbn,'zpl_aI.','r.', ubn,'zpl_aI.','u.', vbn,'zpl_aI.','v.', wbn,'zpl_aI.','w.', ebn,'zpl_aI.','e.')
    enddo
  endif
! write text files with new mesh
  if (nrank.eq.0) then
    open(11,file = 'postproc/xgrid_I.txt')
      write(11,'(F12.6)') ReTau_new
      write(11,'(F12.6)') gridStretchX_new
      do i=1,inew
        write(11,'(i5,2F12.6)') i,xnew(i)
      enddo
    close(11)
  endif
  if (nrank.eq.0) then
    open(11,file = 'postproc/ygrid_I.txt')
      do j=1,jnew
        write(11,'(i5,1F12.6)') j,ynew(j)
      enddo
    close(11)
  endif
  if (nrank.eq.0) then
    open(11,file = 'postproc/zgrid_I.txt')
    write(11,'(F12.6)') z_1_new
    write(11,'(F12.6)') z_2_new
    write(11,'(F12.6)') bumpz1_new
    write(11,'(F12.6)') bumpz2_new
    write(11,'(F12.6)') zplus_min_new
    write(11,'(F12.6)') zplus_max_new
    do k=1,knew
      write(11,'(i5,2F15.8)') k,znew_global(k)
    enddo
    close(11)
  endif
  if (nrank==0) then
    write(*,*)
    write (stdout, *) "o--------------------------------------------------o"
    write(*,*) 'INTERPOLATION done!'
  endif
! print total time 
  if (nrank == 0) print '("Total time = ",f10.3," minutes.")', (MPI_WTIME() - wt_start)/60.0

  
  !===============================================================================================!
  !                                       DEALLOCATION
  !===============================================================================================!
  deallocate(tmp_old)
  deallocate(tmp)
  ! stop the mpi 
  call decomp_2d_finalize
  call mpi_finalize(ierr)
end program
