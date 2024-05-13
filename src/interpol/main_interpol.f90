
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

  integer :: ierr,i,j,k,lenr,h
  real(mytype) :: time
  real(8) :: wt_start
  character*2 :: cha
  CHARACTER(100) :: commLinePar1
  CHARACTER(100) :: commLinePar2
  real(mytype), allocatable, dimension(:) :: tmp, tmp_old

  character(len=30) :: date
  

  ! SCRINS Version number
   real(mytype), parameter                    :: version = 1.1                      ! SCRINS version

!===============================================================================================!
!
!      INITIALIZATION
!
!===============================================================================================!

  call fdate(date) 


! users can specify the number of cores (p_row and p_col) also from the command line
  if (COMMAND_ARGUMENT_COUNT().eq.2) then
    call GET_COMMAND_ARGUMENT(1, commLinePar1)
    call GET_COMMAND_ARGUMENT(2, commLinePar2)
    read(commLinePar1,*) p_row_intp
    read(commLinePar2,*) p_col_intp
  endif

  call mpi_init(ierr)

  ! Reading parameters
  call read_config()
  select case (CASE)
    case("BoundaryLayer")
      call read_initBL_params()
  end select

  ! ---- initialize the first partition
  call decomp_2d_init(imax,jmax,kmax,p_row_intp,p_col_intp)

  ! ---- copy the main partion 
  call get_decomp_info(part1)

  ! ---- get partition for the new grid
  call decomp_info_init(inew,jnew,knew,part2)

!-----------------------------!
! Print a welcome message     !
!-----------------------------!

if (nrank == 0) then

  write (stdout, *)
  write (stdout, *) "o----------------------------------------------------------------------------------o"
  write (stdout, *) "|           ________  ________  ________  ___  ________   ________                 |"
  write (stdout, *) "|          |\   ____\|\   ____\|\   __  \|\  \|\   ___  \|\   ____\                |"
  write (stdout, *) "|          \ \  \___|\ \  \___|\ \  \|\  \ \  \ \  \\ \  \ \  \___|                |"
  write (stdout, *) "|           \ \_____  \ \  \    \ \   _  _\ \  \ \  \\ \  \ \_____  \              |"
  write (stdout, *) "|            \|____|\  \ \  \____\ \  \\  \\ \  \ \  \\ \  \|____|\  \             |"
  write (stdout, *) "|              ____\_\  \ \_______\ \__\\ _\\ \__\ \__\\ \__\____\_\  \            |"
  write (stdout, *) "|             |\_________\|_______|\|__|\|__|\|__|\|__| \|__|\_________\           |"
  write (stdout, *) "|             \|_________|                                  \|_________|           |"
  write (stdout, *) "|                                                                                  |"
  write (stdout, *) "|                                                                                  |"
  write (stdout, *) "|                        SuperCRItical Navier-Stokes solver                        |"
  write (stdout, *) "|                                                                                  |"
  write (stdout, '(A,F4.2,A,A,A)') " | Version ",version,": ",date,"                                     |"                                                                         
  write (stdout, *) "|                                                                                  |" 
  write (stdout, *) "o----------------------------------------------------------------------------------o"
  write (stdout, *)
  write (stdout, '(A,A)') " compiled with ", compiler_version()
  write (stdout, *)
  write (stdout, '(A,I0,A,I0,A,I0,A)') " Using ",nproc, " MPI, ",p_row_intp," number of procs in z and ",p_col_intp," in y"
  write (stdout, *)
  write (stdout, *) "o----------------------------------------------------------------------------------o"
  write (stdout, *)
  write (stdout, *) "|                                     INTERPOL                                     |"
  write (stdout, *)
  write (stdout, *) "o----------------------------------------------------------------------------------o"
  write (stdout, *)

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(tmp_old(kmax))
allocate(tmp(knew))

wt_start = MPI_WTIME()

  ! init the domain of partition 1
  call initDomain(xmesh_type,part1,imax,jmax,kmax,len_x,len_y,len_z,  &
                  ReTau,gridStretchX, &
                  zmesh_type, z_1, z_2, bumpz1, bumpz2, zplus_min, zplus_max, &
                  zpert_1, zpert_2, bumpzpert1, bumpzpert2, zpluspert_min, &
                  xold,yold,zold,zold_global,rbo,ubo,vbo,wbo,ebo)

  call mpi_allreduce(zold_global, tmp_old, kmax,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
  zold_global = tmp_old

  if (nrank==0) then
    write(stdout,* ) 'OLD mesh:                                        '
    write(stdout,* ) 
    write(stdout,* ) 'Simulation domain in z-direction                 ',len_z
    write(stdout,* ) 'Simulation domain in y-direction                 ',len_y
    write(stdout,* ) 'Simulation domain in x-direction                 ',len_x
    write(stdout,* )
    write(stdout,* ) 'Grid points nx (wall-normal)                     ',imax 
    write(stdout,* ) 'Grid points ny (spanwise)                        ',jmax 
    write(stdout,* ) 'Grid points nx (streamwise)                      ',kmax 
    write(stdout,* )
    write(stdout,* ) 'Wall-normal mesh:                                ',xmesh_type
    write(stdout,* ) 'Stretching factor                                ',gridStretchX
    write(stdout,* ) 'Streamwise mesh:                                 ',zmesh_type
    write(stdout,* ) 
  endif

  ! init the domain of partition 2
  call initDomain(xmesh_type_new, part2, inew,jnew,knew,len_x_new, len_y_new, len_z_new, &
                  ReTau_new, gridStretchX_new, &
                  zmesh_type_new, z_1_new, z_2_new, bumpz1_new, bumpz2_new, zplus_min_new, zplus_max_new, &
                  zpert_1, zpert_2, bumpzpert1, bumpzpert2, zpluspert_min, &
                  xnew,ynew,znew,znew_global,rbn,ubn,vbn,wbn,ebn)

  call mpi_allreduce(znew_global, tmp, knew,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
  znew_global = tmp


  if (nrank==0) then
    write(stdout,* ) 'NEW mesh:                                        '
    write(stdout,* ) 
    write(stdout,* ) 'Simulation domain in z-direction                 ',len_z_new
    write(stdout,* ) 'Simulation domain in y-direction                 ',len_y_new
    write(stdout,* ) 'Simulation domain in x-direction                 ',len_x_new
    write(stdout,* )
    write(stdout,* ) 'Grid points nx (wall-normal)                     ',inew 
    write(stdout,* ) 'Grid points ny (spanwise)                        ',jnew
    write(stdout,* ) 'Grid points nx (streamwise)                      ',knew
    write(stdout,* )
    write(stdout,* ) 'Wall-normal mesh:                                ',xmesh_type_new
    write(stdout,* ) 'Stretching factor                                ',gridStretchX_new
    write(stdout,* ) 'Streamwise mesh:                                 ',zmesh_type_new
    write(stdout,* )
  endif

  ! if (nrank==0) then
  !   write(stdout,*) "part_old:", part1%xsz(1), part1%xsz(2), part1%xsz(3)
  !   write(stdout,*) "part_new:", part2%xsz(1), part2%xsz(2), part2%xsz(3)
  !   write(stdout,*)
  ! endif

  ! initialize EOS model
  if (nrank==0) then
    write(stdout,*) "Initializing EoS"
    write(stdout,*)
  endif
  call init_EOSModel()

  ! initialize visc model
  if (nrank==0) then
    write(stdout,*) "Initializing viscosity and conductivity"
    write(stdout,*)
  endif
  call init_VISCModel()

! initialize EOS and VISC
  if (nrank==0) then
    write(stdout,*) 'Initializing EOS and model for transport coefficients'
    write(stdout,*)
  endif
  call init_PARAM_EOS()



  call initSolution(part2,xnew,znew,rbn,ubn,vbn,wbn,ebn)

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

  if(nrank.eq.0) then
    write(stdout,'(A)',ADVANCE='NO') 'Reading old restart files:'
  endif

  call loadRestart(timeStepRead,time,rbo,ubo,vbo,wbo,ebo,nHaloInterpol,part1)
  call saveRestart(timeStepRead,time,rbo,ubo,vbo,wbo,ebo,nHaloInterpol,part1,'bI')


  if (yi_plane.gt.0)  then
    call output2dPlane(part1,nHaloInterpol,timeStepRead, 2, yi_plane, &
         rbo,'ypl_bI.','r.', ubo,'ypl_bI.','u.', vbo,'ypl_bI.','v.', wbo,'ypl_bI.','w.', ebo,'ypl_bI.','e.')
  endif

  if (xi_plane.gt.0)  then  
    call output2dPlane(part1,nHaloInterpol,timeStepRead, 1, xi_plane, &
         rbo,'xpl_bI.','r.', ubo,'xpl_bI.','u.', vbo,'xpl_bI.','v.', wbo,'xpl_bI.','w.', ebo,'xpl_bI.','e.')
  endif

  if (zi_plane.gt.0)  then  
    call output2dPlane(part1,nHaloInterpol,timeStepRead, 3, zi_plane, &
         rbo,'zpl_bI.','r.', ubo,'zpl_bI.','u.', vbo,'zpl_bI.','v.', wbo,'zpl_bI.','w.', ebo,'zpl_bI.','e.')
  endif

  inquire (iolength=lenr) xnew(1)
  if (nrank == 0) then
    open(11,file='postproc/planes/x_I.bin', status='REPLACE', access='direct', recl=inew*lenr)
    write(11,rec=1) xnew(1:inew)
    close(11)

    open(11,file='postproc/planes/y_I.bin', status='REPLACE', access='direct', recl=jnew*lenr)
    write(11,rec=1) ynew(1:jnew)
    close(11)

    open(11,file='postproc/planes/z_I.bin', status='REPLACE', access='direct', recl=knew*lenr)
    write(11,rec=1) znew_global(1:knew)
    close(11)
  endif

    if ((nrank == 0) .AND. (jnew > 1)) then
      open(11,file='restart/x_I.bin', status='REPLACE', access='direct', recl=inew*lenr)
      write(11,rec=1) xnew(1:inew)
      close(11)

      open(11,file='restart/y_I.bin', status='REPLACE', access='direct', recl=jnew*lenr)
      write(11,rec=1) ynew(1:jnew)
      close(11)

      open(11,file='restart/z_I.bin', status='REPLACE', access='direct', recl=knew*lenr)
      write(11,rec=1) znew_global(1:knew) 
      close(11)
    endif


! start the interpolation from source to destination
  write(stdout,* )

  if(nrank==0) write(stdout,*) "interpolating density"
  call interp3D(rbo,rbn, imax,jmax,kmax, inew,jnew,knew)

  if(nrank==0) write(stdout,*) "interpolating u-vel"
  call interp3D(ubo,ubn, imax,jmax,kmax, inew,jnew,knew)

  if(nrank==0) write(stdout,*) "interpolating v-vel"
  call interp3D(vbo,vbn, imax,jmax,kmax, inew,jnew,knew)

  if(nrank==0) write(stdout,*) "interpolating w-vel"
  call interp3D(wbo,wbn, imax,jmax,kmax, inew,jnew,knew)

  if(nrank==0) write(stdout,*) "interpolating total energy"
  call interp3D(ebo,ebn, imax,jmax,kmax, inew,jnew,knew)

  write(stdout,* ) 

  ! write the new interpolated restart file and planes

  if(nrank.eq.0) write(stdout,*) "Writing new restart files" 
  call saveRestart(timeStepSave,time,rbn,ubn,vbn,wbn,ebn,nHaloInterpol,part2,'aI')

  if (yi_plane_new.gt.0)  then
    call output2dPlane(part2,nHaloInterpol,timeStepSave, 2, yi_plane_new, &
         rbn,'ypl_aI.','r.', ubn,'ypl_aI.','u.', vbn,'ypl_aI.','v.', wbn,'ypl_aI.','w.', ebn,'ypl_aI.','e.')
  endif

  if (xi_plane_new.gt.0)  then
    call output2dPlane(part2,nHaloInterpol,timeStepSave, 1, xi_plane_new, &
         rbn,'xpl_aI.','r.', ubn,'xpl_aI.','u.', vbn,'xpl_aI.','v.', wbn,'xpl_aI.','w.', ebn,'xpl_aI.','e.')
  endif

  if (zi_plane_new.gt.0)  then
    call output2dPlane(part2,nHaloInterpol,timeStepSave, 3, zi_plane_new, &
         rbn,'zpl_aI.','r.', ubn,'zpl_aI.','u.', vbn,'zpl_aI.','v.', wbn,'zpl_aI.','w.', ebn,'zpl_aI.','e.')
  endif

  if (nrank == 0) write(*,*) "interpolation DONE!"
  if (nrank == 0) print '("Total time = ",f10.3," minutes.")', (MPI_WTIME() - wt_start)/60.0

  call decomp_2d_finalize()
  call mpi_finalize(ierr)

end program





