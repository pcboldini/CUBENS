! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -                                                                                                             
! This header contains input DNS variables. Flow parameters are to be set in /preproc/initCHA/  
!
! ------------------------------ FLOW INPUT  ------------------------------
!
    Re   = 1600                                     ! Bulk Reynolds number
    dpdz = 2.0                                      ! Initial pressure gradient (adjusted during calculation to ensure u_bulk=1)
!
! --------------------------- BOUNDARY CONDITIONS  ------------------------
!
!   wall BCs:
!   adiab_std,  isoth_std  : standard implementation of wall BC
!   adiab_nrbc, isoth_nrbc : non-reflecting implementation of wall BC

    BC_bot = "isoth_nrbc"                           ! see above

!   freestream BCs:
    BC_top = "isoth_nrbc"                           ! see above

    perBC  = (/.false.,.true.,.true./)              ! if .true., periodic BC
!
! --------------------------------- TIME ----------------------------------
!  
    CFL              = 0.8                          ! CFL number, keep < 1
    dtMax            = -3.385e-4                    ! maximum timestep, if <0 ignored
    nsteps           = 2                            ! simulations steps
    intvCalcCFL      = 1                            ! interval CFL calculation
!
! ---------------------------------- I/O ---------------------------------- 
!     
    intvPrint        = 1                            ! on the screen
    readRestartFile  = -11                          ! if <0, not read
    intvSaveRestart  = 10                           ! interval save restart files
    saveRestartAfter = 10                           ! save restart files after which timestep
    intvSavePlanes   = 1                            ! interval save planes
    savePlanesAfter  = 0                            ! save planes after which timestep
    intvSaveStats    = 1                            ! interval save statistics
    saveStatsAfter   = 1                            ! save statistics after which timestep
    intvReadParam    = -2                           ! read variation file

    yi_plane = (/1/)                                ! y-index cut z-x plane          
    xi_plane = (/1/)                                ! x-index cut z-y plane        
    zi_plane = (/1/)                                ! z-index cut y-x plane
!
! -------------------------------- NUMERICS -------------------------------
! 
    nStencilConv = 3                                ! order/2 convective fluxes
    nStencilVisc = 2                                ! order/2 diffusive fluxes  
!
! --------------------------------- GRID ----------------------------------
!
    imax = 300                                      ! number of points in x
    jmax = 10                                       ! number of points in y
    kmax = 1000                                     ! number of points in z
    
    len_x = 2.0                                     ! wall-normal length
    len_y = 3.0                                     ! spanwise length
    len_z = 10.0                                    ! streamwise length

    xmesh_type = "non_equid"                        ! options: "equid", "non_equid"
    gridStretchX = 5.0                              ! wall clustering in x-direction
    ReTau = 85.0                                    ! Re_tau at the domain inlet

    zmesh_type = "equid"                            ! option: "equid"    
!
! ----------------------------- ONLY FOR POST -----------------------------
!
    p_row_pp  = 1                                   ! keep at 1
    p_col_pp  = 1                                   ! number of streamwise partitions (needs to match the launch command)

    istart_pp = 0                                   ! "restart": start timestep for averaging
    iend_pp   = 10                                  ! "restart": end timestep for averaging
    istep_pp  = 1                                   ! "restart": averaging timestep

    rms_flag = 0                                    ! options: 1 on, 0 off

    fft_flag = 1                                    ! options: 1 on, 0 off 
    index_fft_span = (/0,1,2/)                      ! spanwise modes
!
! --------------------------- ONLY FOR INTERPOL ---------------------------
!
    p_row_intp = 1                                  ! keep at 1 
    p_col_intp = 1                                  ! number of streamwise partitions (need to match the launch command)

    timeStepRead = 10                               ! timestep of reading restart file
    timeStepSave = 11                               ! timestep of writing restart file (can be = timeStepRead but then overwritten)

    yi_plane_new = (/1/)                            ! y-index cut z-x plane
    xi_plane_new = (/25/)                           ! x-index cut z-y plane 
    zi_plane_new = (/-1/)                           ! z-index cut y-x plane

    inew = imax                                     ! new number of points in x
    jnew = 20                                       ! new number of points in y
    knew = kmax                                     ! new number of points in z

    len_x_new = len_x                               ! new wall-normal length
    len_y_new = len_y                               ! new spanwise length

    ReTau_new = ReTau                               ! new Re_tau at the domain inlet

    xmesh_type_new = "non_equid"                    ! options: "equid", "non_equid"
    gridStretchX_new = gridStretchX                 ! new wall clustering in x-direction

    zmesh_type_new = "equid"                        ! option: "equid"           

