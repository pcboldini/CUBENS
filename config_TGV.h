! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -                                                                                                             
! This header contains input DNS variables. 
!
! ------------------------------ FLOW INPUT  ------------------------------
!
    Re   = 1600.0                                   ! Bulk Reynolds number
    Pra  = 0.75                                     ! Prandtl number 
    Ma   = 0.1                                      ! Mach number 
    ig_gam  = 1.4                                   ! ratio of specific heat
    eos_Rgas = 1.0/(ig_gam*Ma**2.0)                 ! specific gas constant
    Ec   = (ig_gam-1.0)*Ma**2                       ! Eckert number 
    Pref = 1/ig_gam/Ma**2                           ! Pressure
    Tinf = 300                                      ! Reference temperature for Sutherland
!   
! --------------------------- BOUNDARY CONDITIONS  ------------------------
!
!   wall BCs:
!   adiab_std,  isoth_std  : standard implementation of wall BC
!   adiab_nrbc, isoth_nrbc : non-reflecting implementation of wall BC
    perBC  = (/.true.,.true.,.true./)               ! if .true., periodic BC
!
! --------------------------------- TIME ----------------------------------
!  
    CFL              = 0.5                          ! CFL number, keep < 1
    dtMax            = -3.385e-4                    ! maximum timestep, if <0 ignored
    nsteps           = 10                            ! simulations steps
    intvCalcCFL      = 1                            ! interval CFL calculation
!
! ---------------------------------- I/O ---------------------------------- 
!     
    intvPrint        = 1                            ! on the screen
    readRestartFile  = -11                          ! if <0, not read
    intvSaveRestart  = 100                           ! interval save restart files
    saveRestartAfter = 100                           ! save restart files after which timestep
    intvSavePlanes   = 1                            ! interval save planes
    savePlanesAfter  = 1                            ! save planes after which timestep
    intvSaveStats    = 100                            ! interval save statistics
    saveStatsAfter   = -1                           ! save statistics after which timestep
    intvReadParam    = -2                           ! read variation file

    xi_plane = (/1/)                         ! y-index cut z-x plane 
    yi_plane = (/1/)                         ! x-index cut z-y plane 
    zi_plane = (/1/)                         ! z-index cut y-x plane 
!
! -------------------------------- NUMERICS -------------------------------
! 
    nStencilConv = 3                                ! order/2 convective fluxes
    nStencilVisc = 2                                ! order/2 diffusive fluxes 
    keep_flag = "pep"                               ! KEEP-scheme: "classic" or "pep" (pressure equilibrium)
!
! --------------------------------- GRID ----------------------------------
!
    imax = 64                                       ! number of points in x
    jmax = 64                                       ! number of points in y
    kmax = 64                                       ! number of points in z
    
    len_x = 2.0*pi_const                            ! wall-normal length
    len_y = 2.0*pi_const                            ! spanwise length
    len_z = 2.0*pi_const                            ! streamwise length

    xmesh_type = "equid"                            ! option: "equid"
    zmesh_type = "equid"                            ! option: "equid"
