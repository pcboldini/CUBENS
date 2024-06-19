! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini and the CUBENS contributors. All rights reserved.
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
    dtMax            = 3.385e-4                    ! maximum timestep, if <0 ignored
    nsteps           = 48000                            ! simulations steps
    intvCalcCFL      = 100                            ! interval CFL calculation
!
! ---------------------------------- I/O ---------------------------------- 
!     
    intvPrint        = 100                            ! on the screen
    readRestartFile  = -11                          ! if <0, not read
    intvSaveRestart  = -10000                           ! interval save restart files
    saveRestartAfter = -10000                           ! save restart files after which timestep
    intvSavePlanes   = -10000                            ! interval save planes
    savePlanesAfter  = -10000                            ! save planes after which timestep
    intvSaveStats    = -10000                            ! interval save statistics
    saveStatsAfter   = -1                           ! save statistics after which timestep
    intvReadParam    = -2                           ! read variation file

    xi_plane = (/pi_const/)                         ! y-index cut z-x plane 
    yi_plane = (/pi_const/)                         ! x-index cut z-y plane 
    zi_plane = (/pi_const/)                         ! z-index cut y-x plane 
!
! -------------------------------- NUMERICS -------------------------------
! 
    nStencilConv = 3                                ! order/2 convective fluxes
    nStencilVisc = 2                                ! order/2 diffusive fluxes 
!
! --------------------------------- GRID ----------------------------------
!
    imax = 128                                       ! number of points in x
    jmax = 128                                       ! number of points in y
    kmax = 128                                       ! number of points in z
    
    len_x = 2.0*pi_const                            ! wall-normal length
    len_y = 2.0*pi_const                            ! spanwise length
    len_z = 2.0*pi_const                            ! streamwise length

    xmesh_type = "equid"                            ! option: "equid"
    zmesh_type = "equid"                            ! option: "equid"
