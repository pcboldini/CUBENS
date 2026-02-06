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
    Re   = 1532.13                            ! Bulk Reynolds number
    Pra  = 0.71                                     ! Prandtl number 
    Ma   = 0.2                                      ! Mach number 
    ig_gam  = 1.4                                   ! ratio of specific heat
    eos_Rgas = 1.0/(ig_gam*Ma**2.0)                 ! specific gas constant
    Ec   = (ig_gam-1.0)*Ma**2                       ! Eckert number 
    Pref = 1/ig_gam/Ma**2                           ! Pressure
    Tinf = 273                                      ! Reference temperature for Sutherland
    Ri_unit = 0.5                                   ! Richardson number
!   
! --------------------------- BOUNDARY CONDITIONS  ------------------------
!
!   bottom BCs:
!   adiab_std,  isoth_std  : standard implementation of wall BC
    
    BC_bot = "adiab_std"                            ! see above

!   top BCs:
    BC_top = "adiab_std"                            ! see above

!   inlet BCs:
    BC_inl = "isoth_std"                            ! options: isoth_std
    Twall_inl   = 1.6                               ! Temperature inlet

!   outlet BCs:
    BC_out = "isoth_std"                            ! options: isoth_std
    Twall_out   = 0.4                               ! Temperature outlet

    perBC  = (/.false.,.true.,.false./)             ! if .true., periodic BC
!
! --------------------------------- TIME ----------------------------------
!  
    CFL              = 0.8                          ! CFL number, keep < 1
    dtMax            = -3.385e-4                    ! maximum timestep, if <0 ignored
    nsteps           = 2000000                      ! simulations steps
    intvCalcCFL      = 1                            ! interval CFL calculation
!
! ---------------------------------- I/O ---------------------------------- 
!     
    intvPrint        = 100                          ! on the screen
    readRestartFile  = -1                           ! if <0, not read
    intvSaveRestart  = 1000000                      ! interval save restart files
    saveRestartAfter = 1                            ! save restart files after which timestep
    intvSavePlanes   = 1000000                      ! interval save planes
    savePlanesAfter  = 1                            ! save planes after which timestep
    intvSaveStats    = -100                         ! interval save statistics
    saveStatsAfter   = -1                           ! save statistics after which timestep
    intvReadParam    = -2                           ! read variation file

    xi_plane = (/-1/)                               ! y-index cut z-x plane 
    yi_plane = (/1/)                                ! x-index cut z-y plane 
    zi_plane = (/-1/)                               ! z-index cut y-x plane 
!
! -------------------------------- NUMERICS -------------------------------
! 
    nStencilConv = 3                                ! order/2 convective fluxes
    nStencilVisc = 2                                ! order/2 diffusive fluxes 
    keep_flag = "classic"                           ! KEEP-scheme: "classic" or "pep" (pressure equilibrium)
!
! --------------------------------- GRID ----------------------------------
!
    imax = 301                                      ! number of points in x
    jmax = 1                                        ! number of points in y
    kmax = 301                                      ! number of points in z
    
    len_x = 1.0                                     ! wall-normal length
    len_y = 0.1                                     ! spanwise length
    len_z = 1.0                                     ! streamwise length

    xmesh_type = "equid"                            ! option: "equid"
    zmesh_type = "equid"                            ! option: "equid"
