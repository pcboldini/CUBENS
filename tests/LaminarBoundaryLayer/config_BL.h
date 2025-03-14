! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -                                                                                                             
! This header contains input DNS variables. Flow parameters are to be set in /preproc/initBL/                                                                                                                                                                                                            
!
! --------------------------- BOUNDARY CONDITIONS  ------------------------
!
!   wall BCs:
!   adiab_std,  isoth_std  : standard implementation of wall BC
!   adiab_nrbc, isoth_nrbc : non-reflecting implementation of wall BC

    BC_bot = "isoth_nrbc"                           ! see above
 
!   freestream BCs:
    BC_top = "free_nrbc"                            ! options: free_nrbc

!   inlet BCs:
    BC_inl = "inlet_nrbc"                           ! options: inlet_nrbc, inlet_std

!   recycling-rescaling:
    BC_inl_rescale = .false.                         ! options: .true. or .false.
    z_recycle = 100.0_mytype                        ! recycle position
    delta_inl = 2.0_mytype                          ! inlet boundary thickness

!   outlet BCs:
    BC_out = "outlet_nrbc"                          ! options: outlet_nrbc (subsonic)

    perBC  = (/.false.,.true.,.false./)             ! if .true., periodic BC
!
! --------------------------------- TIME ----------------------------------
!  
    CFL              = 0.8                          ! CFL number, keep < 1
    dtMax            = -0.1                         ! maximum timestep, if <0 ignored (only applies to pert_bc=0)
    nsteps           = 2360000                       ! simulations steps
    intvCalcCFL      = 100                         ! interval CFL calculation
!
! ---------------------------------- I/O ----------------------------------
!     
    intvPrint        = 1000                         ! on the screen
    readRestartFile  = -11                          ! if <0, not read
    intvSaveRestart  = 2300000                        ! interval save restart files
    saveRestartAfter = 1                            ! save restart files after which timestep
    intvSavePlanes   = 1000                        ! interval save planes
    savePlanesAfter  = 2300000                            ! save planes after which timestep
    intvSaveStats    = -1                            ! interval save statistics
    saveStatsAfter   = -1                            ! save statistics after which timestep
    intvReadParam    = -2                           ! read variation file

    yi_plane = (/1/)                                ! y-index cut z-x plane  
    xi_plane = (/-1/)                             ! x-index cut z-y plane
    zi_plane = (/-1/)                                ! z-index cut y-x plane
!
! ------------------------- FFT (if pert_calc=1) --------------------------         
!
    fft_samples = 15                                ! per period (usually 2 periods enough)
    fft_step    = 1000                              ! per sample
!
! -------------------------------- NUMERICS -------------------------------
! 
    nStencilConv = 3                                ! order/2 convective fluxes
    nStencilVisc = 2                                ! order/2 diffusive fluxes
    keep_flag = "pep"
!
! --------------------------------- GRID ----------------------------------
!
    imax = 500                                     ! number of points in x
    jmax = 1                                       ! number of points in y
    kmax = 3000                                    ! number of points in z
    
    len_x = 40.0                                    ! wall-normal length
    len_y = 9.63                                    ! spanwise length

    Redelta_start = 700                          ! start Re_delta
    Redelta_end = 1500                            ! end Re_delta

    xmesh_type = "non_equid"                        ! options: "equid", "non_equid"
    gridStretchX = 8.0                              ! wall clustering in x-direction
    ReTau = 54                                    ! Re_tau at the domain inlet

    zmesh_type = "equid"                            ! options: "equid", "non_equid"

    zplus_max = 20.0                                ! non_equid: zplus in laminar region

    zplus_min = 10.0                                ! non_equid: zplus in turbulent region 
    z_1 = 0.42                                      ! non_equid: location of first bump in z-direction
    z_2 = 0.92                                      ! non_equid: location of second bump in z-direction
    bumpz1 = 0.1                                    !non_equid: width of first bump in z-direction
    bumpz2 = 0.04                                   ! non_equid: width of second bump in z-direction
    
    zpluspert_min = 10.0                            ! non_equid: zplus in strip region 
    zpert_1 = 0.12                                  ! non_equid: location of first strip bump   
    zpert_2 = 0.18                                  ! non_equid: location of second strip bump  
    bumpzpert1 = 0.04                               ! non_equid: width of first strip bump 
    bumpzpert2 = 0.04                               ! non_equid: width of second strip bump           
!
! -------------------------------- SPONGE ---------------------------------
!
!   Inlet: 
    spInlLen = 20.0                                  ! length in z
    spInlStr = 0.5                                  ! strength
    spInlExp = 2.0                                  ! exponent

!   Outlet: 
    spOutLen = 50.0                                 ! length in z
    spOutStr = 0.5                                  ! strength
    spOutExp = 2.0                                  ! exponent

!   Free stream: 
    spTopLen = 2.0                                  ! length in x
    spTopStr = 0.5                                  ! strength
    spTopExp = 2.0                                  ! exponent
!
! ---------------------------- BLOWING/SUCTION  ---------------------------
!
    pert_calc=1                                     ! options: 1 on, 0 off

    pert_zLen = 10                                  ! length perturbation
    pert_ReMid = 780                                ! Re_delta for mid point perturbation

    beta0 = 2.0*pi_const/len_y                      ! fundamental spanwise wavenumber

    pert_ampl = (/5.0e-7,0.0e-7/)                   ! amplitude 2-D and 3-D wave
    pert_F    = (/80.000e-6,1.000e-6/)            ! frequency 2-D and 3-D wave

    ! Sayadi et al., JFM 724, 2013

    pert_beta = (/beta0/)                           ! multiple spanwise modes

    ! Franko & Lele, JFM 730, 2013
    
!    pert_beta = (/0.0_mytype,beta0/)               ! multiple spanwise modes
    pert_zSig  = 0.7                                ! Gaussian standard deviation
    pert_ySig  = 0.7                                ! spanwise length
!
! ----------------------------- ONLY FOR POST -----------------------------
!
    p_row_pp  = 1                                   ! keep at 1
    p_col_pp  = 1                                   ! number of streamwise partitions (needs to match the launch command)
    avg_flag = "restart"                            ! options: "restart" or "stats"

    istart_pp = 5                                   ! "restart": start timestep for averaging
    iend_pp   = 10                                  ! "restart": end timestep for averaging
    istep_pp  = 1                                   ! "restart": averaging timestep

    stats_step      = (/10/)                        ! "stats": timestep at which averaging is obtained from ./simulate
    stats_time_rate = (/1/)                         ! "stats": time rate for different statistics timesteps

    rms_flag = 1                                    ! options: 1 on, 0 off
    istart_rms = 5                                  ! start RMS calculation       
    iend_rms   = 10                                 ! end RMS calculation         
    istep_rms  = 1                                  ! step RMS calculation         
    index_rms_xpl = (/1/)                           ! index x-plane for RMS calculation              
    index_rms_zpl = (/446/)                         ! index z-plane for RMS calculation  

    fft_flag = 1                                    ! options: 1 on, 0 off 
    index_fft_span = (/0,1,2/)                      ! spanwise modes

    ! in use when avg_flag = 'stats'; else istart_pp,iend_pp,istep_pp are used
    istart_fft = 5                                  ! "stats": start timestep FFT          
    iend_fft   = 10                                 ! "stats": end timestep FFT
    istep_fft  = 1                                  ! "stats": timestep for FFT                        
    index_fft_xpl = 45                              ! "stats": index x-plane for FFT, e.g. FFT at constant x             
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
    Redelta_end_new = Redelta_end                   ! new end local Reynolds number 

    ReTau_new = ReTau                               ! new Re_tau at the domain inlet

    xmesh_type_new = "non_equid"                    ! options: "equid", "non_equid"
    gridStretchX_new = gridStretchX                 ! new wall clustering in x-direction

    zmesh_type_new = "equid"                        ! options: "equid", "non_equid"

    z_1_new = 0.25                                  ! new location of first bump in z-direction
    z_2_new = 0.90                                  ! new location of second bump in z-direction 
    bumpz1_new = 0.2                                ! new width of first bump in z-direction
    bumpz2_new = 0.04                               ! new width of second bump in z-direction
    zplus_min_new = 10.0                            ! new zplus in turbulent region 
    zplus_max_new = 20.0                            ! new zplus in laminar region                

