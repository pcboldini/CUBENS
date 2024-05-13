
!----------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                                  !
!      config.h                                                                                                                    !
!                                                                                                                                  !
!      DESCRIPTION:                                                                                                                !
!      ------------                                                                                                                !
!      This header contains general DNS variables.                                                                                 !
!                                                                                                                                  !
!                                                                                                                                  !
!      CHANGELOG:                                                                                                                  !
!      ----------                                                                                                                  !
!         xx.xx.2021: header created (Rene)                                                                                        !
!         31.12.2021: header reorganised (Pietro)                                                                                  !
!         09.02.2022: multiple EOS added (Pietro)                                                                                  ! 
!         06.04.2022: flag USE_EOS removed, EOS selection in initBL (Pietro)                                                       ! 
!                                                                                                                                  !    
!----------------------------------------------------------------------------------------------------------------------------------!

! ------------ FLOW INPUT  ------------------------------
!
    CASE  = "BoundaryLayer" ! initBoundaryLayer, initField_CHIT
!
! ------------ BOUNDARY CONDITIONS  ------------------------------
!
!   wall BCs:
!   adiab_std,  isoth_std  : standard implementation of wall BC
!   adiab_nrbc, isoth_nrbc : non-reflecting implementation of wall BC
!   temp_nrbc, temp_std : wall heating implementation of wall BC

    BC_bot = "adiab_nrbc" ! 
 
!   freestream BCs:
    BC_top = "free_nrbc" ! free_nrbc

!   inlet BCs:
    BC_inl = "inlet_nrbc" ! inlet_nrbc, inlet_std, inlet_lst
    alphaLST = (0.061251_mytype,0.013715_mytype) ! only with inlet_lst
    epsilonLST = 1.0e-2_mytype ! only with inlet_lst

!   outlet BCs:
    BC_out = "outlet_nrbc" ! outlet_nrbc 

    perBC  = (/.false.,.true.,.false./)
!
! ------------ TIME & OUTPUT FILING  ------------------
!  
    CFL              = 0.8
    nsteps           = 11
    readRestartFile  = -11
    intvSaveRestart  = 2
    intvSavePlanes   = 11
    savePlanesAfter  = 6
    saveRestartAfter = 6
    intvCalcCFL      = 3
    intvReadParam    = -2
    intvPrint        = 2
!
! ------------ FFT ------------------         
!
    fft_samples = 30 ! per period (usually 2 periods enough)
    fft_step    = 1000 ! per sample
!
! ----------- MPI ------------------
!
    p_row = 1 ! y-direction (spanwise)
    p_col = 1 ! z-direction (streamwise)
!
! ---------- NUMERICS --------------
! 
    nStencilConv = 3
    nStencilVisc = 2    
!
! ------------ GRID  ------------------
!
    imax = 300
    jmax = 1
    kmax = 100
    
    len_x = 20.0 ! wall-normal
    len_y = 9.63  ! spanwise

    Redelta_start = 250.0 ! start Re_delta
    Redelta_end = 500   ! end Re_delta

    xmesh_type = "non_equid" ! equid, non_equid
    gridStretchX = 5.0  ! clustering in x-direction
    ReTau = 39.295139312744141 ! Re_tau at the inlet

    zmesh_type = "equid" ! equid, non_equid

    zplus_max = 20.0 ! zplus in laminar region

    zplus_min = 10.0 ! zplus in turbulent region 
    z_1 = 0.42 ! location of first bump in z-direction
    z_2 = 0.92 ! location of second bump in z-direction
    bumpz1 = 0.1 ! width of first bump in z-direction
    bumpz2 = 0.04 ! width of second bump in z-direction
    
    zpluspert_min = 10.0 ! zplus in strip region 
    zpert_1 = 0.12 ! location of first strip bump   
    zpert_2 = 0.18 ! location of second strip bump  
    bumpzpert1 = 0.04 ! width of first strip bump 
    bumpzpert2 = 0.04 ! width of second strip bump 
     
!
! ------------ OUTPUT CUT PLANES  ------------------
!    
    yi_plane = 1
    xi_plane = 25
    zi_plane = -1
!
! ------------ SPONGE & FILTER & MASS DIFFUSION  ------------------
!
!   Inlet: 
    spInlLen = 20.0 ! length in z
    spInlStr = 0.5 ! strength
    spInlExp = 2.0 ! exponent=2

!   Outlet: 
    spOutLen = 20.0 ! length in z
    spOutStr = 0.5 ! strength
    spOutExp = 2.0 ! exponent=2

!   Free stream: 
    spTopLen = 1.0 ! length in x
    spTopStr = 0.5 ! strength
    spTopExp = 2.0 ! exponent=2

!
! ------------ BLOWING/SUCTION  ------------------------
!
    pert_calc=0 ! 1 on, 0 off, 0 for inlet_lst

    pert_zLen = 10
    pert_ReMid = 320 ! as Re_delta

    beta0 = 2.0*pi_const/len_y

    pert_ampl = (/5.0e-7,5.0e-7/)
    pert_F    = (/100.000e-6,62.000e-6/)

!Sayadi

    pert_beta = (/beta0/)

!Franko&Lele
    
!    pert_beta = (/0.0_mytype,beta0/)
    pert_zSig  = 0.7
    pert_ySig  = 0.7   
!
! ------------ WALL HEATING  ------------------------
!
    Twall_new = 1.30     ! with respect to Tinf
    ReTem1 = 550
    ReTem2 = 700
    deltatem1 = 0.02
    deltatem2 = 0.02 
!
! ----------- POSTPRO ------------------
!
    p_row_pp  = 1 ! fixed at 1 so far
    p_col_pp  = 1
    istart_pp = 0
    iend_pp   = 10
    istep_pp  = 1

    fft_flag = 1 ! 1 on, 0 off
    index_fft_span = (/0,1,2/)

    rms_flag = 0 ! 1 on, 0 off
!
! ------------ INTERPOLATION ------------------------
!
    p_row_intp = 1 ! fixed at 1 so far
    p_col_intp = 1 

    timeStepRead = 10 
    timeStepSave = 11

    yi_plane_new = 1
    xi_plane_new = 25
    zi_plane_new = -1

    inew = imax
    jnew = 20
    knew = kmax

    len_x_new = len_x
    len_y_new = len_y
    Redelta_end_new = Redelta_end

    ReTau_new = ReTau

    xmesh_type_new = "non_equid" ! equid, non_equid
    gridStretchX_new = gridStretchX  ! clustering in x-direction

    zmesh_type_new = "equid" ! equid, non_equid

    z_1_new = 0.25 ! location of first bump in z-direction
    z_2_new = 0.90 ! location of second bump in z-direction 
    bumpz1_new = 0.2 ! width of first bump in z-direction
    bumpz2_new = 0.04 ! width of second bump in z-direction
    zplus_min_new = 10.0 ! zplus in turbulent region 
    zplus_max_new = 20.0 ! zplus in laminar region                

