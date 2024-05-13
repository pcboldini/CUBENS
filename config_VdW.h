
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
    INIT  = "initField_CHIT" ! initBoundaryLayer, initField_CHIT

    USE_EOS   = "PR"
    USE_VISC  = "Chung"   ! only PowerLaw

    Re   = 1600.0
    Pra  = 0.75
    Ma   = 0.1
    Ec   = 0.004
    
! ------------ BOUNDARY CONDITIONS  ------------------------------
!
!   wall BCs:
!   adiab_std,  isoth_std  : standard implementation of wall BC
!   adiab_nrbc, isoth_nrbc : non-reflecting implementation of wall BC

    perBC  = (/.true.,.true.,.true./)

! ------------ TIME & OUTPUT FILING  ------------------
!  
    CFL              = 0.5
    dtMax            = 3.385e-4
    nsteps           = 3
    readRestartFile  = -10000
    intvSaveRestart  = 10000
    intvSavePlanes   = 100
    saveRestartAfter = 100
    savePlanesAfter  = 100
    intvCalcCFL      = 100
    intvReadParam    = 100
    intvPrint        = 1

    ! ------------ set reference free-stream values for computation ------------------
    Tref =0.9000000000
    Pref =1.1000000000
    Rhoref =1.8047114354
    Cpref =8.0239456825
    SOSref =2.7658411983
    Rref =2.666667
! ------------ set viscosity and conductivity ------------------
    Muref =1.6400035224e-03
    Kref =3.5314649926e-02

! ----------- MPI ------------------
!
    p_row = 1 ! y-direction (spanwise)
    p_col = 1 ! z-direction (streamwise)


! ------------ GRID  ------------------
!
    imax = 256
    jmax = 256
    kmax = 256
    
    len_x = 2.0*pi  ! wall-normal
    len_y = 2.0*pi  ! spanwise
    len_z = 2.0*pi  ! spanwise

    xmesh_type = "equid" ! equid, non_equid  
    zmesh_type = "equid" ! equid, non_equid  

    gridStretchX = 5.0  ! clustering in x-direction


! ------------ OUTPUT CUT PLANES  ------------------
!    
    xi_plane = pi
    yi_plane = pi
    zi_plane = pi



