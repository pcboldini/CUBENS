
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
!
! ------------ FLOW INPUT  ------------------------------
!
    CASE  = "TGV" ! BoundaryLayer, TGV, Channel

    USE_EOS   = "IG"
    USE_VISC  = "Sutherland"   ! only PowerLaw

    Re   = 1600.0
    Pra  = 0.75
    Ma   = 0.1
    ig_gam  = 1.4
    eos_Rgas = 1.0/(ig_gam*Ma**2.0)
    Ec   = (ig_gam-1.0)*Ma**2
    Pref = 1/ig_gam/Ma**2
    Tinf = 300
!   
! ------------ BOUNDARY CONDITIONS  ------------------------------
!
!   wall BCs:
!   adiab_std,  isoth_std  : standard implementation of wall BC
!   adiab_nrbc, isoth_nrbc : non-reflecting implementation of wall BC
    perBC  = (/.true.,.true.,.true./)
!
! ------------ TIME & OUTPUT FILING  ------------------
!  
    CFL              = 0.5
    dtMax            = 3.385e-4
    nsteps           = 10
    readRestartFile  = -1
    intvSaveRestart  = 10
    intvSavePlanes   = 10
    saveRestartAfter = 10
    savePlanesAfter  = 10
    intvCalcCFL      = 10
    intvReadParam    = 100
    intvPrint        = 1
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
    imax = 64
    jmax = 64
    kmax = 64
    
    len_x = 2.0*pi_const  ! wall-normal
    len_y = 2.0*pi_const  ! spanwise
    len_z = 2.0*pi_const  ! spanwise

    xmesh_type = "equid" ! equid, non_equid  
    zmesh_type = "equid" ! equid, non_equid  

    gridStretchX = 5.0  ! clustering in x-direction
!
! ------------ OUTPUT CUT PLANES  ------------------
!    
    xi_plane = pi_const
    yi_plane = pi_const
    zi_plane = pi_const



