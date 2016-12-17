!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: MOD10A1_obsMod
! \label(MOD10A1_obsMod)
!
! !INTERFACE:
module MOD10A1_obsMod
! 
! !USES:   
  use ESMF
  use map_utils

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!
!   This subroutine provides the observation plugin for reading the 
!   MOD10A1 fractional snow cover product (from Terra). Note that 
!   this is a resampled data at 1km. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  7 Feb 2011   Sujay Kumar  Initial Specification
! 
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: MOD10A1_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: Mod10a1obs
!EOP
  type, public :: mod10a1obsdec
     character*100           :: odir
     integer                 :: modis_nc, modis_nr
     real, allocatable           :: rlat(:)
     real, allocatable           :: rlon(:)
     integer, allocatable        :: n11(:)
     integer, allocatable        :: n12(:)
     integer, allocatable        :: n21(:)
     integer, allocatable        :: n22(:)     
     real,    allocatable        :: w11(:)
     real,    allocatable        :: w12(:)
     real,    allocatable        :: w21(:)
     real,    allocatable        :: w22(:)
     logical                 :: startFlag
     type(proj_info)         :: mod_proj
  end type mod10a1obsdec

  type(mod10a1obsdec),save:: mod10a1obs(2)

contains
  
!BOP
! 
! !ROUTINE: MOD10A1_obsInit
! \label{MOD10A1_obsInit}
!
! !INTERFACE: 
  subroutine MOD10A1_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_logMod
    use LVT_histDataMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!
!   This subroutine initializes and sets up the data structures required
!   for reading the MOD10A1 data, including the setup of spatial interpolation
!   weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer              :: status
    real                 :: gridDesci(50)
    integer              :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real                 :: upgmt
    real                 :: cornerlat1, cornerlat2, cornerlon1, cornerlon2

    call ESMF_ConfigGetAttribute(LVT_Config, mod10a1obs(i)%odir, &
         label='MOD10A1 observation directory:', rc=status)
    call LVT_verify(status, 'MOD10A1 observation directory: not defined')

    call LVT_update_timestep(LVT_rc, 86400)
    allocate(mod10a1obs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(mod10a1obs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    allocate(mod10a1obs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(mod10a1obs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(mod10a1obs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(mod10a1obs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
    allocate(mod10a1obs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(mod10a1obs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(mod10a1obs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(mod10a1obs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))

    cornerlat1 = max(-59.995,nint((LVT_rc%gridDesc(4)+59.995)/0.01)*0.01-59.995-2*0.01)
    cornerlon1 = max(-179.995,nint((LVt_rc%gridDesc(5)+179.995)/0.01)*0.01-179.995-2*0.01)
    cornerlat2 = min(89.995,nint((LVT_rc%gridDesc(7)+59.995)/0.01)*0.01-59.995+2*0.01)
    cornerlon2 = min(179.995,nint((LVT_rc%gridDesc(8)+179.995)/0.01)*0.01-179.995+2*0.01)
    
    mod10a1obs(i)%modis_nr = nint((cornerlat2-cornerlat1)/0.01)+1
    mod10a1obs(i)%modis_nc = nint((cornerlon2-cornerlon1)/0.01)+1

    gridDesci = 0 
    gridDesci(1) = 0 
    gridDesci(2) = mod10a1obs(i)%modis_nc
    gridDesci(3) = mod10a1obs(i)%modis_nr
    gridDesci(4) = cornerlat1
    gridDesci(5) = cornerlon1
    gridDesci(7) = cornerlat2
    gridDesci(8) = cornerlon2
    gridDesci(6) = 128
    gridDesci(9) = 0.01
    gridDesci(10) = 0.01
    gridDesci(20) = 64


    call bilinear_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr, &
         mod10a1obs(i)%rlat, mod10a1obs(i)%rlon,&
         mod10a1obs(i)%n11, mod10a1obs(i)%n12, &
         mod10a1obs(i)%n21, mod10a1obs(i)%n22, & 
         mod10a1obs(i)%w11, mod10a1obs(i)%w12, &
         mod10a1obs(i)%w21, mod10a1obs(i)%w22)

    mod10a1obs(i)%startflag = .true. 

    call map_set(PROJ_LATLON, gridDesci(4), gridDesci(5), &
         0.0, gridDesci(9), gridDesci(10), 0.0, &
         mod10a1obs(i)%modis_nc, mod10a1obs(i)%modis_nr, mod10a1obs(i)%mod_proj)
  end subroutine MOD10A1_obsinit


end module MOD10A1_obsMod
