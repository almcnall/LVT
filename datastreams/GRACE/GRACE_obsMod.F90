!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: GRACE_obsMod
! \label(GRACE_obsMod)
!
! !INTERFACE:
module GRACE_obsMod
! 
! !USES: 
  use ESMF

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  18 May 2011   Sujay Kumar  Initial Specification
! 
!EOP
!

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GRACE_obsinit !Initializes structures for reading GRACE data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GRACEobs !Object to hold GRACE observation attributes
!EOP

  type, public :: fluxnetdec
     character*100           :: odir
     character*50            :: config
     logical                 :: startFlag
     integer                 :: useRawData
     integer                 :: useRealSensorData
     integer                 :: yr
     integer                 :: mo
     integer                 :: da

     integer                 :: nc
     integer                 :: nr
     real, allocatable       :: rlat(:)
     real, allocatable       :: rlon(:)
     integer, allocatable    :: n11(:)
     integer, allocatable    :: n12(:)
     integer, allocatable    :: n21(:)
     integer, allocatable    :: n22(:)     
     real,    allocatable    :: w11(:)
     real,    allocatable    :: w12(:)
     real,    allocatable    :: w21(:)
     real,    allocatable    :: w22(:)
  end type fluxnetdec
     
  type(fluxnetdec), save :: GRACEObs(2)

contains
  
!BOP
! 
! !ROUTINE: GRACE_obsInit
! \label{GRACE_obsInit}
!
! !INTERFACE: 
  subroutine GRACE_obsinit(i)
! 
! !USES: 
    use LVT_coreMod,   only : LVT_rc, LVT_Config
    use LVT_histDataMod
    use LVT_logMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 

! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine initializes and sets up the data structures required
!   for reading the GRACE data, including the computation of spatial 
!   interpolation weights. The GRACE data is provides in the 
!   EASE grid projection. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: i
    integer               :: status
    real                  :: gridDesci(50)

    call ESMF_ConfigGetAttribute(LVT_Config, GRACEObs(i)%odir, &
         label='GRACE data directory: ',rc=status)
    call LVT_verify(status, 'GRACE data directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, GRACEObs(i)%config, &
         label='GRACE configuration: ',rc=status)
    call LVT_verify(status, 'GRACE configuration: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, GRACEObs(i)%useRawData, &
         label='GRACE process raw (anomaly) data: ',rc=status)
    call LVT_verify(status, 'GRACE process raw (anomaly) data: not defined')

    if(GRACEobs(i)%userawdata.gt.0) then 
       call ESMF_ConfigGetAttribute(LVT_Config, GRACEObs(i)%useRealSensorData, &
            label='GRACE use data from real sensor: ',rc=status)
       call LVT_verify(status, 'GRACE use data from real sensor: not defined')


       allocate(GRACEObs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(GRACEObs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(GRACEObs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(GRACEObs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(GRACEObs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(GRACEObs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
       allocate(GRACEObs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(GRACEObs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(GRACEObs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(GRACEObs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
       
       gridDesci    = 0
       gridDesci(1) = 0
       gridDesci(2) = 360
       gridDesci(3) = 180
       gridDesci(4) = -89.5
       gridDesci(5) = -179.5
       gridDesci(7) = 89.5
       gridDesci(8) = 179.5
       gridDesci(6) = 128
       gridDesci(9) = 1.0
       gridDesci(10) = 1.0
       gridDesci(20) = 64
       
       GRACEObs(i)%nc = 360
       GRACEObs(i)%nr = 180
       
       call bilinear_interp_input(gridDesci,LVT_rc%gridDesc,            &
            LVT_rc%lnc*LVT_rc%lnr,                &
            GRACEObs(i)%rlat, GRACEObs(i)%rlon,     &
            GRACEObs(i)%n11, GRACEObs(i)%n12,       &
            GRACEObs(i)%n21, GRACEObs(i)%n22,       &
            GRACEObs(i)%w11, GRACEObs(i)%w12,       &
            GRACEObs(i)%w21, GRACEObs(i)%w22)
       
    endif

    GRACEobs(i)%yr = -1
!    GRACEobs(i)%mo = LVT_rc%mo
    GRACEobs(i)%mo = -1
    GRACEobs(i)%da = LVT_rc%da
    GRACEobs(i)%startFlag = .true. 

!    if(LVT_rc%tavgInterval.lt.2592000) then 
!       write(LVT_logunit,*) 'The time averaging interval must be greater than'
!       write(LVT_logunit,*) 'equal to a month since the GRACE data is monthly'
!       call LVT_endrun()
!    endif

    call LVT_update_timestep(LVT_rc, 86400)

  end subroutine GRACE_obsinit


end module GRACE_obsMod
