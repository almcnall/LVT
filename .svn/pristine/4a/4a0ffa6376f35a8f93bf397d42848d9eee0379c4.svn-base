!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: FLUXNET_obsMod
! \label(FLUXNET_obsMod)
!
! !INTERFACE:
module FLUXNET_obsMod
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
  PUBLIC :: FLUXNET_obsinit !Initializes structures for reading FLUXNET data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: FLUXNETobs !Object to hold FLUXNET observation attributes
!EOP

  type, public :: fluxnetdec
     character*100           :: odir
     real, allocatable           :: rlat(:)
     real, allocatable           :: rlon(:)
     integer, allocatable        :: n11(:)
     real,    allocatable        :: qle(:,:)
     real,    allocatable        :: qh(:,:)
     integer                 :: yr
     integer                 :: mo
  end type fluxnetdec
     
  type(fluxnetdec), save :: FLUXNETObs(2)

contains
  
!BOP
! 
! !ROUTINE: FLUXNET_obsInit
! \label{FLUXNET_obsInit}
!
! !INTERFACE: 
  subroutine FLUXNET_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 

! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine initializes and sets up the data structures required
!   for reading the FLUXNET data, including the computation of spatial 
!   interpolation weights. The FLUXNET data is provides in the 
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

    call ESMF_ConfigGetAttribute(LVT_Config, FLUXNETobs(I)%odir, &
         label='FLUXNET data directory: ',rc=status)
    call LVT_verify(status, 'FLUXNET data directory: not defined')

    allocate(FLUXNETobs(I)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(FLUXNETobs(I)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    allocate(FLUXNETobs(I)%n11(LVT_rc%lnc*LVT_rc%lnr))

    allocate(FLUXNETobs(I)%qle(LVT_rc%lnc*LVT_rc%lnr,12))
    allocate(FLUXNETobs(I)%qh(LVT_rc%lnc*LVT_rc%lnr,12))

    gridDesci = 0
    
    !filling the items needed by the interpolation library
    gridDesci(1) = 0  !input is EASE grid
    gridDesci(2) = 720
    gridDesci(3) = 291
    gridDesci(4) = -55.25
    gridDesci(5) = -179.75
    gridDesci(7) = 89.75
    gridDesci(8) = 179.75
    gridDesci(6) = 128
    gridDesci(9) = 0.50
    gridDesci(10) = 0.50
    gridDesci(20) = 64
 
    call neighbor_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr, &
         FLUXNETobs(I)%rlat, FLUXNETobs(I)%rlon,&
         FLUXNETobs(I)%n11)

    FLUXNETobs(I)%yr = -1
    FLUXNETobs(I)%mo = LVT_rc%mo

    call LVT_update_timestep(LVT_rc, 2592000)

  end subroutine FLUXNET_obsinit


end module FLUXNET_obsMod
