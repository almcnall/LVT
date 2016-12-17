!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: GIMMS_NDVIobsMod
! \label(GIMMS_NDVIobsMod)
!
! !INTERFACE:
module GIMMS_NDVIobsMod
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
!  7 Mar 2015   Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GIMMS_NDVIobsinit !Initializes structures for reading MOD16A2 data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GIMMSNDVIobs !Object to hold GIMMSNDVI observation attributes
!EOP

  type, public :: gimmsndvidec
     character*100           :: odir
     integer                 :: nc, nr
     real,    allocatable    :: rlat(:)
     real,    allocatable    :: rlon(:)
     integer, allocatable    :: n11(:)
     real                    :: gridDesc(50)
     integer                 :: yr
     integer                 :: mo
     logical                 :: startFlag
     real                    :: datares
  end type gimmsndvidec
     
  type(gimmsndvidec), save :: GIMMSNDVIObs(2)

contains
  
!BOP
! 
! !ROUTINE: GIMMS_NDVIobsInit
! \label{GIMMS_NDVIobsInit}
!
! !INTERFACE: 
  subroutine GIMMS_NDVIobsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_logMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine initializes and sets up the data structures required
!   for reading the GIMMS NDVI data, including the computation of spatial 
!   interpolation weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: status

    call ESMF_ConfigGetAttribute(LVT_Config, GIMMSNDVIobs(i)%odir, &
         label='GIMMS NDVI data directory:',rc=status)
    call LVT_verify(status, 'GIMMS NDVI data directory: not defined')

    call LVT_update_timestep(LVT_rc, 2592000)

    gimmsndviobs(i)%gridDesc = 0
        
    gimmsndviobs(i)%nc = 4320
    gimmsndviobs(i)%nr = 2160

    !filling the items needed by the interpolation library
    gimmsndviobs(i)%gridDesc(1) = 0  
    gimmsndviobs(i)%gridDesc(2) = gimmsndviobs(i)%nc
    gimmsndviobs(i)%gridDesc(3) = gimmsndviobs(i)%nr
    gimmsndviobs(i)%gridDesc(4) = -89.9583 ! dx/2,dy/2 = 1/24 = 0.04167
    gimmsndviobs(i)%gridDesc(5) = -179.9583
    gimmsndviobs(i)%gridDesc(7) = 89.9583
    gimmsndviobs(i)%gridDesc(8) = 179.9583
    gimmsndviobs(i)%gridDesc(6) = 128
    gimmsndviobs(i)%gridDesc(9) = 0.0833
    gimmsndviobs(i)%gridDesc(10) = 0.0833
    gimmsndviobs(i)%gridDesc(20) = 64

    gimmsndviobs(i)%datares  = 0.0833

    if(LVT_isAtAfinerResolution(gimmsndviobs(i)%datares)) then
       
       allocate(gimmsndviobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(gimmsndviobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(gimmsndviobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       
       call neighbor_interp_input(gimmsndviobs(i)%gridDesc, &
            LVT_rc%gridDesc, &
            LVT_rc%lnc*LVT_rc%lnr, &
            gimmsndviobs(i)%rlat, &
            gimmsndviobs(i)%rlon, &
            gimmsndviobs(i)%n11)
    else
       allocate(gimmsndviobs(i)%n11(gimmsndviobs(i)%nc*gimmsndviobs(i)%nr))
       call upscaleByAveraging_input(gimmsndviobs(i)%gridDesc,&
            LVT_rc%gridDesc,gimmsndviobs(i)%nc*gimmsndviobs(i)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,gimmsndviobs(i)%n11)
    endif

    gimmsndviobs(i)%yr = -1
    gimmsndviobs(i)%mo = LVT_rc%mo
    gimmsndviobs(i)%startFlag = .false. 

  end subroutine GIMMS_NDVIobsinit


end module GIMMS_NDVIobsMod
