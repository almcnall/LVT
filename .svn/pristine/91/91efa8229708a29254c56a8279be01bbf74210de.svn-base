!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: GL6JULES_obsMod
! \label(GL6JULES_obsMod)
!
! !INTERFACE:
module GL6JULES_obsMod
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
!  8 July 2015   Sujay Kumar  Initial Specification
! 
!EOP
!

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GL6JULES_obsinit !Initializes structures for reading GL6JULES data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GL6JULESobs !Object to hold GL6JULES observation attributes
!EOP

  type, public :: fluxnetdec
     character*100               :: odir
     real,    allocatable        :: qle_c(:,:,:)
     real,    allocatable        :: time_val(:)
     type(ESMF_Time)             :: refTime
     integer                     :: ntimes
     integer                     :: yr
  end type fluxnetdec
     
  type(fluxnetdec), save :: GL6JULESObs(2)

contains
  
!BOP
! 
! !ROUTINE: GL6JULES_obsInit
! \label{GL6JULES_obsInit}
!
! !INTERFACE: 
  subroutine GL6JULES_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
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
!   for reading the GL6JULES data, including the computation of spatial 
!   interpolation weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: i 
    integer               :: status

    call ESMF_ConfigGetAttribute(LVT_Config, GL6JULESobs(i)%odir, &
         label='GL6 JULES data directory: ',rc=status)
    call LVT_verify(status, 'GL6 JULES data directory: not defined')

    call ESMF_TimeSet(GL6JULESobs(i)%refTime,  yy=1980, &
         mm = 1, &
         dd = 1,&
         h = 0,&
         m = 0,&
         s = 0,&
         calendar = LVT_calendar, &
         rc=status)
    call LVT_verify(status, 'ESMF_TimeSet error in GL6JULES_obsInit')

    GL6JULESobs(i)%yr = -1

  end subroutine GL6JULES_obsinit

end module GL6JULES_obsMod
