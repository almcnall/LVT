!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: LISda_obsMod
! \label(LISda_obsMod)
!
! !INTERFACE:
module LISda_obsMod
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the use of a LIS model simulation output as 
!  "observations". 
!  
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LISda_obsInit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: lisdaobs
!
!EOP
  
  type, public :: daobsdec
     character*100  :: odir
     integer        :: scal
  end type daobsdec

  type(daobsdec)  :: lisdaobs(2)

contains

!BOP
! 
! !ROUTINE: LISda_obsInit
! \label{LISda_obsInit}
!
! !INTERFACE: 
  subroutine LISda_obsInit(i)
! 
! !USES: 
    use ESMF
    use LVT_coreMod,    only : LVT_rc, LVT_config
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod,     only : LVT_verify

    implicit none
!
! !INPUT PARAMETERS: 
    integer,     intent(IN) :: i   ! index of the observation type
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! This routine initializes the structures required for the handling of a 
! land surface model output (from a LIS simulation) as observations.  
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    real                    :: run_dd(8)
    integer                 :: t
    integer                 :: ts
    type(ESMF_Config)       :: modelSpecConfig
    character*20            :: domain
    character*10            :: time
    integer                 :: rc

    call ESMF_ConfigGetAttribute(LVT_config,lisdaobs(i)%odir, &
         label="LIS DAOBS output directory:",rc=rc)
    call LVT_verify(rc,'LIS DAOBS output directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,time,&
         label="LIS DAOBS output interval:",rc=rc)
    call LVT_verify(rc,'LIS DAOBS output interval: not defined')

    call LVT_parseTimeString(time,ts)

    call LVT_update_timestep(LVT_rc, ts)

    call ESMF_ConfigGetAttribute(LVT_config,lisdaobs(i)%scal, &
         label="LIS DAOBS use scaled obs:",rc=rc)
    call LVT_verify(rc,'LIS DAOBS use scaled obs: not defined')
    
  end subroutine LISda_obsInit
  
end module LISda_obsMod
