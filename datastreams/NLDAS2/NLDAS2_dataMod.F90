!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
!
! !MODULE: NLDAS2_dataMod
!
! !INTERFACE:
!
! !USES:
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!   This subroutine provides the observation plugin for reading the
!   operational NLDAS-2 output.
!
! !FILES USED:
!
! !REVISION HISTORY:
!  09 Dec 2010  Sujay Kumar, Initial Specification
!  27 Jan 2014: David Mocko, Updates for NLDAS-2 SAC
!  11 Dec 2014: David Mocko, Added additional NLDAS-2 variables
!
!EOP
module NLDAS2_dataMod
! !USES:
  use ESMF
  
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: NLDAS2_datainit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: NLDAS2data
!EOP
  type, public            :: nldas2datadec
     character*120           :: odir
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
     real,    allocatable    :: vic_depth1(:,:)
     real,    allocatable    :: vic_depth2(:,:)
     real,    allocatable    :: vic_depth3(:,:)
     integer                 :: nc
     integer                 :: nr
     character*10            :: lsm
     logical                 :: smvol
     type(ESMF_TimeInterval) :: ts
     character*50            :: anlys_data_class
  end type nldas2datadec

  type(nldas2datadec),save:: nldas2data(2)

contains
!BOP
!
! !ROUTINE: NLDAS2_dataInit
! \label{NLDAS2_dataInit}
!
! !INTERFACE:
  subroutine NLDAS2_datainit(i)
!
! !USES:
    use LVT_coreMod
    use LVT_logMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    
    implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!   This subroutine initializes and sets up the data structures required
!   for reading the NLDAS-2 data, including the setup of spatial interpolation
!   weights.
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP
!BOP
! !ARGUMENTS:
    integer,  intent(IN) :: i
!EOP
    integer              :: status
    real                 :: gridDesci(50)
    integer              :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real                 :: upgmt
    integer              :: ftn
    character*50         :: vic_d1file,vic_d2file,vic_d3file
    
    call ESMF_ConfigGetAttribute(LVT_Config, nldas2data(i)%odir,        &
         label='NLDAS2 data directory:', rc=status)
    call LVT_verify(status, 'NLDAS2 data directory: not defined')
    
    call ESMF_ConfigGetAttribute(LVT_Config, nldas2data(i)%lsm,         &
         label='NLDAS2 land surface model:', rc=status)
    call LVT_verify(status, 'NLDAS2 land surface model: not defined')
    
    call ESMF_ConfigGetAttribute(LVT_Config, nldas2data(i)%smvol,       &
         label='NLDAS2 soil moisture volumetric:', rc=status)
    call LVT_verify(status,                                          &
         'NLDAS2 soil moisture volumetric: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, nldas2data(i)%anlys_data_class, &
         label='NLDAS2 analysis data class:', rc=status)
    call LVT_verify(status,                                          &
         'NLDAS2 analysis data class: not defined')
    
    call LVT_update_timestep(LVT_rc, 3600)
    call ESMF_TimeIntervalSet(nldas2data(i)%ts,s=3600,rc=status)
    call LVT_verify(status)
    
    allocate(nldas2data(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(nldas2data(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
    allocate(nldas2data(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(nldas2data(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(nldas2data(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(nldas2data(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
    allocate(nldas2data(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(nldas2data(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(nldas2data(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(nldas2data(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
    
    gridDesci    = 0
    gridDesci(1) = 0
    gridDesci(2) = 464
    gridDesci(3) = 224
    gridDesci(4) = 25.0625
    gridDesci(5) = -124.9375
    gridDesci(7) = 52.9375
    gridDesci(8) = -67.0625
    gridDesci(6) = 128
    gridDesci(9) = 0.125
    gridDesci(10) = 0.125
    gridDesci(20) = 64
    
    nldas2data(i)%nc = 464
    nldas2data(i)%nr = 224
    
    call bilinear_interp_input(gridDesci,LVT_rc%gridDesc,            &
         LVT_rc%lnc*LVT_rc%lnr,                &
         nldas2data(i)%rlat, nldas2data(i)%rlon,     &
         nldas2data(i)%n11, nldas2data(i)%n12,       &
         nldas2data(i)%n21, nldas2data(i)%n22,       &
         nldas2data(i)%w11, nldas2data(i)%w12,       &
         nldas2data(i)%w21, nldas2data(i)%w22)

    if (nldas2data(i)%lsm.eq."VIC") then 
       call ESMF_ConfigGetAttribute(LVT_Config, vic_d1file,          &
            label='NLDAS2 VIC soil depth1 file:', rc=status)
       call LVT_verify(status,                                       &
            'NLDAS2 VIC soil depth1 file: not defined')

       call ESMF_ConfigGetAttribute(LVT_Config, vic_d2file,          &
            label='NLDAS2 VIC soil depth2 file:', rc=status)
       call LVT_verify(status,                                       &
            'NLDAS2 VIC soil depth2 file: not defined')

       call ESMF_ConfigGetAttribute(LVT_Config, vic_d3file,          &
            label='NLDAS2 VIC soil depth3 file:', rc=status)
       call LVT_verify(status,                                       &
            'NLDAS2 VIC soil depth3 file: not defined')
       
       allocate(nldas2data(i)%vic_depth1(nldas2data(i)%nc,nldas2data(i)%nr))
       allocate(nldas2data(i)%vic_depth2(nldas2data(i)%nc,nldas2data(i)%nr))
       allocate(nldas2data(i)%vic_depth3(nldas2data(i)%nc,nldas2data(i)%nr))
       
       ftn = LVT_getNextUnitNumber()
       open(ftn,file=vic_d1file,form='unformatted',access='direct',  &
            recl=nldas2data(i)%nc*nldas2data(i)%nr*4)
       read(ftn,rec=1) nldas2data(i)%vic_depth1
       call LVT_releaseUnitNumber(ftn)
       ftn = LVT_getNextUnitNumber()
       open(ftn,file=vic_d2file,form='unformatted',access='direct',  &
            recl=nldas2data(i)%nc*nldas2data(i)%nr*4)
       read(ftn,rec=1) nldas2data(i)%vic_depth2
       call LVT_releaseUnitNumber(ftn)
       ftn = LVT_getNextUnitNumber()
       open(ftn,file=vic_d3file,form='unformatted',access='direct',  &
            recl=nldas2data(i)%nc*nldas2data(i)%nr*4)
       read(ftn,rec=1) nldas2data(i)%vic_depth3
       call LVT_releaseUnitNumber(ftn)
    endif

  end subroutine NLDAS2_datainit
  
end module NLDAS2_dataMod
