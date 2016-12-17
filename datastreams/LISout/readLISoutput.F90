!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readLISoutput
! \label(readLISoutput)
!
! !INTERFACE:
subroutine readLISoutput(source)
! !USES:   
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_fileIOMod
  use LVT_statsDataMod
  use LVT_LISoutputHandlerMod

  implicit none
  
  integer, intent(in)    :: source
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
! 
!EOP
  
  character*100    :: fname 
  logical          :: file_exists
  real             :: obsData(LVT_rc%lnc, LVT_rc%lnr)
  
  integer          :: t
  integer          :: ftn
  type(LVT_metadataEntry), pointer :: obs
  integer          :: c,r
    
  if(LVT_rc%chkTS) then 
     if(LVT_LIS_rc(source)%anlys_data_class.eq."LSM") then 
        call LVT_create_output_filename(LVT_LIS_rc(source)%nest, source,fname, &
             'SURFACEMODEL',&
             writeint = LVT_rc%ts, wout = LVT_LIS_rc(source)%format,&
             style=LVT_LIS_rc(source)%style, &
             odir=LVT_LIS_rc(source)%odir)
     elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Routing") then 
        call LVT_create_output_filename(LVT_LIS_rc(source)%nest, source, fname, &
             'ROUTING',&
             writeint = LVT_rc%ts, wout = LVT_LIS_rc(source)%format,&
             style=LVT_LIS_rc(source)%style,&
             odir=LVT_LIS_rc(source)%odir)
     elseif(LVT_LIS_rc(source)%anlys_data_class.eq."RTM") then 
        call LVT_create_output_filename(LVT_LIS_rc(source)%nest, source, fname, &
             'RTM',&
             writeint = LVT_rc%ts, wout = LVT_LIS_rc(source)%format,&
             style=LVT_LIS_rc(source)%style, &
             odir=LVT_LIS_rc(source)%odir)
     elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Irrigation") then 
        call LVT_create_output_filename(LVT_LIS_rc(source)%nest, source, fname, &
             'IRRIGATION',&
             writeint = LVT_rc%ts, wout = LVT_LIS_rc(source)%format,&
             style=LVT_LIS_rc(source)%style,&
             odir=LVT_LIS_rc(source)%odir)
     endif
     inquire(file=trim(fname),exist=file_exists)
       
     if(file_exists) then 
        ! The following line is modified by Shugong Wang: LVT should be LIS
        !write(LVT_logunit,*) 'Reading LVT output ',trim(fname)
        write(LVT_logunit,*) '[INFO] Reading LIS output ',trim(fname)
        call LVT_readLISModelOutput(trim(fname),source, LVT_LIS_rc(source)%format,&
             wopt=LVT_LIS_rc(source)%wopt)
     else
!        LVT_stats%datamask = 0.0  !initially set all to false. 
!toggle these print/stop statements on/off
!if LVT needs to work with missing LIS files
        write(LVT_logunit,*) '[WARN] LIS file ',trim(fname),' does not exist'
!          write(LVT_logunit,*) 'Program stopping.. '
!          stop
     endif
  endif

end subroutine readLISoutput

