!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readLISdaAsObs
! \label(readLISdaAsObs)
!
! !INTERFACE:
subroutine readLISdaAsObs(source)
! !USES:   
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif
  use LVT_coreMod,      only : LVT_rc
  use LVT_histDataMod
  use LVT_logMod
  use LISda_obsMod,    only : lisdaobs

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in) :: source 
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
  integer          :: c,r

  character*100      :: cdate, cdate1

   
   write(unit=cdate, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
        LVT_rc%dyr(source), LVT_rc%dmo(source), &
        LVT_rc%dda(source), LVT_rc%dhr(source),LVT_rc%dmn(source)

   write(unit=cdate1, fmt='(i4.4, i2.2)') &
        LVT_rc%dyr(source), LVT_rc%dmo(source)
     
   fname = trim(lisdaobs(source)%odir)//'/' & 
        //trim(cdate1)//'/'&
        //trim(cdate)//'.1gs4r'
     
  inquire(file=trim(fname),exist=file_exists)
  
  write(LVT_logunit,*) '[INFO] reading DA obs output ',trim(fname)
  if(file_exists) then
     ftn = 11
     open(ftn,file=trim(fname), form='unformatted')
     read(ftn) obsData
     if(lisdaobs(source)%scal.eq.1) then 
        read(ftn) obsData
     endif
     close(ftn)

  else
     write(LVT_logunit,*) '[WARN] Warning: DAobs file ',trim(fname),' does not exist'
     obsData = -9999.0
  endif
  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source, obsData,vlevel=1,units="m3/m3")

end subroutine readLISdaAsObs

