!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readGRACEObs
! \label{readGRACEObs}
!
! !INTERFACE: 
subroutine readGRACEObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,     only : LVT_rc
  use LVT_logMod
  use LVT_histDataMod
  use LVT_timeMgrMod
  use GRACE_obsMod, only : GRACEObs

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  NOTES: 
!   The GRACE output is available at monthly intervals. So 
!   the comparisons against model data should use at least a 
!   24 hour (1day) averaging interval. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  10 Dec 2010: Sujay Kumar, Initial Specification
! 
!EOP

  character*100           :: filename
  logical                 :: file_exists
  integer                 :: c,r,t,kk
  integer                 :: ftn
  integer                 :: iret
  real                    :: tws(LVT_rc%lnc,LVT_rc%lnr)
  real,  allocatable      :: tws_in(:)
  logical*1, allocatable  :: bitmap_in(:)
  logical*1               :: bitmap_out(LVT_rc%lnc*LVT_rc%lnr)
  type(ESMF_Time)         :: cTime
  type(ESMF_TimeInterval) :: tw
  integer                 :: yr,mo,da,hr,mn,ss
  integer                 :: status

  tws = LVT_rc%udef

  if(GRACEobs(source)%useRawData.eq.1) then 
     if(GRACEobs(source)%useRealSensorData.eq.1) then 
        
     else
        if(GRACEobs(source)%mo.ne.LVT_rc%d_nmo(source).and.&
             LVT_rc%dda(source).eq.1) then 
!             LVT_rc%dda(source).eq.15) then 

           if(GRACEobs(source)%startFlag) then 
              GRACEobs(source)%yr = LVT_rc%dyr(source)
              GRACEobs(source)%mo = LVT_rc%dmo(source)
              GRACEobs(source)%startFlag = .false. 
           endif
           
           call create_GRACE_Raw_filename(Graceobs(source)%odir, &
                GRACEobs(source)%yr,&
                GRACEobs(source)%mo, filename)

           inquire(file=trim(filename),exist=file_exists) 
           
           if(file_exists) then 
              write(LVT_logunit,*) '[INFO] Reading GRACE file ',trim(filename)
              
              allocate(tws_in(GRACEObs(source)%nc*GRACEobs(source)%nr))
              allocate(bitmap_in(GRACEObs(source)%nc*GRACEobs(source)%nr))
              
              ftn = LVT_getNextUnitNumber()
              open(ftn,file=trim(filename),form='unformatted')
              read(ftn) tws_in
              call LVT_releaseUnitNumber(ftn)
              
              bitmap_in = .true. 
              do t=1,GRACEObs(source)%nc*GRACEobs(source)%nr
                 if(tws_in(t).eq.-9999.0) then 
                    bitmap_in(t) = .false.
                 endif
              enddo
              
              call bilinear_interp(LVT_rc%gridDesc,bitmap_in,tws_in,    &
                   bitmap_out,tws,GRACEobs(source)%nc*GRACEobs(source)%nr,&
                   LVT_rc%lnc*LVT_rc%lnr,   &
                   GRACEobs(source)%rlat,GRACEobs(source)%rlon, &
                   GRACEobs(source)%w11,GRACEobs(source)%w12,   &
                   GRACEobs(source)%w21,GRACEobs(source)%w22,   &
                   GRACEobs(source)%n11,GRACEobs(source)%n12,   &
                   GRACEobs(source)%n21,GRACEobs(source)%n22,   &
                   LVT_rc%udef,iret)

           else
              tws = LVT_rc%udef
           endif
           
           GRACEobs(source)%yr = LVT_rc%d_nyr(source)
           GRACEobs(source)%mo = LVT_rc%d_nmo(source)
                      
        else
           tws  = LVT_rc%udef
        endif

     endif
  else
     if(GRACEobs(source)%config.eq."default".or.&
          GRACEobs(source)%config.eq."follow-on") then 
        if(GRACEobs(source)%mo.ne.LVT_rc%d_nmo(source).and.&
             LVT_rc%dda(source).eq.1) then 
!             LVT_rc%dda(source).eq.15) then 
           
           if(GRACEobs(source)%startFlag) then 
              GRACEobs(source)%yr = LVT_rc%dyr(source)
              GRACEobs(source)%mo = LVT_rc%dmo(source)
              GRACEobs(source)%startFlag = .false. 
           endif
           
           call create_GRACE1_filename(Graceobs(source)%odir, &
                GRACEobs(source)%yr,&
                GRACEobs(source)%mo, filename)
           
           inquire(file=trim(filename),exist=file_exists) 
           
           if(file_exists) then 
              write(LVT_logunit,*) '[INFO] Reading GRACE file ',trim(filename)
              
              ftn = LVT_getNextUnitNumber()
              open(ftn,file=trim(filename),form='unformatted')
              read(ftn) tws
              call LVT_releaseUnitNumber(ftn)
              
           else
              tws = LVT_rc%udef
           endif
           
           GRACEobs(source)%yr = LVT_rc%d_nyr(source)
           GRACEobs(source)%mo = LVT_rc%d_nmo(source)
           
        else
           tws  = LVT_rc%udef
        endif
     else
        call ESMF_TimeSet(cTime, yy = LVT_rc%dyr(source), &
             mm = LVT_rc%dmo(source), &
             dd = LVT_rc%dda(source), &
             h  = LVT_rc%dhr(source), &
             m  = LVT_rc%dmn(source), & 
             s  = LVT_rc%dss(source), &
             calendar = LVT_calendar, & 
             rc = status)
        call LVT_verify(status)
        !7 day
        call ESMF_TimeIntervalSet(tw,s=604800,rc=status)
        ctime = ctime - tw
        
        call ESMF_TimeGet(ctime,yy=yr,mm=mo,dd=da,calendar=LVT_calendar,&
             rc=status)
        call LVT_verify(status)
        
        if(GRACEobs(source)%da.ne.da) then 
           if(GRACEobs(source)%startFlag) then 
              GRACEobs(source)%yr = yr
              GRACEobs(source)%mo = mo
              GRACEobs(source)%da = da
              GRACEobs(source)%startFlag = .false. 
           endif
           
           call create_GRACE2_filename(Graceobs(source)%odir, yr, &
                mo, da, filename)
           
           inquire(file=trim(filename),exist=file_exists) 
           
           if(file_exists) then 
              write(LVT_logunit,*) '[INFO] Reading GRACE file ',trim(filename)
              
              ftn = LVT_getNextUnitNumber()
              open(ftn,file=trim(filename),form='unformatted')
              read(ftn) tws
              call LVT_releaseUnitNumber(ftn)
              
           else
              tws = LVT_rc%udef
           endif
           GRACEobs(source)%yr = yr
           GRACEobs(source)%mo = mo
           GRACEobs(source)%da = da
           
        endif
     endif
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_TWS,source,tws,vlevel=1,units="mm")
  
end subroutine readGRACEObs

!BOP
! 
! !ROUTINE: create_GRACE1_filename
! \label{create_GRACE1_filename}
!
! !INTERFACE: 
subroutine create_GRACE1_filename(odir,yr,mo,filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for GRACE_LH data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GRACE base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the GRACE_LH file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr,mo
  character(len=*)             :: filename
!EOP

  character*6             :: fdate
  
  write(unit=fdate, fmt='(i4.4,i2.2)') yr,mo
  
  filename = trim(odir)//'/GRACE_obs_'//trim(fdate)//'.bin'
  
end subroutine create_GRACE1_filename


!BOP
! 
! !ROUTINE: create_GRACE2_filename
! \label{create_GRACE2_filename}
!
! !INTERFACE: 
subroutine create_GRACE2_filename(odir,yr,mo,da, filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for GRACE_LH data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GRACE base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the GRACE_LH file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr,mo,da
  character(len=*)             :: filename
!EOP

  character*8           :: fdate
  
  write(unit=fdate, fmt='(i4.4,i2.2,i2.2)') yr,mo,da
  
  filename = trim(odir)//'/GRACE_obs_'//trim(fdate)//'.bin'
  
end subroutine create_GRACE2_filename


!BOP
! 
! !ROUTINE: create_GRACE_raw_filename
! \label{create_GRACE_raw_filename}
!
! !INTERFACE: 
subroutine create_GRACE_raw_filename(odir,yr,mo,filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for GRACE_LH data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GRACE base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the GRACE_LH file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr,mo
  character(len=*)             :: filename
!EOP

  character*6             :: fdate
  
  write(unit=fdate, fmt='(i4.4,i2.2)') yr,mo
  
  filename = trim(odir)//'/simGRACE_'//trim(fdate)//'0112.bin'
  
end subroutine create_GRACE_raw_filename
