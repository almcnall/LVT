!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !ROUTINE: readUSGSGWwellobs
! \label{readUSGSGWwellobs}
!
! !INTERFACE: 
subroutine readUSGSGWwellobs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,    only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_logMod,     only : LVT_logunit, LVT_getNextUnitNumber, &
       LVT_releaseUnitNumber, LVT_verify, LVT_endrun
  use LVT_timeMgrMod, only : LVT_calendar, LVT_tick
  use USGSGWwell_obsMod,    only : usgsgwwellobs
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
  integer,    intent(in) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for ground water
! well data. 
! LVT expects the data to be organized per calendar year, with 
! each file containing a daily data. Each reported observation is
! assumed to be time averaged. 
! 
! At the start of the each year, the whole year's data is read and
! stored. At other times, LVT simply indexes into the stored
! arrays for retrieving the required values. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  07 Nov 2013: Sujay Kumar, Initial Specification
! 
!EOP
  integer, parameter     :: CONST_VAL = 0.0
  integer                :: i,j,t,c,r,jj
  integer                :: stn_col, stn_row
  real                   :: col,row
  character*100          :: usgsgwwellfilename
  logical                :: file_exists
  logical                :: readflag
  integer                :: ftn, ios,ios1
  integer                :: yr,mo,da
  integer                :: status
  type(ESMF_Time)        :: usgsgwwelltime, usgsgwwelltime1
  integer                :: stnindex,tind
  real                   :: offset
  real                   :: gw_data, prcp_data
  character*100          :: line
  character*10           :: datestring
  integer                :: kk, iloc
  real                   :: gw(LVT_rc%lnc,LVT_rc%lnr)

  gw  = LVT_rc%udef

  if((LVT_rc%dyr(source).ne.usgsgwwellobs(source)%yr)) then 
     
     usgsgwwellobs(source)%yr = LVT_rc%dyr(source)
     usgsgwwellobs(source)%gw = LVT_rc%udef

     call ESMF_TimeSet(usgsgwwellobs(source)%startTime,  yy=usgsgwwellobs(source)%yr, &
          mm = LVT_rc%dmo(source), &
          dd = 1, &
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status, 'error in setting usgsgwwell start time')

     do i=1,usgsgwwellobs(source)%nstns
        call create_USGSGWwell_filename(usgsgwwellobs(source)%odir,  &
             usgsgwwellobs(source)%stnid(i),LVT_rc%dyr(source), usgsgwwellfilename)
        
        inquire(file=trim(usgsgwwellfilename),exist=file_exists)
        
        if(file_exists) then 
           write(LVT_logunit,*) '[INFO] Reading Well GW data ',trim(usgsgwwellfilename)
           ftn=LVT_getNextUnitNumber()
           open(ftn,file=trim(usgsgwwellfilename),form='formatted')
           
           readflag = .true. 
           do while(readflag) 
              read(ftn,'(a)',iostat=ios) line
              if(ios.ne.0) then 
                 readflag = .false. 
                 exit
              endif
              iloc = index(line,",")
              read(line(1:iloc-1),*) yr
              line = line(iloc+1:len(line))
              
              iloc = index(line,",")
              read(line(1:iloc-1),*) mo
              line = line(iloc+1:len(line))
              
              iloc = index(line,",")
              read(line(1:iloc-1),*) da
              line = line(iloc+1:len(line))
              
              read(line(1:iloc-1),*,iostat=ios1) gw_data
              if(gw_data.lt.0) gw_data = -9999.0
              if(ios1.ne.0) gw_data = -9999.0
              
              if(ios.ne.0) then 
                 readflag = .false. 
              else
                 call ESMF_TimeSet(usgsgwwelltime,yy=yr,&
                      mm=mo, dd=da, h=0,calendar=LVT_calendar,&
                      rc=status)
                 call LVT_verify(status,'ESMF_TimeSet in readUsgsgwwellobs(Source)')
                 
                 t = nint((usgsgwwelltime-usgsgwwellobs(source)%starttime)/&
                      usgsgwwellobs(source)%timestep)+1
                 
                 if(t.ge.1.and.t.le.366) then 
                    if(gw_data.ge.0.and.usgsgwwellobs(source)%stnqa(i).ge.2) then 
                       usgsgwwellobs(source)%gw(i,t) = CONST_VAL - gw_data*&
                            usgsgwwellobs(source)%stnsy(i)*0.3048 ! m to ft
                    else
                       usgsgwwellobs(source)%gw(i,t) = LVT_rc%udef
                    endif
                 endif
                 
              endif
              if(ios.ne.0) then 
                 readflag = .false. 
              endif
           enddo
           call LVT_releaseUnitNumber(ftn)
        endif
     enddo
  endif

  call ESMF_TimeSet(usgsgwwelltime1,yy=LVT_rc%dyr(source), &
       mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), h=LVT_rc%dhr(source), &
       m=LVT_rc%dmn(source), &
       s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
  call LVT_verify(status, 'usgsgwwelltime1 set failed')

  offset = (usgsgwwelltime1-usgsgwwellobs(source)%starttime)/&
       usgsgwwellobs(source)%timestep

  if((nint(offset)-offset).eq.0) then 
     tind = nint((usgsgwwelltime1-usgsgwwellobs(source)%starttime)/&
          usgsgwwellobs(source)%timestep)+1
     
     do i=1,usgsgwwellobs(source)%nstns
        call latlon_to_ij(LVT_domain%lvtproj, &
             usgsgwwellobs(source)%stnlat(i),usgsgwwellobs(source)%stnlon(i),&
             col,row)
        stn_col = nint(col)
        stn_row = nint(row)

        if(stn_col.gt.0.and.stn_row.gt.0.and.tind.ge.0.and.&
             stn_col.le.LVT_rc%lnc.and.stn_row.le.LVT_rc%lnr) then 
           if(usgsgwwellobs(source)%gw(i,tind).ne.LVT_rc%udef) then 
              gw(stn_col, stn_row) =&
                   usgsgwwellobs(source)%gw(i,tind)
           endif
        endif
     enddo
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_WATERTABLED,source, gw,vlevel=1,units="m")
 
end subroutine readUSGSGWwellobs

!BOP
! 
! !ROUTINE: create_USGSGWwell_filename
! \label(create_USGSGWwell_filename)
!
! !INTERFACE:
subroutine create_USGSGWwell_filename(odir, stnid, yr,fname)
! 
! !USES:   
  use LVT_String_Utility
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in)  :: odir
  character(len=*), intent(in)  :: stnid
  integer,          intent(in)  :: yr
! 
! !OUTPUT PARAMETERS:
  character(len=*), intent(out) :: fname
!
! !DESCRIPTION: 
! 
! This routine creates a filename for the ground water well
! 
!  The arguments are: 
!  \begin{description}
!   \item[stnid] Station ID 
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  character*4             :: fyr

  write(fyr, '(i4.4)' ) yr

  fname = trim(odir)//'/'//trim(stnid)//'_'//trim(fyr)//'.txt'
  
end subroutine create_USGSGWwell_filename
