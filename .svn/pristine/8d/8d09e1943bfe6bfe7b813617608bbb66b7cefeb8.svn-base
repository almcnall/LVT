!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_ensPercentileMod
! \label(LVT_ensPercentileMod)
!
! !INTERFACE:
module LVT_ensPercentileMod
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_timeMgrMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
  use LVT_TSMod
  use LVT_logMod
  use LVT_CIMod
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the computation of the percentiles for a given variable
!  against a corresponding climatology. 
!  
!  The code makes two passes through the data. During the first pass, it
!  generates the percentile climatology and during the second pass, it computes
!  daily percentiles using the percentile climatology. Note that code 
!  requires the use of a time averaging interval of a day. 
! 
!  The general methodology for assembling percentile climatology is as follows:
!  For each grid point, the variable values for the calendar day 
!  across all years is assembled. For example, assume that a 10 year simulation
!  is used for the analysis (say 2000-2010) and we are computing percentiles 
!  for soil moisture. For day 1, 10 soil moisture values are assembled by 
!  using the Jan 1 values for year 2000, 2001, 2002,... to 2010. Similarly 
!  10 soil moisture values are assembled for 365 days. 
!
!  Once the percentile climatology is assembled, this array is then sorted 
!  in the ascending order. During the second pass through the data, each 
!  day's soil moisture value is then ranked against the percentile climatology. 
!  For example, soil moisture value from Jan 1, 2000 is ranked against
!  the 10 sorted percentile climatology values and the percentile is determined.
! 
!  The particular implementation used here follows an extension of the 
!  above-mentioned strategy, where a moving window of 5 days is employed 
!  to improve the sampling density while calculating percentiles. While
!  assembling the percentile climatology, instead of using a single day 
!  across all years, we use 5 days (2 previous days, current day and 2 
!  next days) from each year. For the 10 year example, we will have 
!  5 days for Jan 3 2000, using Jan 1 to Jan 5 values. As a result 
!  across the 10 years, each day will have 5x10 = 50 values. This sample
!  is then used to compute the percentiles as before. 
!
!
!  NOTE: Because of the long temporal dimension, the memory becomes limiting
!  to do these calculations all in memory. LVT gets around this issue by 
!  performing intermediate I/O.   
!
! !FILES USED:
!
! !REVISION HISTORY: 
!  22 Mar 2012    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initensPercentile
  public :: LVT_diagnoseensPercentile
  public :: LVT_computeensPercentile
  public :: LVT_writeMetric_ensPercentile
  public :: LVT_resetMetric_ensPercentile
  public :: LVT_writerestart_ensPercentile
  public :: LVT_readrestart_ensPercentile

!EOP

  type, public :: pctiledec
     integer          :: climo_calc_only
     integer          :: nyears
     integer          :: nwindow
     integer          :: nsize_total
     real             :: d_classes(5)
     integer, allocatable :: ftn_ts_darea(:)
  end type pctiledec

  type(pctiledec) :: LVT_ensPctile_struc
  private



contains
!BOP
!
! !ROUTINE: LVT_initensPercentile
! \label{LVT_initensPercentile}
! 
! !INTERFACE: 
  subroutine LVT_initensPercentile(model, obs, stats,metric)
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
!  This routine initializes the data structures required
!  for the ensPercentile computations. 
! 
!EOP

    type(ESMF_Time)         :: startTime
    type(ESMF_Time)         :: stopTime
    type(ESMF_TimeInterval) :: timeStep
    integer                 :: i,k,t
    character*100           :: filename
    character*3             :: ftime
    integer                 :: rc
    integer                 :: ftn
    real, allocatable         :: pctile_model(:,:)

    if(LVT_metrics%ensPercentile%selectOpt.eq.1.or.&
         LVT_metrics%ensPercentile%timeOpt.eq.1) then 

       LVT_ensPctile_struc%d_classes(1) = 0.30
       LVT_ensPctile_struc%d_classes(2) = 0.20
       LVT_ensPctile_struc%d_classes(3) = 0.10
       LVT_ensPctile_struc%d_classes(4) = 0.05
       LVT_ensPctile_struc%d_classes(5) = 0.02

       if(LVT_rc%tavgInterval.ne.86400) then
          write(LVT_logunit,*) 'EnsPercentile computations are only supported'
          write(LVT_logunit,*) 'with a daily average interval...'
          call LVT_endrun()
       endif
       LVT_ensPctile_struc%nwindow = 5
       call LVT_computeTimeSpanInYears(LVT_ensPctile_struc%nyears)

       LVT_ensPctile_struc%nsize_total = LVT_ensPctile_struc%nwindow * & 
            LVT_ensPctile_struc%nyears
       
!       allocate(stats%pctile_model_final(365,LVT_LIS_rc(1)%ntiles,&
!            LVT_ensPctile_struc%nsize_total))
!open 365 files for each day 
       allocate(pctile_model(LVT_LIS_rc(1)%ntiles,&
            LVT_ensPctile_struc%nsize_total))
       pctile_model = LVT_rc%udef

       if(LVT_rc%startmode.eq."coldstart") then 
          do k=1,365
             write(unit=ftime,fmt='(I3.3)') k
             filename = trim(LVT_rc%statsodir)//'/'//trim(ftime)//&
                  '_pctile_climo.bin'
             ftn = LVT_getNextUnitNumber()
             open(ftn,file=trim(filename),form='unformatted')
             write(ftn) pctile_model
             call LVT_releaseUnitNumber(ftn)
          enddo
       else
          do k=1,365
             write(unit=ftime,fmt='(I3.3)') k
             filename = trim(LVT_rc%statsodir)//'/'//&
                  trim(ftime)//'_pctile_climo.bin'
             ftn = LVT_getNextUnitNumber()
             open(ftn,file=trim(filename),form='unformatted')
             read(ftn) pctile_model
             call LVT_releaseUnitNumber(ftn)             
             
             do t=1,LVT_LIS_rc(1)%ntiles
                call compute_sorted_sizes(&
                     LVT_ensPctile_struc%nsize_total,&
                     pctile_model(t,:),&
                     stats%enspctile%value_model_nsize(t,k))
             enddo
          enddo
       endif

       deallocate(pctile_model)
       
       allocate(stats%enspctile%value_model_nsize(LVT_LIS_rc(1)%ntiles, &
            365))
       stats%enspctile%value_model_nsize = 0 
       
       allocate(stats%enspctile%value_model(LVT_LIS_rc(1)%ntiles,model%selectNlevs,1))
       stats%enspctile%value_model  = 0
       allocate(stats%enspctile%value_count_model(LVT_LIS_rc(1)%ntiles,model%selectNlevs,1))
       stats%enspctile%value_count_model  = 0

       allocate(stats%enspctile%value_model_ci(model%selectNlevs, 1))

!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------
       
       metric%npass = 2   
       metric%obsData = .true. 
       
       allocate(LVT_ensPctile_struc%ftn_ts_darea(LVT_rc%ntslocs))
       
       call system('mkdir -p '//(LVT_rc%statsodir))       

       do i=1,LVT_rc%ntslocs
          
          filename = trim(LVT_rc%statsodir)//'/'//&
               'PercentArea_ensPercentile_'//&
               trim(LVT_TSobj(i)%tslocname)//&
               '.dat'

          LVT_ensPctile_struc%ftn_ts_darea(i) = LVT_getNextUnitNumber()

          open(LVT_ensPctile_struc%ftn_ts_darea(i),file=(filename),&
               form='formatted')

       enddo

!-------------------------------------------------------------------------
! LVT can be used to compute/establish the climatology alone. If
! LVT is used to compute ensPercentiles based on an established climo, 
! then the 'restart' mode must be used. 
!-------------------------------------------------------------------------

       call ESMF_ConfigGetAttribute(LVT_config,&
            LVT_ensPctile_struc%climo_calc_only,&
            label="Compute only the climatology for percentiles:",rc=rc)
       call LVT_verify(rc,'Compute only the climatology for percentiles: option not specified in the config file')
       if(LVT_ensPctile_struc%climo_calc_only.eq.1) then 
          metric%npass = 1
       endif
       metric%stdevFlag = .false. 
    endif

  end subroutine LVT_initensPercentile

!BOP
! 
! !ROUTINE: LVT_diagnoseensPercentile
! \label{LVT_diagnoseensPercentile}
!
! !INTERFACE: 
  subroutine LVT_diagnoseensPercentile(pass)
! 
! !USES:     

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the computations for 
!   calculating the ensPercentile of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelensPercentile](\ref{diagnoseSingleModelensPercentile})
!     updates the ensPercentile computation for a single variable. This routine
!     stores the precip values to be used later for computing ensPercentile, 
!     during the first pass through the data. 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
    integer       :: pass
!EOP
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(LVT_rc%startmode.eq."coldstart") then 
       if(pass.eq.1) then 
          if(LVT_metrics%ensPercentile%selectOpt.eq.1.or.&
               LVT_metrics%ensPercentile%timeOpt.eq.1) then 
             
             call LVT_getDataStream1Ptr(model)
             call LVT_getDataStream2Ptr(obs)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(model))   
                call diagnoseSingleensPercentileParams(model, obs, stats,&
                     LVT_metrics%ensPercentile)
                
                model => model%next
                obs => obs%next
                stats => stats%next
                
             enddo
          endif
       elseif(pass.eq.2) then 
          if(LVT_metrics%ensPercentile%selectOpt.eq.1.or.&
               LVT_metrics%ensPercentile%timeOpt.eq.1) then 
             
             call LVT_getDataStream1Ptr(model)
             call LVT_getDataStream2Ptr(obs)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(model))   
                
                call diagnoseSingleensPercentile(model, obs, stats,&
                     LVT_metrics%ensPercentile)
                
                model => model%next
                obs => obs%next
                stats => stats%next
                
             enddo
          endif
          
       endif
    elseif(LVT_rc%startmode.eq."restart") then 
       if(LVT_metrics%ensPercentile%selectOpt.eq.1.or.&
            LVT_metrics%ensPercentile%timeOpt.eq.1) then 
          
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)
          
          do while(associated(model))   
             
             call diagnoseSingleensPercentile(model, obs, stats,&
                  LVT_metrics%ensPercentile)
             
             model => model%next
             obs => obs%next
             stats => stats%next
             
          enddo
       endif
    endif
  end subroutine LVT_diagnoseensPercentile


!BOP
! 
! !ROUTINE: diagnoseSingleensPercentileParams
! \label{diagnoseSingleensPercentileParams}
!
! !INTERFACE: 
  subroutine diagnoseSingleensPercentileParams(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the ensPercentile computation of the 
!   specified variable. 
!
!  The arguments are: 
!
!  \begin{description}
!   \item[model] model variable object
!   \item[stats] object to hold the updated statistics
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: model
    type(LVT_metadataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!EOP
    integer                 :: t,k,tind
    type(ESMF_Time)         :: startTime, currTime
    type(ESMF_TimeInterval) :: ts
    character*3             :: ftime
    character*100           :: filename
    integer                 :: tindex, day1,day2, day3, day4, day5
    integer                 :: yr_index
    integer                 :: p1,p2,p3,p4,p5
    logical                 :: skip_flag, leap_year
    integer                 :: ftn
    integer                 :: rc    
    real                    :: min_value

    real                    :: pctile_model_day1(LVT_LIS_rc(1)%ntiles,&
         LVT_ensPctile_struc%nsize_total)
    real                    :: pctile_model_day2(LVT_LIS_rc(1)%ntiles,&
         LVT_ensPctile_struc%nsize_total)
    real                    :: pctile_model_day3(LVT_LIS_rc(1)%ntiles,&
         LVT_ensPctile_struc%nsize_total)
    real                    :: pctile_model_day4(LVT_LIS_rc(1)%ntiles,&
         LVT_ensPctile_struc%nsize_total)
    real                    :: pctile_model_day5(LVT_LIS_rc(1)%ntiles,&
         LVT_ensPctile_struc%nsize_total)


    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then        

       call LVT_getYearIndex(LVT_rc, yr_index, leap_year)

       skip_flag = .false.
       if(.not.leap_year) then 
          tindex = LVT_rc%doy
       else
          if(LVT_rc%doy.gt.60) then 
             tindex = LVT_rc%doy - 1
          elseif(LVT_rc%doy.eq.60) then  !29th
             skip_flag = .true. 
          else
             tindex = LVT_rc%doy
          endif
       endif
       
       if(.not.skip_flag) then 
          if(tindex-2.lt.1) then 
             day1 = 365+tindex-2
          else                      
             day1 = tindex - 2
          endif
          if(tindex-1.lt.1) then 
             day2 = 365+tindex-1
          else                      
             day2 = tindex - 1
          endif
          
          day3 = tindex
          
          if(tindex+1.gt.365) then 
             day4 = tindex+1-365
          else
             day4 = tindex+1
          endif
          
          if(tindex+2.gt.365) then 
             day5 = tindex+2-365
          else
             day5 = tindex+2
          endif
                      
          p1  = (yr_index -1)*LVT_ensPctile_struc%nwindow +1
          p2  = (yr_index -1)*LVT_ensPctile_struc%nwindow +2
          p3  = (yr_index -1)*LVT_ensPctile_struc%nwindow +3
          p4  = (yr_index -1)*LVT_ensPctile_struc%nwindow +4
          p5  = (yr_index -1)*LVT_ensPctile_struc%nwindow +5
!open 5 files corresponding to the 5 days. 

          write(unit=ftime,fmt='(I3.3)') day1
          filename = trim(LVT_rc%statsodir)//'/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          read(ftn) pctile_model_day1
          call LVT_releaseUnitNumber(ftn)

          write(unit=ftime,fmt='(I3.3)') day2
          filename = trim(LVT_rc%statsodir)//'/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          read(ftn) pctile_model_day2
          call LVT_releaseUnitNumber(ftn)

          write(unit=ftime,fmt='(I3.3)') day3
          filename = trim(LVT_rc%statsodir)//'/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          read(ftn) pctile_model_day3
          call LVT_releaseUnitNumber(ftn)

          write(unit=ftime,fmt='(I3.3)') day4
          filename = trim(LVT_rc%statsodir)//'/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          read(ftn) pctile_model_day4
          call LVT_releaseUnitNumber(ftn)
         
          write(unit=ftime,fmt='(I3.3)') day5
          filename = trim(LVT_rc%statsodir)//'/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          read(ftn) pctile_model_day5
          call LVT_releaseUnitNumber(ftn)

          if(trim(model%short_name).eq."SWE") then 
             min_value = 1 !1mm
          elseif(trim(model%short_name).eq."Qle".or.&
               trim(model%short_name).eq."Evap") then 
             min_value = -50
          else
             min_value = 0.0 !soil moisture
          endif

          do t=1,LVT_LIS_rc(1)%ntiles
             if(model%count(t,1).gt.0) then
                if(metric%selectOpt.eq.1) then 
                   if(model%value(t,1).ne.LVT_rc%udef.and.&
                        model%value(t,1).gt.min_value) then                    
!ignore if the data goes outside the allotted time window                      
                      if(p1.le.LVT_ensPctile_struc%nsize_total) then 
                         pctile_model_day1(t,p1)  = & 
                              model%value(t,1)
                      endif
                      
                      
                      if(p2.le.LVT_ensPctile_struc%nsize_total) then 
                         pctile_model_day2(t,p2)   = & 
                              model%value(t,1)
                      endif

                      if(p3.le.LVT_ensPctile_struc%nsize_total) then 
                         pctile_model_day3(t,p3)   = & 
                              model%value(t,1)
                      endif


                      if(p4.le.LVT_ensPctile_struc%nsize_total) then 
                         pctile_model_day4(t,p4) = &
                              model%value(t,1)
                      endif

                      if(p5.le.LVT_ensPctile_struc%nsize_total) then 
                         pctile_model_day5(t,p5) = &
                              model%value(t,1)
                      endif
                   endif
                endif
             endif
          enddo
          
          write(unit=ftime,fmt='(I3.3)') day1
          filename = trim(LVT_rc%statsodir)//'/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          write(ftn) pctile_model_day1
          call LVT_releaseUnitNumber(ftn)

          write(unit=ftime,fmt='(I3.3)') day2
          filename = trim(LVT_rc%statsodir)//'/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          write(ftn) pctile_model_day2
          call LVT_releaseUnitNumber(ftn)

          write(unit=ftime,fmt='(I3.3)') day3
          filename = trim(LVT_rc%statsodir)//'/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          write(ftn) pctile_model_day3
          call LVT_releaseUnitNumber(ftn)

          write(unit=ftime,fmt='(I3.3)') day4
          filename = trim(LVT_rc%statsodir)//'/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          write(ftn) pctile_model_day4
          call LVT_releaseUnitNumber(ftn)
         
          write(unit=ftime,fmt='(I3.3)') day5
          filename = trim(LVT_rc%statsodir)//'/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          write(ftn) pctile_model_day5
          call LVT_releaseUnitNumber(ftn)
       endif
       
    endif
  end subroutine diagnoseSingleensPercentileParams

!BOP
! 
! !ROUTINE: diagnoseSingleensPercentile
! \label{diagnoseSingleensPercentile}
!
! !INTERFACE: 
  subroutine diagnoseSingleensPercentile(model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine diagnoses the ensPercentile values for each grid cell. 
!  The arguments are: 
!
!  \begin{description}
!    \item[model] model variable object
!    \item[stats] object to hold the updated statistics
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    logical                 :: alarm
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!EOP

    integer  :: t,l,k
    real                    :: min_value
    character*3             :: ftime
    character*100           :: filename
    integer  :: yr_index, tindex
    logical  :: leap_year, skip_flag
    real     :: pvalue
    integer                 :: ftn
    real     :: pctile_model(LVT_LIS_rc(1)%ntiles,&
         LVT_ensPctile_struc%nsize_total)

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then 

       call LVT_getYearIndex(LVT_rc, yr_index, leap_year)

       skip_flag = .false.
       if(.not.leap_year) then 
          tindex = LVT_rc%doy
       else
          if(LVT_rc%doy.gt.60) then 
             tindex = LVT_rc%doy - 1
          elseif(LVT_rc%doy.eq.60) then  !29th
             skip_flag = .true. 
          else
             tindex = LVT_rc%doy
          endif
       endif

       if(.not.skip_flag) then 
!read the day's file to compute ensPercentile values. 
             write(unit=ftime,fmt='(I3.3)') tindex
             filename = trim(LVT_rc%statsodir)//'/'//&
                  trim(ftime)//'_pctile_climo.bin'
             ftn = LVT_getNextUnitNumber()
             open(ftn,file=trim(filename),form='unformatted')
             read(ftn) pctile_model
             call LVT_releaseUnitNumber(ftn)             

             if(trim(model%short_name).eq."SWE") then 
                min_value = 1 !1mm
             elseif(trim(model%short_name).eq."Qle".or.&
                  trim(model%short_name).eq."Evap") then 
                min_value = -50
             else
                min_value = 0.0 !soil moisture
             endif
             
             do t=1,LVT_LIS_rc(1)%ntiles
                if(model%value(t,1).gt.min_value) then 
                   call ensPercentile_value(&
                        model%value(t,1),&
                        stats%enspctile%value_model_nsize(t,tindex),&                  
                        pctile_model(t,:),&
                        pvalue)
                   
                   stats%enspctile%value_model(t,1,1) = & 
                        stats%enspctile%value_model(t,1,1) + pvalue
                   stats%enspctile%value_count_model(t,1,1) = &
                        stats%enspctile%value_count_model(t,1,1) +1
                endif
             enddo
          endif
       endif
     end subroutine diagnoseSingleensPercentile

!BOP
! 
! !ROUTINE: LVT_computeensPercentile
! \label{LVT_computeensPercentile}
!
! !INTERFACE: 
  subroutine LVT_computeensPercentile(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the ensPercentile values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelensPercentile](\ref{computeSingleModelensPercentile})
!     computes the ensPercentile values for a single variable
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
    integer               :: pass
    logical               :: alarm
!EOP
    integer         :: i 
    integer         :: index
    type(ESMF_Time) :: currTime
    integer         :: rc
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(LVT_rc%startmode.eq."coldstart") then 
       if(pass.eq.1) then 
          if(LVT_metrics%ensPercentile%selectOpt.eq.1.or.&
               LVT_metrics%ensPercentile%timeOpt.eq.1) then 
             call LVT_getDataStream1Ptr(model)
             call LVT_getDataStream2Ptr(obs)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(model))             
                call computeSingleensPercentileparams(alarm,&
                     model,obs,stats,LVT_metrics%ensPercentile)
                
                model => model%next
                obs => obs%next
                stats => stats%next
                
             enddo
          endif
       elseif(pass.eq.2) then        
          
          if(LVT_metrics%ensPercentile%selectOpt.eq.1.or.&
               LVT_metrics%ensPercentile%timeOpt.eq.1) then 
             if(alarm) then 
                if(LVT_metrics%ensPercentile%timeOpt.eq.1.and.&
                     LVT_metrics%ensPercentile%extractTS.eq.1) then 
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%ensPercentile%ftn_ts_loc(i),&
                           200,advance='no') &
                           LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                           LVT_rc%hr,'',LVT_rc%mn, '' 
                      write(LVT_ensPctile_struc%ftn_ts_darea(i),&
                           200,advance='no') &
                           LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                           LVT_rc%hr,'',LVT_rc%mn, '' 
                   enddo
                endif
             endif
200          format(I4, a1, I2.2, a1, I2.2, a1, I2.2, a1, I2.2,a1)
          
             call LVT_getDataStream1Ptr(model)
             call LVT_getDataStream2Ptr(obs)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(model)) 
                call computeSingleensPercentile(alarm,&
                     model,obs,stats,LVT_metrics%ensPercentile)
                
                model => model%next
                obs => obs%next
                stats => stats%next
                
             enddo
             
             if(alarm) then 
                if(LVT_metrics%ensPercentile%timeOpt.eq.1.and.&
                     LVT_metrics%ensPercentile%extractTS.eq.1) then 
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%ensPercentile%ftn_ts_loc(i),fmt='(a1)') ''
                      write(LVT_ensPctile_struc%ftn_ts_darea(i),fmt='(a1)') ''
                   enddo
                endif
             endif
          endif
       endif
    elseif(LVT_rc%startmode.eq."restart") then 

       if(LVT_metrics%ensPercentile%selectOpt.eq.1.or.&
            LVT_metrics%ensPercentile%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%ensPercentile%timeOpt.eq.1.and.&
                  LVT_metrics%ensPercentile%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ensPercentile%ftn_ts_loc(i),&
                        200,advance='no') &
                        LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                        LVT_rc%hr,'',LVT_rc%mn, '' 
                   write(LVT_ensPctile_struc%ftn_ts_darea(i),&
                        200,advance='no') &
                        LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                        LVT_rc%hr,'',LVT_rc%mn, '' 
                enddo
             endif
          endif
          
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)
          
          do while(associated(model)) 
             call computeSingleensPercentile(alarm,&
                  model,obs,stats,LVT_metrics%ensPercentile)
             
             model => model%next
             obs => obs%next
             stats => stats%next
             
          enddo
          
          if(alarm) then 
             if(LVT_metrics%ensPercentile%timeOpt.eq.1.and.&
                  LVT_metrics%ensPercentile%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ensPercentile%ftn_ts_loc(i),fmt='(a1)') ''
                   write(LVT_ensPctile_struc%ftn_ts_darea(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_computeensPercentile
  

!BOP
! 
! !ROUTINE: computeSingleensPercentileparams
! \label{computeSingleensPercentileparams}
!
! !INTERFACE: 
  subroutine computeSingleensPercentileparams(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the ensPercentile values for each grid cell. 
!  The arguments are: 
!
!  \begin{description}
!    \item[model] model variable object
!    \item[stats] object to hold the updated statistics
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    logical                 :: alarm
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!EOP

    character*100           :: filename
    character*3             :: ftime
    integer                 :: ftn
    real     :: pctile_model(LVT_LIS_rc(1)%ntiles,&
         LVT_ensPctile_struc%nsize_total)
    integer  :: t,l,k

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do k=1,365
!read each day's file, sort it and write it back. 
             write(unit=ftime,fmt='(I3.3)') k
             filename = trim(LVT_rc%statsodir)//'/'//&
                  trim(ftime)//'_pctile_climo.bin'
             ftn = LVT_getNextUnitNumber()
             open(ftn,file=trim(filename),form='unformatted')
             read(ftn) pctile_model
             call LVT_releaseUnitNumber(ftn)             

             do t=1,LVT_LIS_rc(1)%ntiles
                call compress_and_sort(&
                     LVT_ensPctile_struc%nsize_total,&
                     pctile_model(t,:),&
                     stats%enspctile%value_model_nsize(t,k))
             enddo

             write(unit=ftime,fmt='(I3.3)') k
             filename = trim(LVT_rc%statsodir)//'/'//&
                  trim(ftime)//'_pctile_climo.bin'
             ftn = LVT_getNextUnitNumber()
             open(ftn,file=trim(filename),form='unformatted')
             write(ftn) pctile_model

             call LVT_releaseUnitNumber(ftn)                          
          enddo

       endif
    endif

  end subroutine computeSingleensPercentileparams

!BOP
!
! !ROUTINE: compress_and_sort
! \label{compress_and_sort}
!
! !INTERFACE: 
  subroutine compress_and_sort(size_in, value_in, size_out)
! !USES: 
    use LVT_SortingMod, only : LVT_sort
! !ARGUMENTS: 
    integer             :: size_in
    real                :: value_in(size_in)
    integer             :: size_out
! 
! !DESCRIPTION: 
!   This subroutine compresses the input array after removing
!   undefined values and sorts it in the ascending order. 
! 
!EOP

    real                :: value_out(size_in)
    integer             :: i 

!remove undefs, compute size_out

    value_out = LVT_rc%udef

    size_out = 0 
    do i=1,size_in
       if(value_in(i).ne.LVT_rc%udef) then 
          size_out = size_out + 1
          value_out(size_out) = value_in(i)
       endif
    enddo
!sort value_out and set it back to value_in
    call LVT_sort(value_out, size_out)
    
    value_in = value_out
    
  end subroutine compress_and_sort

!BOP
!
! !ROUTINE: compute_sorted_sizes
! \label{compute_sorted_sizes}
!
! !INTERFACE: 
  subroutine compute_sorted_sizes(size_in, value_in, size_out)
! !USES: 
    use LVT_SortingMod, only : LVT_sort
! !ARGUMENTS: 
    integer             :: size_in
    real                :: value_in(size_in)
    integer             :: size_out
! 
! !DESCRIPTION: 
!   This subroutine compresses the input array after removing
!   undefined values and sorts it in the ascending order. 
! 
!EOP

    integer             :: i 

!remove undefs, compute size_out

    size_out = 0 
    do i=1,size_in
       if(value_in(i).ne.LVT_rc%udef) then 
          size_out = size_out + 1
       endif
    enddo
  end subroutine compute_sorted_sizes

!BOP
! 
! !ROUTINE: computeSingleensPercentile
! \label{computeSingleensPercentile}
!
! !INTERFACE: 
  subroutine computeSingleensPercentile(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the ensPercentile values for each grid cell. 
!  The arguments are: 
!
!  \begin{description}
!    \item[model] model variable object
!    \item[stats] object to hold the updated statistics
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    logical                 :: alarm
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!EOP

    real     :: pctile_model(LVT_LIS_rc(1)%ntiles,&
         LVT_ensPctile_struc%nsize_total)
    integer  :: t,l,k,i,kk,m,tid,stid
    integer  :: sumv(5, LVT_rc%ntslocs)
    real     :: ensPercentile_darea(5, LVT_rc%ntslocs)
    integer  :: count_ensPercentile_model(LVT_LIS_rc(1)%ntiles,model%selectNlevs,1)
    integer  :: count_ensPercentile_obs(LVT_LIS_rc(1)%ntiles,obs%selectNlevs,1)

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then 

       do t=1,LVT_LIS_rc(1)%ntiles
          do k=1,model%selectNlevs
             if(stats%enspctile%value_count_model(t,k,1).gt.0) then 
                stats%enspctile%value_model(t,k,1) = &
                     stats%enspctile%value_model(t,k,1)/&
                     stats%enspctile%value_count_model(t,k,1)
             else
                stats%enspctile%value_model(t,k,1) = LVT_rc%udef
             endif
          enddo
       enddo

       count_ensPercentile_model = 1

       if(metric%extractTS.eq.1) then 
          call LVT_writeTSinfo(metric%ftn_ts_loc,&
               model,&
               LVT_LIS_rc(1)%ntiles,&
               stats%enspctile%value_model,&
               count_ensPercentile_model)
       endif

       do k=1,model%selectNlevs
          call LVT_computeCI(stats%enspctile%value_model(:,k,1),&
               LVT_LIS_rc(1)%ntiles,&
               LVT_rc%pval_CI, stats%enspctile%value_model_ci(k,1))
       enddo

          
       do k=1,model%selectNlevs
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_LIS_rc(1)%ntiles, stats%enspctile%value_model)
       enddo

       do k=1,model%selectNlevs
          do i=1,LVT_rc%ntslocs
             sumv(:,i) = 0
             do kk=1,LVT_TSobj(i)%npts
                tid = LVT_TSobj(i)%ts_tindex(kk)
                if(tid.ne.-1) then 
                   stid = (tid-1)*LVT_LIS_rc(1)%nensem
                   do m=1,LVT_LIS_rc(1)%nensem
                      t = stid+m
                      if(stats%enspctile%value_model(t,k,1).ne.LVT_rc%udef) then 
                         if(stats%enspctile%value_model(t,k,1).lt.&
                              LVT_ensPctile_struc%d_classes(1)) then 
                            sumv(1,i) = sumv(1,i) + 1
                         endif
                         if(stats%enspctile%value_model(t,k,1).lt.&
                              LVT_ensPctile_struc%d_classes(2)) then 
                            sumv(2,i) = sumv(2,i) + 1
                         endif
                         if(stats%enspctile%value_model(t,k,1).lt.&
                              LVT_ensPctile_struc%d_classes(3)) then 
                            sumv(3,i) = sumv(3,i) + 1
                         endif
                         if(stats%enspctile%value_model(t,k,1).lt.&
                              LVT_ensPctile_struc%d_classes(4)) then 
                            sumv(4,i) = sumv(4,i) + 1
                         endif
                         if(stats%enspctile%value_model(t,k,1).lt.&
                              LVT_ensPctile_struc%d_classes(5)) then 
                            sumv(5,i) = sumv(5,i) + 1
                         endif
                      endif
                   enddo
                endif
             enddo
             if(LVT_TSobj(i)%npts.gt.0) then
                ensPercentile_darea(1,i) = real(sumv(1,i))*100.0/&
                     real(LVT_TSobj(i)%npts)
                ensPercentile_darea(2,i) = real(sumv(2,i))*100.0/&
                     real(LVT_TSobj(i)%npts)
                ensPercentile_darea(3,i) = real(sumv(3,i))*100.0/&
                     real(LVT_TSobj(i)%npts)
                ensPercentile_darea(4,i) = real(sumv(4,i))*100.0/&
                     real(LVT_TSobj(i)%npts)
                ensPercentile_darea(5,i) = real(sumv(5,i))*100.0/&
                     real(LVT_TSobj(i)%npts)
             else
                ensPercentile_darea(:,i) = LVT_rc%udef
             endif
          enddo
       enddo
       
       do i=1,LVT_rc%ntslocs
          write(LVT_ensPctile_struc%ftn_ts_darea(i),203,advance='no') &
               (ensPercentile_darea(k,i),k=1,5)
       enddo

203 format(5E14.6)


    endif
  end subroutine computeSingleensPercentile
  
  subroutine ensPercentile_value(& 
       value_in, &
       nsize, &
       value_sorted, &
       ensPercentile)

    real    :: value_in
    integer :: nsize
    real    :: value_sorted(nsize)
    real    :: ensPercentile

    integer :: k 
    integer :: nx

    nx = 1
    k = 1
    do while(value_sorted(k).lt.value_in & 
         .and.k.le.nsize)  
       k = k + 1
    enddo

    nx = k
    ensPercentile = real(nx)/real(nsize+1)
  end subroutine ensPercentile_value
!BOP
! 
! !ROUTINE: LVT_writeMetric_ensPercentile
! \label(LVT_writeMetric_ensPercentile)
!
! !INTERFACE:
  subroutine LVT_writeMetric_ensPercentile(pass,final,vlevels,stats,obs)
! 
! !USES:   
    use LVT_statsMod, only : LVT_writeSummaryStats
    use LVT_pluginIndices
    implicit none
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

    integer                 :: pass
    integer                 :: final
    integer                 :: vlevels
    type(LVT_statsEntry)    :: stats
    type(LVT_metaDataEntry) :: obs

    integer                 :: k,l,tind
!    real                    :: ensPercentile_model(LVT_LIS_rc(1)%ntiles,1,1)
!    real                    :: ensPercentile_obs(LVT_LIS_rc(1)%ntiles,1,1)
!    integer                 :: count_ensPercentile_model(LVT_LIS_rc(1)%ntiles,1,1)
!    integer                 :: count_ensPercentile_obs(LVT_LIS_rc(1)%ntiles,1,1)

!    count_ensPercentile_model = 1
!    count_ensPercentile_obs = 1

    if(pass.eq.LVT_metrics%ensPercentile%npass) then 
       if(final.ne.1) then
           if(stats%selectOpt.eq.1) then 
              do k=1,vlevels
                if(LVT_metrics%ensPercentile%timeOpt.eq.1) then 
                   call LVT_writevar_gridded(LVT_metrics%ensPercentile%ftn_ts, &
                        stats%enspctile%value_model(:,k,1),&
                        stats%vid_ts(LVT_ensPercentileid,1),1.0,dim1=k)
                endif
             enddo

          endif
       endif
    endif

  end subroutine LVT_writeMetric_ensPercentile

!BOP
! 
! !ROUTINE: LVT_resetMetric_ensPercentile
! \label(LVT_resetMetric_ensPercentile)
!
! !INTERFACE:
  subroutine LVT_resetMetric_ensPercentile
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine resets the arrays that stores the precipitation
!   values. The reset is done at every stats output writing interval
!   to get the arrays reinitialized for the next set of time series
!   computations.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer                :: i,k,l
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)

    do while(associated(model))
       if(stats%selectOpt.eq.1) then 
          if(LVT_metrics%ensPercentile%timeOpt.eq.1) then 
             stats%enspctile%value_model  = 0
             stats%enspctile%value_count_model = 0 
          endif
       endif

       model => model%next
       obs   => obs%next
       stats => stats%next

    enddo

  end subroutine LVT_resetMetric_ensPercentile

!BOP
! 
! !ROUTINE: LVT_writerestart_ensPercentile
! 
! !INTERFACE:
  subroutine LVT_writerestart_ensPercentile(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for ensPercentile metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    integer                          :: k,l
    type(LVT_metaDataEntry), pointer :: model
    type(LVT_metaDataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats

#if 0 
    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)

    do while(associated(model))
       if(LVT_metrics%ensPercentile%selectOpt.eq.1) then 
          do l=1,LVT_ensPctile_struc%nsize_total
             call LVT_writevar_restart(ftn,&
                  LVT_ensPctile_struc%model_final(:,l),tileflag=1)
          enddo
          do l=1,LVT_ensPctile_struc%nasc
             call LVT_writevar_restart(ftn,&
                  LVT_ensPctile_struc%model_mu(:,l),tileflag=1)
             call LVT_writevar_restart(ftn,&
                  LVT_ensPctile_struc%model_sigma(:,l),tileflag=1)
          enddo
       end if

       model => model%next
       obs   => obs%next
       stats => stats%next
    enddo
#endif
  end subroutine LVT_writerestart_ensPercentile

!BOP
! 
! !ROUTINE: LVT_readrestart_ensPercentile
! 
! !INTERFACE:
  subroutine LVT_readrestart_ensPercentile(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for ensPercentile metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    type(LVT_metaDataEntry), pointer :: model
    type(LVT_metaDataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats
    integer              :: k,l
#if 0 

    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(model%short_name.eq."RootMoist") then 
          if(LVT_metrics%ensPercentile%selectOpt.eq.1) then 
             do l=1,LVT_ensPctile_struc%nsize_total
                call LVT_readvar_restart(ftn,&
                     LVT_ensPctile_struc%model_final(:,l),tileflag=1)
             enddo
             do l=1,LVT_ensPctile_struc%nasc
                call LVT_readvar_restart(ftn,&
                     LVT_ensPctile_struc%model_mu(:,l),tileflag=1)
                call LVT_readvar_restart(ftn,&
                     LVT_ensPctile_struc%model_sigma(:,l),tileflag=1)
             enddo
          end if
       end if
       
       model => model%next
       obs   => obs%next
       stats => stats%next
    enddo
#endif
  end subroutine LVT_readrestart_ensPercentile

end module LVT_ensPercentileMod
