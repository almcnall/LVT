!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_FARMod
! \label(LVT_FARMod)
!
! !INTERFACE:
module LVT_FARMod
! 
! !USES:   
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
  use LVT_TSMod
  use LVT_logMod
  use LVT_CIMod

  implicit none

  private
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
!  !DESCRIPTION: 
!   This module handles the False Alarm Ratio (FAR) 
!   computations by comparing 
!   the LIS output to the specified observations. 
!
!     FAR = Missed negatives / (Missed negatives + Hits)
! 
! !FILES USED:
!
!  !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initFAR
  public :: LVT_diagnoseFAR
  public :: LVT_computeFAR 
  public :: LVT_writeMetric_FAR
  public :: LVT_resetMetric_FAR
  public :: LVT_writerestart_FAR
  public :: LVT_readrestart_FAR

contains
  
!BOP
! 
! !ROUTINE: LVT_initFAR
! \label{LVT_initFAR}
!
! !INTERFACE: 
  subroutine LVT_initFAR(model,obs,stats,metric)
! 
! !USES:   
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
!BOP
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!EOP
    if(metric%selectOpt.eq.1.or.metric%timeopt.eq.1) then 
       allocate(stats%far%value_total(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))
       allocate(stats%far%count_value_total(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))   
       stats%far%value_total = 0.0
       stats%far%count_value_total = 0

       allocate(stats%far%value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%far%value_ci = LVT_rc%udef    

       allocate(stats%far%value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       allocate(stats%far%count_value_ts(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       

       stats%far%value_ts = 0.0
       stats%far%count_value_ts = 0 

       allocate(stats%far%tavg_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       allocate(stats%far%tavg_count_value_ts(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       

       stats%far%tavg_value_ts = 0.0
       stats%far%tavg_count_value_ts = 0 

       if(metric%computeSC.eq.1) then 
          allocate(stats%far%value_asc(LVT_rc%ngrid, model%selectNlevs,LVT_rc%nasc))
          allocate(stats%far%count_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%far%value_asc = 0.0
          stats%far%count_value_asc = 0
       endif
       if(metric%computeADC.eq.1) then 
          allocate(stats%far%value_adc(LVT_rc%ngrid, model%selectNlevs,LVT_rc%nadc))
          allocate(stats%far%count_value_adc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nadc))
          stats%far%value_adc = 0.0
          stats%far%count_value_adc = 0
       endif
    endif
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1
    metric%obsData = .false. 
    metric%stdevFlag = .false. 

  end subroutine LVT_initFAR
  
!BOP
! 
! !ROUTINE: LVT_diagnoseFAR
! \label{LVT_diagnoseFAR}
!
! !INTERFACE: 
  subroutine LVT_diagnoseFAR(pass)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the FAR calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleFAR](\ref{diagnoseSingleFAR})
!     updates the FAR computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none
    integer                 :: pass
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%far%selectOpt.eq.1.or.&
            LVT_metrics%far%timeOpt.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleFAR(obs, model, stats, &
                  LVT_metrics%far)

             model => model%next
             obs => obs%next
             stats => stats%next
             
          enddo
       endif
    endif

  end subroutine LVT_diagnoseFAR

!BOP
! 
! !ROUTINE: diagnoseSingleFAR
! \label{diagnoseSingleFAR}
!
! !INTERFACE: 
  subroutine diagnoseSingleFAR(obs, model, stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the FAR computation (updates the running 
!  sum calculations of the squared error) 
!  The arguments are: 
!
!  \begin{description}
!   \item[obs] observation object
!   \item[model] model variable object
!   \item[stats] object to hold the updated statistics
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none

    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry) :: stats
    type(LVT_metricEntry)   :: metric

    integer    :: t,k
    real       :: aval,bval

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             if(trim(obs%units).eq.trim(model%units)) then 
                if(obs%count(t,k).ne.0.and. &
                     model%count(t,k).ne.0) then                    
                   aval = 0.0
                   bval = 0.0                  
                   if(metric%selectOpt.eq.1.or.metric%timeOpt.eq.1) then
                      if(model%value(t,k).gt.metric%threshold.and.&
                           obs%value(t,k).gt.metric%threshold) then 
                         aval = aval + 1
                      endif
                      if(model%value(t,k).gt.metric%threshold.and.&
                           obs%value(t,k).le.metric%threshold) then 
                         bval = bval + 1
                      endif
                      if(aval+bval.ne.0) then 
                         stats%far%value_ts(t,k,1)  = stats%far%value_ts(t,k,1) + &
                              bval/(aval+bval)
                         stats%far%count_value_ts(t,k,1) = &
                              stats%far%count_value_ts(t,k,1) + 1
                         stats%far%value_total(t,k,1)  = stats%far%value_total(t,k,1) + &
                              bval/(aval+bval)
                         stats%far%count_value_total(t,k,1) = &
                              stats%far%count_value_total(t,k,1) + 1
                      endif
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         aval = 0.0
                         bval = 0.0
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then
                            if(model%value(t,k).gt.metric%threshold.and.&
                                 obs%value(t,k).gt.metric%threshold) then 
                               aval = aval+1
                            endif
                            if(model%value(t,k).gt.metric%threshold.and.&
                                 obs%value(t,k).le.metric%threshold) then 
                               bval = bval+1
                            endif
                            if(aval+bval.ne.0) then 
                               stats%far%value_ts(t,k,2) = bval/(aval+bval)
                               stats%far%count_value_ts(t,k,2) = &
                                    stats%far%count_value_ts(t,k,2) + 1
                               stats%far%value_total(t,k,2)  = stats%far%value_total(t,k,2) + &
                                    bval/(aval+bval)
                               stats%far%count_value_total(t,k,2) = &
                                    stats%far%count_value_total(t,k,2) + 1
                            endif

                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then 
                            if(model%value(t,k).gt.metric%threshold.and.&
                                 obs%value(t,k).gt.metric%threshold) then 
                               aval = aval+1
                            endif
                            if(model%value(t,k).gt.metric%threshold.and.&
                                 obs%value(t,k).le.metric%threshold) then 
                               bval = bval+1
                            endif 
                            if(aval+bval.ne.0) then 
                               stats%far%value_ts(t,k,3) = bval/(aval+bval)
                               stats%far%count_value_ts(t,k,3) = &
                                    stats%far%count_value_ts(t,k,3) + 1
                               stats%far%value_total(t,k,3)  = stats%far%value_total(t,k,3) + &
                                    bval/(aval+bval)
                               stats%far%count_value_total(t,k,3) = &
                                    stats%far%count_value_total(t,k,3) + 1
                            endif
                         endif
                      endif
                   endif                   
                endif
             else
                write(LVT_logunit,*) 'For variable ',trim(model%standard_name)
                write(LVT_logunit,*) 'observations are in ',trim(obs%units)
                write(LVT_logunit,*) 'and LIS output is in ',trim(model%units)
                write(LVT_logunit,*) 'please add the support of ',&
                     trim(model%units), ' in the observation plugin'
                call LVT_endrun
             endif
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleFAR


!BOP
! 
! !ROUTINE: LVT_computeFAR
! \label{LVT_computeFAR}
!
! !INTERFACE: 
  subroutine LVT_computeFAR(pass,alarm)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute FAR values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleFAR](\ref{computeSingleFAR})
!     updates the FAR computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     FAR computation has been reached
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: i,pass

    logical     :: alarm
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
    if(LVT_metrics%far%selectOpt.eq.1.or.LVT_metrics%far%timeOpt.eq.1) then 
       if(alarm) then 
          if(LVT_metrics%far%timeOpt.eq.1.and.&
               LVT_metrics%far%extractTS.eq.1) then 
             do i=1,LVT_rc%ntslocs
                write(LVT_metrics%far%ftn_ts_loc(i),200,advance='no') &
                     LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                     LVT_rc%hr,'',LVT_rc%mn, '' 
             enddo
          endif
       endif
200    format(I4, a1, I2.2, a1, I2.2, a1, I2.2, a1, I2.2,a1)

       call LVT_getDataStream1Ptr(model)
       call LVT_getDataStream2Ptr(obs)
       call LVT_getstatsEntryPtr(stats)
       
       do while(associated(model))
          
          call computeSingleFAR(alarm,obs,model,stats, &
               LVT_metrics%far)

          model => model%next
          obs => obs%next
          stats => stats%next
       enddo
       
       if(alarm) then 
          if(LVT_metrics%far%timeOpt.eq.1.and.&
               LVT_metrics%far%extractTS.eq.1) then 
             do i=1,LVT_rc%ntslocs
                write(LVT_metrics%far%ftn_ts_loc(i),fmt='(a1)') ''
             enddo
          endif
       endif
    endif
 endif
  end subroutine LVT_ComputeFAR

!BOP
! 
! !ROUTINE: computeSingleFAR
! \label{computeSingleFAR}
!
! !INTERFACE: 
  subroutine computeSingleFAR(alarm,obs, model,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the FAR values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     FAR computation has been reached
!    \item[obs] observation object
!    \item[model] model variable object
!    \item[stats] object to hold the updated statistics
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none
    logical                 :: alarm
    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer  :: t,l,k,tind

    if((metric%selectOpt.eq.1.or.metric%timeOpt.eq.1)) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%far%count_value_ts(t,k,l).ne.0) then 
                      stats%far%value_ts(t,k,l) = stats%far%value_ts(t,k,l)/&
                           (stats%far%count_value_ts(t,k,l))

                      stats%far%tavg_value_ts(t,k,l) = &
                           stats%far%tavg_value_ts(t,k,l) + & 
                           stats%far%value_ts(t,k,l)
                      stats%far%tavg_count_value_ts(t,k,l) = & 
                           stats%far%tavg_count_value_ts(t,k,l) + 1
                   else
                      stats%far%value_ts(t,k,l) = LVT_rc%udef
                   endif
                enddo                   
                if(metric%computeSC.eq.1) then 
                   if(stats%far%value_ts(t,k,1).ne.LVT_rc%udef) then 
                      call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                           tind)
                      stats%far%value_asc(t,k,tind) = &
                           stats%far%value_asc(t,k,tind)+ &
                           stats%far%value_ts(t,k,1)
                      stats%far%count_value_asc(t,k,tind) = &
                           stats%far%count_value_asc(t,k,tind)+ 1
                   endif
                endif
                if(metric%computeADC.eq.1) then 
                   if(stats%far%value_ts(t,k,1).ne.LVT_rc%udef) then 
                      call LVT_getADCTimeIndex(tind)
                      stats%far%value_adc(t,k,tind) = &
                           stats%far%value_adc(t,k,tind)+ &
                           stats%far%value_ts(t,k,1)
                      stats%far%count_value_adc(t,k,tind) = &
                           stats%far%count_value_adc(t,k,tind)+ 1
                   endif                
                endif
             enddo
          enddo
          if(alarm) then 
             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%far%tavg_count_value_ts(t,k,l).ne.0) then 
                         stats%far%tavg_value_ts(t,k,l) = &
                              stats%far%tavg_value_ts(t,k,l)/&
                              (stats%far%tavg_count_value_ts(t,k,l))
                      else
                         stats%far%tavg_value_ts(t,k,l) = LVT_rc%udef
                      end if
                   enddo
                enddo
             enddo
             if(metric%extractTS.eq.1) then 
                call LVT_writeTSinfo(metric%ftn_ts_loc,&
                     model,&
                     LVT_rc%ngrid,&
                     stats%far%value_ts,&
                     stats%far%count_value_ts)
             endif
          endif
       endif
    endif


    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if((stats%far%count_value_total(t,k,l).ne.0.and.&
                        stats%far%count_value_total(t,k,l)&
                        .gt.LVT_rc%obsCountThreshold)) then 
                      
                      stats%far%value_total(t,k,l) = stats%far%value_total(t,k,l)/&
                           (stats%far%count_value_total(t,k,l))
                   else
                      stats%far%value_total(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1,LVT_rc%nasc
                      if(stats%far%count_value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then
                         stats%far%value_asc(t,k,l) = stats%far%value_asc(t,k,l)/&
                              stats%far%count_value_asc(t,k,l)
                      else
                         stats%far%value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
                if(metric%computeADC.eq.1) then 
                   do l=1,LVT_rc%nadc
                      if(stats%far%count_value_adc(t,k,l).gt.&
                           LVT_rc%ADCCountThreshold) then
                         stats%far%value_adc(t,k,l) = stats%far%value_adc(t,k,l)/&
                              stats%far%count_value_adc(t,k,l)
                      else
                         stats%far%value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
             enddo
          enddo
          
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%far%value_total(:,k,l),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%far%value_ci(k,l))
             enddo
          enddo
       
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid,stats%far%value_total)      
          if(metric%computeSC.eq.1) then 
             call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%far%value_asc,stats%far%count_value_asc)
          endif
          if(metric%computeADC.eq.1) then 
             call LVT_writeAvgDiurnalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%far%value_adc,stats%far%count_value_adc)
          endif
       endif
    endif

  end subroutine computeSingleFAR

!BOP
! 
! !ROUTINE: LVT_writeMetric_FAR
! \label{LVT_writeMetric_FAR}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_FAR(pass,final,vlevels,stats,obs)
! 
! !USES:
    use LVT_statsMod, only : LVT_writeSummaryStats
    use LVT_pluginIndices
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
!BOP
! !ARGUMENTS: 
    integer                 :: pass
    integer                 :: final
    integer                 :: vlevels
    type(LVT_statsEntry)   :: stats
    type(LVT_metaDataEntry) :: obs
!EOP
    integer                 :: l,tind
    integer                 :: k

    if(pass.eq.LVT_metrics%far%npass) then 
       if(final.ne.1) then 
          if(stats%selectOpt.eq.1) then 
             do k=1,vlevels
                if(LVT_metrics%far%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_gridded(LVT_metrics%far%ftn_ts, &
                           stats%far%value_ts(:,k,l),stats%vid_ts(LVT_FARid,1),1)
                      call LVT_writevar_gridded(LVT_metrics%far%ftn_ts, &
                           real(stats%far%count_value_ts(:,k,l)),&
                           stats%vid_count_ts(LVT_FARid,1),1)
                   enddo
                endif
             enddo
          endif
       else
          if(pass.eq.1) then 
             if(LVT_metrics%far%selectOpt.eq.1) then
                if(stats%selectOpt.eq.1) then 
                   do k=1,vlevels
                      do l=1,LVT_rc%strat_nlevels
                         
                         call LVT_writevar_gridded(LVT_metrics%far%ftn_total, &
                              stats%far%value_total(:,k,l),stats%vid_total(LVT_FARid,1))
                         call LVT_writevar_gridded(LVT_metrics%far%ftn_total, &
                              real(stats%far%count_value_total(:,k,l)),&
                              stats%vid_count_total(LVT_FARid,1))
                         
                      enddo
                      if(LVT_metrics%far%computeSC.eq.1) then 
                         do tind = 1,LVT_rc%nasc
                            call LVT_writevar_gridded(LVT_metrics%far%ftn_total,&
                                 stats%far%value_asc(:,k,tind),&
                                 stats%vid_sc_total(tind,LVT_FARid,1))
                         enddo
                      endif
                      
                      if(LVT_metrics%far%computeADC.eq.1) then 
                         do tind = 1,LVT_rc%nadc
                            call LVT_writevar_gridded(LVT_metrics%far%ftn_total,&
                                 stats%far%value_adc(:,k,tind),&
                                 stats%vid_adc_total(tind,LVT_FARid,1))
                         enddo
                      endif
                      
                      call LVT_writeSummaryStats(&
                           LVT_metrics%far%ftn_summ,&
                           LVT_metrics%far%short_name,&
                           LVT_rc%ngrid,&
                           stats%far%value_total(:,k,:), &
                           stats%far%count_value_total(:,k,:),stats%standard_name,&
                           stats%far%value_ci(k,:))
                   enddo
                endif
             endif
          endif
       end if
    endif
  end subroutine LVT_writeMetric_FAR


!BOP
! 
! !ROUTINE: LVT_resetMetric_FAR
! \label(LVT_resetMetric_FAR)
!
! !INTERFACE:!
  subroutine LVT_resetMetric_FAR(alarm)
 
! !INPUT PARAMETERS: 
    logical           :: alarm
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

    integer                 :: vlevels
    integer                 :: i,k,l
    type(LVT_metadataEntry), pointer :: model
    type(LVT_statsEntry)   , pointer :: stats


    call LVT_getDataStream1Ptr(model)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(stats%selectOpt.eq.1) then 
          do k=1,model%selectNlevs
             if(LVT_metrics%far%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   stats%far%value_ts(:,k,l) = 0.0
                   stats%far%count_value_ts(:,k,l)=0 
                   if(alarm) then 
                      stats%far%tavg_value_ts(:,k,l) = 0.0
                      stats%far%tavg_count_value_ts(:,k,l)=0 
                   endif
                enddo
             endif
          enddo
       endif

       model => model%next
       stats => stats%next

    enddo

  end subroutine LVT_resetMetric_FAR

!BOP
! 
! !ROUTINE: LVT_writerestart_FAR
! 
! !INTERFACE:
  subroutine LVT_writerestart_FAR(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for FAR metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%far%selectOpt.eq.1) then 
       
       print*, 'The writerestart method is not implemented for FAR'
       stop
    end if
    
  end subroutine LVT_writerestart_FAR

!BOP
! 
! !ROUTINE: LVT_readrestart_FAR
! 
! !INTERFACE:
  subroutine LVT_readrestart_FAR(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for FAR metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%far%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for FAR'
       stop
    end if
    
  end subroutine LVT_readrestart_FAR

end module LVT_FARMod

