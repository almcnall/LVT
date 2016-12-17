!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_ACCMod
! \label(LVT_ACCMod)
!
! !INTERFACE:
module LVT_ACCMod
! 
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
  use LVT_TSMod
  use LVT_logMod
  use LVT_CIMod

  implicit none
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the Accuracy measure (ACC)
!  computations by comparing 
!  the LIS output to the specified observations. 
! 
!  ACC = (Hits + Correct Negatives)/ Total
! 
! !FILES USED:
!
!  !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP
  private

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initACC
  public :: LVT_diagnoseACC
  public :: LVT_computeACC 
  public :: LVT_writeMetric_ACC
  public :: LVT_resetMetric_ACC
  public :: LVT_writerestart_ACC
  public :: LVT_readrestart_ACC

contains
  
!BOP
! 
! !ROUTINE: LVT_initACC
! \label{LVT_initACC}
!
! !INTERFACE:
  subroutine LVT_initACC(model,obs,stats,metric)
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
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    if(metric%selectOpt.eq.1) then 
       allocate(stats%acc%value_total(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))
       stats%acc%value_total = 0.0
       allocate(stats%acc%count_value_total(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       
       stats%acc%count_value_total = 0
       allocate(stats%acc%value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%acc%value_ci = LVT_rc%udef
    endif
    
    if(metric%timeopt.eq.1) then 
       allocate(stats%acc%value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%acc%value_ts = 0.0
       allocate(stats%acc%count_value_ts(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       
       stats%acc%count_value_ts = 0 

       allocate(stats%acc%tavg_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%acc%tavg_value_ts = 0.0
       allocate(stats%acc%tavg_count_value_ts(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       
       stats%acc%tavg_count_value_ts = 0 

       if(metric%computeSC.eq.1) then 
          allocate(stats%acc%value_asc(LVT_rc%ngrid, model%selectNlevs,LVT_rc%nasc))
          stats%acc%value_asc = 0.0
          allocate(stats%acc%count_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%acc%count_value_asc = 0
       endif
       if(metric%computeADC.eq.1) then 
          allocate(stats%acc%value_adc(LVT_rc%ngrid, model%selectNlevs,LVT_rc%nadc))
          stats%acc%value_adc = 0.0
          allocate(stats%acc%count_value_adc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nadc))
          stats%acc%count_value_adc = 0
       endif
    endif
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1
    metric%obsData = .false. 
    metric%stdevFlag = .false. 

  end subroutine LVT_initACC
  
!BOP
! 
! !ROUTINE: LVT_diagnoseACC
! \label{LVT_diagnoseACC}
!
! !INTERFACE: 
  subroutine LVT_diagnoseACC(pass)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the ACC calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleACC](\ref{diagnoseSingleACC})
!     updates the ACC computation for a single variable 
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
       if(LVT_metrics%acc%selectOpt.eq.1.or.&
            LVT_metrics%acc%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleACC(obs,model,stats,&
                  LVT_metrics%acc)

             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo
       endif
    endif

  end subroutine LVT_diagnoseACC

!BOP
! 
! !ROUTINE: diagnoseSingleACC
! \label{diagnoseSingleACC}
!
! !INTERFACE: 
  subroutine diagnoseSingleACC(obs, model, stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the ACC computation (updates the running 
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

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             if(obs%count(t,k).ne.0.and. &
                  model%count(t,k).ne.0) then      
                if(metric%selectOpt.eq.1) then
                   if((model%value(t,k).gt.metric%threshold.and.&
                        obs%value(t,k).gt.metric%threshold).or.&
                        (model%value(t,k).le.metric%threshold).and.&
                        (obs%value(t,k).le.metric%threshold))then 
                      stats%acc%value_total(t,k,1) = stats%acc%value_total(t,k,1)+1
                   endif
                   stats%acc%count_value_total(t,k,1) = &
                        stats%acc%count_value_total(t,k,1) + 1
                   
                   if(LVT_rc%strat_nlevels.gt.1) then 
                      if(LVT_stats%strat_var(t,k).gt.&
                           LVT_rc%strat_var_threshold) then 
                         if((model%value(t,k).gt.metric%threshold.and.&
                              obs%value(t,k).gt.metric%threshold).or.&
                              (model%value(t,k).le.metric%threshold).and.&
                              (obs%value(t,k).le.metric%threshold))then 
                            stats%acc%value_total(t,k,2) = stats%acc%value_total(t,k,2)+1
                         endif
                         stats%acc%count_value_total(t,k,2) = &
                              stats%acc%count_value_total(t,k,2) + 1
                         
                         
                      elseif(LVT_stats%strat_var(t,k).le.&
                           LVT_rc%strat_var_threshold) then 
                         
                         if((model%value(t,k).gt.metric%threshold.and.&
                              obs%value(t,k).gt.metric%threshold).or.&
                              (model%value(t,k).le.metric%threshold).and.&
                              (obs%value(t,k).le.metric%threshold))then 
                            stats%acc%value_total(t,k,3) = stats%acc%value_total(t,k,3)+1
                         endif
                         stats%acc%count_value_total(t,k,3) = &
                              stats%acc%count_value_total(t,k,3) + 1
                         
                      endif
                   endif
                endif
                if(metric%timeOpt.eq.1) then 
                   if((model%value(t,k).gt.metric%threshold.and.&
                        obs%value(t,k).gt.metric%threshold).or.&
                        (model%value(t,k).le.metric%threshold).and.&
                        (obs%value(t,k).le.metric%threshold))then 
                      stats%acc%value_ts(t,k,1) = stats%acc%value_ts(t,k,1)+1
                   endif
                   stats%acc%count_value_ts(t,k,1) = &
                        stats%acc%count_value_ts(t,k,1) + 1
                   if(LVT_rc%strat_nlevels.gt.1) then 
                      if(LVT_stats%strat_var(t,k).gt.&
                           LVT_rc%strat_var_threshold) then
                         if((model%value(t,k).gt.metric%threshold.and.&
                              obs%value(t,k).gt.metric%threshold).or.&
                              (model%value(t,k).le.metric%threshold).and.&
                              (obs%value(t,k).le.metric%threshold))then 
                            stats%acc%value_ts(t,k,2) = stats%acc%value_ts(t,k,2)+1
                         endif
                         stats%acc%count_value_ts(t,k,2) = &
                              stats%acc%count_value_ts(t,k,2) + 1
                      elseif(LVT_stats%strat_var(t,k).le.&
                           LVT_rc%strat_var_threshold) then 
                         if((model%value(t,k).gt.metric%threshold.and.&
                              obs%value(t,k).gt.metric%threshold).or.&
                              (model%value(t,k).le.metric%threshold).and.&
                              (obs%value(t,k).le.metric%threshold))then 
                            stats%acc%value_ts(t,k,3) = stats%acc%value_ts(t,k,3)+1
                         endif
                         stats%acc%count_value_ts(t,k,3) = &
                              stats%acc%count_value_ts(t,k,3) + 1
                      endif
                   endif
                endif
             endif
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleACC


!BOP
! 
! !ROUTINE: LVT_computeACC
! \label{LVT_computeACC}
!
! !INTERFACE: 
  subroutine LVT_computeACC(pass,alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute ACC values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleACC](\ref{computeSingleACC})
!     updates the ACC computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     ACC computation has been reached
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    implicit none

    integer               :: i, pass
    logical               :: alarm
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%acc%selectOpt.eq.1.or.LVT_metrics%acc%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%acc%timeOpt.eq.1.and.&
                  LVT_metrics%acc%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%acc%ftn_ts_loc(i),200,advance='no') &
                        LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                        LVT_rc%hr,'',LVT_rc%mn, '' 
                enddo
             endif
          endif
200       format(I4, a1, I2.2, a1, I2.2, a1, I2.2, a1, I2.2,a1)
          
          
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call computeSingleACC(alarm,obs,model,stats, &
                  LVT_metrics%acc)

             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%acc%timeOpt.eq.1.and.&
                  LVT_metrics%acc%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%acc%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_ComputeACC

!BOP
! 
! !ROUTINE: computeSingleACC
! \label{computeSingleACC}
!
! !INTERFACE: 
  subroutine computeSingleACC(alarm,obs, model,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the ACC values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     ACC computation has been reached
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

    if(metric%timeOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%acc%count_value_ts(t,k,l).ne.0) then 
                      stats%acc%value_ts(t,k,l) = stats%acc%value_ts(t,k,l)/&
                           stats%acc%count_value_ts(t,k,l)
                      stats%acc%tavg_value_ts(t,k,l) = &
                           stats%acc%tavg_value_ts(t,k,l)+ &
                           stats%acc%value_ts(t,k,l)
                      stats%acc%tavg_count_value_ts(t,k,l) = & 
                           stats%acc%tavg_count_value_ts(t,k,l) + 1
                   else
                      stats%acc%value_ts(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   if(stats%acc%count_value_ts(t,k,1).ne.0) then 
                      call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                           tind)
                      stats%acc%value_asc(t,k,tind) = &
                           stats%acc%value_asc(t,k,tind)+ &
                           stats%acc%value_ts(t,k,1)
                      stats%acc%count_value_asc(t,k,tind) = &
                           stats%acc%count_value_asc(t,k,tind)+ 1
                   endif
                endif
                if(metric%computeADC.eq.1) then 
                   if(stats%acc%count_value_ts(t,k,1).ne.0) then 
                      call LVT_getADCTimeIndex(tind)
                      stats%acc%value_adc(t,k,tind) = &
                           stats%acc%value_adc(t,k,tind)+ &
                           stats%acc%value_ts(t,k,1)
                      stats%acc%count_value_adc(t,k,tind) = &
                           stats%acc%count_value_adc(t,k,tind)+ 1
                   endif                
                endif                
             enddo
          enddo

          if(alarm) then 
             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%acc%tavg_count_value_ts(t,k,l).ne.0) then 
                         stats%acc%tavg_value_ts(t,k,l) = &
                              stats%acc%tavg_value_ts(t,k,l)/ & 
                              stats%acc%tavg_count_value_ts(t,k,l)
                      else
                         stats%acc%tavg_value_ts(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                enddo
             enddo
             if(metric%extractTS.eq.1) then 
                call LVT_writeTSinfo(metric%ftn_ts_loc,&
                     model,&
                     LVT_rc%ngrid,&
                     stats%acc%tavg_value_ts,&
                     stats%acc%tavg_count_value_ts)
             endif
          endif
       endif
    endif


    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%acc%count_value_total(t,k,l).ne.0.and.&
                        stats%acc%count_value_total(t,k,l)&
                        .gt.LVT_rc%obsCountThreshold) then 
                      
                      stats%acc%value_total(t,k,l) =&
                           stats%acc%value_total(t,k,l)/&
                           stats%acc%count_value_total(t,k,l)
                   else
                      stats%acc%value_total(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1,LVT_rc%nasc
                      if(stats%acc%count_value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then
                         stats%acc%value_asc(t,k,l) =&
                              stats%acc%value_asc(t,k,l)/&
                              stats%acc%count_value_asc(t,k,l)
                      else
                         stats%acc%value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
                if(metric%computeADC.eq.1) then 
                   do l=1,LVT_rc%nadc
                      if(stats%acc%count_value_adc(t,k,l).gt.&
                           LVT_rc%ADCCountThreshold) then
                         stats%acc%value_adc(t,k,l) = &
                              stats%acc%value_adc(t,k,l)/&
                              stats%acc%count_value_adc(t,k,l)
                      else
                         stats%acc%value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
             enddo
          enddo
          
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%acc%value_total(:,k,l),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%acc%value_ci(k,l))
             enddo
          enddo
       
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid,stats%acc%value_total)      
          if(metric%computeSC.eq.1) then 
             call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%acc%value_asc,stats%acc%count_value_asc)
          endif
          if(metric%computeADC.eq.1) then 
             call LVT_writeAvgDiurnalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%acc%value_adc,stats%acc%count_value_adc)
          endif
       endif
    endif

  end subroutine computeSingleACC

!BOP
! 
! !ROUTINE: LVT_writeMetric_ACC
! \label{LVT_writeMetric_ACC}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_ACC(pass,final,vlevels,stats,obs)
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
! !ARGUMENTS: 
    integer                 :: pass
    integer                 :: final
    integer                 :: vlevels
    type(LVT_statsEntry)    :: stats
    type(LVT_metaDataEntry) :: obs
!EOP
    integer                 :: l,tind
    integer                 :: k

    if(pass.eq.LVT_metrics%acc%npass) then 
       if(final.ne.1) then 
          if(stats%selectOpt.eq.1) then 
             do k=1,vlevels
                if(LVT_metrics%acc%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_gridded(LVT_metrics%acc%ftn_ts, &
                           stats%acc%value_ts(:,k,l),stats%vid_ts(LVT_ACCid,1),1)
                      call LVT_writevar_gridded(LVT_metrics%acc%ftn_ts, &
                           real(stats%acc%count_value_ts(:,k,l)),&
                           stats%vid_count_ts(LVT_ACCid,1),1)
                   enddo
                endif
             enddo
          endif
       else
          if(pass.eq.1) then 
             if(LVT_metrics%acc%selectOpt.eq.1) then
                if(stats%selectOpt.eq.1) then 
                   do k=1,vlevels
                      do l=1,LVT_rc%strat_nlevels
                         
                         call LVT_writevar_gridded(LVT_metrics%acc%ftn_total, &
                              stats%acc%value_total(:,k,l),stats%vid_total(LVT_ACCid,1))
                         call LVT_writevar_gridded(LVT_metrics%acc%ftn_total, &
                              real(stats%acc%count_value_total(:,k,l)),&
                              stats%vid_count_total(LVT_ACCid,1))
                         
                      enddo
                      if(LVT_metrics%acc%computeSC.eq.1) then 
                         do tind = 1,LVT_rc%nasc
                            call LVT_writevar_gridded(LVT_metrics%acc%ftn_total,&
                                 stats%acc%value_asc(:,k,tind),&
                                 stats%vid_sc_total(tind,LVT_ACCid,1))
                         enddo
                      endif
                      
                      if(LVT_metrics%acc%computeADC.eq.1) then 
                         do tind = 1,LVT_rc%nadc
                            call LVT_writevar_gridded(LVT_metrics%acc%ftn_total,&
                                 stats%acc%value_adc(:,k,tind),&
                                 stats%vid_adc_total(tind,LVT_ACCid,1))
                         enddo
                      endif
                      
                      call LVT_writeSummaryStats(&
                           LVT_metrics%acc%ftn_summ,&
                           LVT_metrics%acc%short_name,&
                           LVT_rc%ngrid,&
                           stats%acc%value_total(:,k,:), &
                           stats%acc%count_value_total(:,k,:),stats%standard_name,&
                           stats%acc%value_ci(k,:))
                   enddo
                endif
             endif
          endif
       end if
    endif
  end subroutine LVT_writeMetric_ACC

!BOP
! 
! !ROUTINE: LVT_resetMetric_ACC
! \label{LVT_resetMetric_ACC}
!
! !INTERFACE: 
  subroutine LVT_resetMetric_ACC(alarm)
! 
! !USES:
    use LVT_coreMod,   only : LVT_rc
!
! !INPUT PARAMETERS: 
    logical,    intent(in) :: alarm
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
    integer                :: i,k,l,index

    type(LVT_metadataEntry), pointer :: model
    type(LVT_statsEntry)   , pointer :: stats

    
    call LVT_getDataStream1Ptr(model)
    call LVT_getstatsEntryPtr(stats)

    do while(associated(model))
       if(stats%selectOpt.eq.1) then 
          do k=1,model%selectNlevs
             if(LVT_metrics%acc%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   stats%acc%value_ts(:,k,l) = 0.0
                   stats%acc%count_value_ts(:,k,l)=0 
                   if(alarm) then 
                      stats%acc%tavg_value_ts(:,k,l) = 0.0
                      stats%acc%tavg_count_value_ts(:,k,l)=0 
                   endif
                enddo
             endif
          enddo
       endif
       model => model%next
       stats => stats%next
    enddo

  end subroutine LVT_resetMetric_ACC


!BOP
! 
! !ROUTINE: LVT_writerestart_ACC
! 
! !INTERFACE:
  subroutine LVT_writerestart_ACC(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for ACC metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%acc%selectOpt.eq.1) then 
       
       print*, 'The writerestart method is not implemented for ACC'
       stop
    end if
    
  end subroutine LVT_writerestart_ACC

!BOP
! 
! !ROUTINE: LVT_readrestart_ACC
! 
! !INTERFACE:
  subroutine LVT_readrestart_ACC(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for ACC metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%acc%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for ACC'
       stop
    end if
    
  end subroutine LVT_readrestart_ACC
end module LVT_ACCMod
