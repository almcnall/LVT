!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_ensLLMod
! label(LVT_ensLLMod)
!
! !INTERFACE:
module LVT_ensLLMod
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

  PRIVATE
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
!  !DESCRIPTION: 
!   This module handles the computations required to compute 
!   ensemble likelihood values of 
!   of desired variables from the LIS output. 
!   
!   NOTES: 
!   * The LIS output should be written in a tile space format
!  
! !FILES USED:
!
!  !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initensLL
  public :: LVT_diagnoseensLL
  public :: LVT_computeensLL
  public :: LVT_writeMetric_ensLL
  public :: LVT_resetMetric_ensLL
  public :: LVT_writerestart_ensLL
  public :: LVT_readrestart_ensLL

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: LVT_ensLL_nbins 
!EOP
  integer, parameter :: LVT_ensLL_nbins = 5000

contains
  subroutine LVT_initensLL(model, obs, stats,metric)
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
! 
!EOP
    if(metric%selectOpt.eq.1) then 
       allocate(stats%ensll%value_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensll%value_total = 0.0
       allocate(stats%ensll%count_value_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensll%count_value_total = 0
       allocate(stats%ensll%value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%ensll%value_ci = LVT_rc%udef
       
       if(metric%computeSC.eq.1) then 
          allocate(stats%ensll%value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%ensll%value_asc = 0.0
          allocate(stats%ensll%count_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%ensll%count_value_asc = 0
       endif
       if(metric%computeADC.eq.1) then 
          allocate(stats%ensll%value_adc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nadc))
          stats%ensll%value_adc = 0.0
          allocate(stats%ensll%count_value_adc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nadc))
          stats%ensll%count_value_adc = 0
       endif
    endif

    if(metric%timeOpt.eq.1) then 
       allocate(stats%ensll%value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensll%value_ts = 0.0
       allocate(stats%ensll%count_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensll%count_value_ts = 0

       allocate(stats%ensll%tavg_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensll%tavg_value_ts = 0.0
       allocate(stats%ensll%tavg_count_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensll%tavg_count_value_ts = 0
       
    endif

!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1    
    metric%obsData = .true. 
    metric%stdevFlag = .false. 
  end subroutine LVT_initensLL

!BOP
! 
! !ROUTINE: LVT_diagnoseensLL
! \label{LVT_diagnoseensLL}
!
! !INTERFACE: 
  subroutine LVT_diagnoseensLL(pass)
! 
! !USES:     

    implicit none

    integer       :: pass
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the computations for 
!   calculating the std of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelensLL](\ref{diagnoseSingleModelensLL})
!     updates the std computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%ensll%selectOpt.eq.1.or.&
            LVT_metrics%ensll%timeOpt.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call diagnoseSingleensLL(model, obs, stats, &
                  LVT_metrics%ensll)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseensLL

!BOP
! 
! !ROUTINE: diagnoseSingleensLL
! \label{diagnoseSingleensLL}
!
! !INTERFACE: 
  subroutine diagnoseSingleensLL(model, obs, stats,metric)
! 
! !USES: 
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the std computation of the 
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
    integer       :: t,k,tind,g,m
    real          :: ll,binsize
    real          :: bins(LVT_ensLL_nbins),sum_bins
    integer       :: binval,obs_binval
    integer       :: nval 
    real          :: sx, sxx, mean_v, std_v

    LVT_rc%ensLLType = 1

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1.and.&
         obs%selectNlevs.ge.1) then       
       do g=1,LVT_rc%ngrid
          do k=1, model%selectNlevs
             if(trim(obs%units).eq.trim(model%units)) then 
                if(metric%selectOpt.eq.1) then 
                   if(LVT_rc%ensLLType.eq.1) then 
                      ll = 0 
                      sx = 0 
                      sxx = 0 
                      nval = 0 
                      do m=1,LVT_LIS_rc(1)%nensem
                         t = (g-1)*LVT_LIS_rc(1)%nensem+m
                         if(model%count(t,k).ne.0.and.obs%count(g,k).ne.0)  then 
                            sx = sx + model%value(t,k)
                            sxx = sxx + model%value(t,k)**2                            
                            nval = nval + 1
                         endif
                      enddo
                      if(nval.gt.0) then 
                         mean_v = sx/nval
                         std_v  = sxx/nval - mean_v**2
                         if(std_v.ge.0) then 
                            std_v = sqrt(std_v)
                            ll = exp(-(obs%count(g,k)-mean_v)**2/(2*std_v**2))/&
                                 sqrt(2*3.14*std_v**2)
                         else
                            ll = LVT_rc%udef
                         endif
                      else
                         ll = LVT_rc%udef
                      endif
                   elseif(LVT_rc%ensLLType.eq.2) then 
                      ll = 0 
                      binsize =(model%valid_max(1) - model%valid_min(1))/LVT_ensLL_nbins
                      bins = 0 
                      nval = 0 
                      do m=1,LVT_LIS_rc(1)%nensem
                         t = (g-1)*LVT_LIS_rc(1)%nensem+m
                         if(model%count(t,k).ne.0.and.obs%count(g,k).ne.0)  then 
                            binval =nint((model%value(t,k)-model%valid_min(1))/binsize)+1
                            if(binval.lt.1) binval = 1
                            if(binval.gt.LVT_ensLL_nbins) binval = LVT_ensLL_nbins
                            bins(binval) = bins(binval) +1 
                            nval = nval+1
                         endif
                      enddo
                      if(nval.gt.0) then 
                         sum_bins = sum(bins)
                         bins = bins/sum_bins
                         obs_binval = nint((obs%value(t,k)-model%valid_min(1))/binsize)+1
                         if(obs_binval.lt.1) obs_binval = 1
                         if(obs_binval.gt.LVT_ensLL_nbins) obs_binval = LVT_ensLL_nbins
                         ll = bins(obs_binval)
                      else
                         ll = LVT_rc%udef
                      endif                   
                   endif
                   if(ll.ne.LVT_rc%udef) then 
                      stats%ensll%value_total(g,k,1) = & 
                           stats%ensll%value_total(g,k,1) + ll
                      stats%ensll%count_value_total(g,k,1) = & 
                           stats%ensll%count_value_total(g,k,1) + 1
                      
                      if(metric%timeOpt.eq.1) then 
                         stats%ensll%value_ts(g,k,1) = & 
                              stats%ensll%value_ts(g,k,1) + ll
                         stats%ensll%count_value_ts(g,k,1) = &
                              stats%ensll%count_value_ts(g,k,1) + 1
                         
                      endif
                      if(metric%computeSC.eq.1) then 
                         call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                              tind)      
                         stats%ensll%value_asc(g,k,tind) = & 
                              stats%ensll%value_asc(g,k,tind) + ll
                         stats%ensll%count_value_asc(g,k,tind) = & 
                              stats%ensll%count_value_asc(g,k,tind) + 1
                      endif
                      if(metric%computeADC.eq.1) then 
                         call LVT_getADCTimeIndex(tind)      
                         stats%ensll%value_adc(g,k,tind) = & 
                              stats%ensll%value_adc(g,k,tind) + ll
                         stats%ensll%count_value_adc(g,k,tind) = & 
                              stats%ensll%count_value_adc(g,k,tind) + 1
                      endif
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(g,k).gt.&
                              LVT_rc%strat_var_threshold) then
                            stats%ensll%value_total(g,k,2) = &
                                 stats%ensll%value_total(g,k,2) + ll
                            stats%ensll%count_value_total(g,k,2) = &
                                 stats%ensll%count_value_total(g,k,2) + 1
                            
                            if(metric%timeOpt.eq.1) then 
                               stats%ensll%value_ts(g,k,2) = & 
                                    stats%ensll%value_ts(g,k,2) + ll
                               stats%ensll%count_value_ts(g,k,2) = &
                                    stats%ensll%count_value_ts(g,k,2) + 1
                               
                            endif
                         elseif(LVT_stats%strat_var(g,k).le.&
                              LVT_rc%strat_var_threshold) then
                            stats%ensll%value_total(g,k,3) = &
                                 stats%ensll%value_total(g,k,3) + ll
                            stats%ensll%count_value_total(g,k,3) = &
                                 stats%ensll%count_value_total(g,k,3) + 1
                               
                            if(metric%timeOpt.eq.1) then 
                               stats%ensll%value_ts(g,k,3) = & 
                                    stats%ensll%value_ts(g,k,3) + ll
                               stats%ensll%count_value_ts(g,k,3) = &
                                    stats%ensll%count_value_ts(g,k,3) + 1
                               
                            endif
                         endif
                      endif
                   endif
                endif
             endif
          enddo
       enddo
    endif
  end subroutine diagnoseSingleensLL

!BOP
! 
! !ROUTINE: LVT_computeensLL
! \label{LVT_computeensLL}
!
! !INTERFACE: 
  subroutine LVT_computeensLL(pass,alarm)
! 
! !USES: 

    implicit none

    integer               :: pass
    logical               :: alarm
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the std values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelensLL](\ref{computeSingleModelensLL})
!     computes the std values for a single variable
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer     :: i 
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats
    
    if(pass.eq.1) then 
       if(LVT_metrics%ensll%selectOpt.eq.1.or.&
            LVT_metrics%ensll%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%ensll%timeOpt.eq.1.and.&
                  LVT_metrics%ensll%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ensll%ftn_ts_loc(i),200,advance='no') &
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
             call computeSingleensLL(alarm,&
                  model, obs, stats,LVT_metrics%ensll)

             model => model%next
             obs => obs%next
             stats => stats%next             
          enddo
          
          if(alarm) then 
             if(LVT_metrics%ensll%timeOpt.eq.1.and.&
                  LVT_metrics%ensll%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ensll%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_computeensLL
  

!BOP
! 
! !ROUTINE: computeSingleensLL
! \label{computeSingleensLL}
!
! !INTERFACE: 
  subroutine computeSingleensLL(alarm,model,obs,stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the std values
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

    if(metric%timeOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%ensll%count_value_ts(t,k,l).gt.0) then 
                      stats%ensll%value_ts(t,k,l) = stats%ensll%value_ts(t,k,l)/&
                           stats%ensll%count_value_ts(t,k,l)                   

                      stats%ensll%tavg_value_ts(t,k,l) = &
                           stats%ensll%tavg_value_ts(t,k,l)  + & 
                           stats%ensll%value_ts(t,k,l)
                      stats%ensll%tavg_count_value_ts(t,k,l) = &
                           stats%ensll%tavg_count_value_ts(t,k,l)  + 1

                   else
                      stats%ensll%value_ts(t,k,l) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
          
          if(alarm) then 

             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%ensll%tavg_count_value_ts(t,k,l).gt.0) then 
                         stats%ensll%tavg_value_ts(t,k,l) = &
                              stats%ensll%tavg_value_ts(t,k,l) / & 
                              stats%ensll%tavg_count_value_ts(t,k,l) 
                      endif
                   enddo
                enddo
             enddo
             if(metric%extractTS.eq.1) then 
                call LVT_writeTSinfo(metric%ftn_ts_loc,&
                     model,&
                     LVT_rc%ngrid,&
                     stats%ensll%tavg_value_ts,&
                     stats%ensll%tavg_count_value_ts)
             endif
          endif
       endif
    endif
       
    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%ensll%count_value_total(t,k,l).gt.&
                        LVT_rc%obsCountThreshold) then 
                      stats%ensll%value_total(t,k,l) = &
                           stats%ensll%value_total(t,k,l)/&
                              stats%ensll%count_value_total(t,k,l)           
                   else
                      stats%ensll%value_total(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1, LVT_rc%nasc
                      if(stats%ensll%count_value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then 
                         stats%ensll%value_asc(t,k,l) = &
                              stats%ensll%value_asc(t,k,l)/&
                              stats%ensll%count_value_asc(t,k,l)           
                      else
                         stats%ensll%value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif

                if(metric%computeADC.eq.1) then 
                   do l=1, LVT_rc%nadc
                      if(stats%ensll%count_value_adc(t,k,l).gt.&
                           LVT_rc%ADCCountThreshold) then 
                         stats%ensll%value_adc(t,k,l) = &
                              stats%ensll%value_adc(t,k,l)/&
                              stats%ensll%count_value_adc(t,k,l)           
                      else
                         stats%ensll%value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
             enddo
          enddo
             
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%ensll%value_total(:,k,l),&
                     LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%ensll%value_ci(k,l))
             enddo
          enddo
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid,stats%ensll%value_total)
          
          if(metric%computeSC.eq.1) then 
             call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%ensll%value_asc,stats%ensll%count_value_asc)
          endif
          if(metric%computeADC.eq.1) then 
             call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%ensll%value_adc,stats%ensll%count_value_adc)
          endif
       endif
    endif
  end subroutine computeSingleensLL


  subroutine LVT_writeMetric_ensLL(pass,final,vlevels,stats,obs)
! !USES:
    use LVT_statsMod
    use LVT_pluginIndices
! !ARGUMENTS:
    implicit none

    integer                 :: pass
    integer                 :: final
    integer                 :: vlevels
    type(LVT_statsEntry)    :: stats
    type(LVT_metaDataEntry) :: obs

    real                    :: dummy
    integer                 :: k,l,tind

    if(pass.eq.LVT_metrics%ensll%npass) then 
       if(final.ne.1) then
           if(stats%selectOpt.eq.1) then 
             do k=1,vlevels
                if(LVT_metrics%ensll%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_gridded(LVT_metrics%ensll%ftn_ts, &
                           stats%ensll%value_ts(:,k,l),&
                           stats%vid_ts(LVT_ensLLid,1),1)
                      call LVT_writevar_gridded(LVT_metrics%ensll%ftn_ts, &
                           real(stats%ensll%count_value_ts(:,k,l)),&
                           stats%vid_count_ts(LVT_ensLLid,1),1)
                   enddo
                endif
             enddo
          endif
       else
          if(pass.eq.1) then 
             if(stats%selectOpt.eq.1) then 
                do k=1,vlevels
                   if(LVT_metrics%ensll%selectOpt.eq.1) then 
                      do l=1,LVT_rc%strat_nlevels                      
                         call LVT_writevar_gridded(&
                              LVT_metrics%ensll%ftn_total, &
                              stats%ensll%value_total(:,k,l),&
                              stats%vid_total(LVT_ensLLid,1))
                         call LVT_writevar_gridded(&
                              LVT_metrics%ensll%ftn_total, &
                              real(stats%ensll%count_value_total(:,k,l)),&
                              stats%vid_count_total(LVT_ensLLid,1))
                         
                      enddo
                   
                      if(LVT_metrics%ensll%computeSC.eq.1) then 
                         do tind = 1,LVT_rc%nasc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%ensll%ftn_total,&
                                 stats%ensll%value_asc(:,k,tind),&
                                 stats%vid_sc_total(tind,LVT_ensLLid,1))
                         enddo
                      endif
                      if(LVT_metrics%ensll%computeADC.eq.1) then 
                         do tind = 1,LVT_rc%nadc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%ensll%ftn_total,&
                                 stats%ensll%value_adc(:,k,tind),&
                                 stats%vid_adc_total(tind,LVT_ensLLid,1))
                         enddo
                      endif
                      call LVT_writeSummaryStats(&
                           LVT_metrics%ensll%ftn_summ,&
                           LVT_metrics%ensll%short_name,&
                           LVT_rc%ngrid,&
                           stats%ensll%value_total(:,k,:), &
                           stats%ensll%count_value_total(:,k,:),&
                           stats%standard_name,&
                           stats%ensll%value_ci(k,:))
                   endif
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_writeMetric_ensLL

!BOP
! 
! !ROUTINE: LVT_reset_Metric_ensLL
! \label(LVT_reset_Metric_ensLL)
!
! !INTERFACE:
  subroutine LVT_resetMetric_ensLL(alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
    logical              :: alarm
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

    integer                :: i,k,l
    type(LVT_metadataEntry), pointer :: model
    type(LVT_statsEntry)   , pointer :: stats

    call LVT_getDataStream1Ptr(model)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))

       if(stats%selectOpt.eq.1) then 
          do k=1,model%selectNlevs
             if(LVT_metrics%ensll%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   stats%ensll%value_ts(:,k,l) = 0.0
                   stats%ensll%count_value_ts(:,k,l)=0 
                   if(alarm) then 
                      stats%ensll%tavg_value_ts(:,k,l) = 0.0
                      stats%ensll%tavg_count_value_ts(:,k,l)=0 
                   endif
                enddo
             endif
          enddo
       endif
       
       model => model%next
       stats => stats%next
    enddo

  end subroutine LVT_resetMetric_ensLL


!BOP
! 
! !ROUTINE: LVT_writerestart_ensLL
! 
! !INTERFACE:
  subroutine LVT_writerestart_ensLL(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for ensLL metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ensll%selectOpt.eq.1) then 
       
       print*, 'The writerestart method is not implemented for ensLL'
       stop
    end if
    
  end subroutine LVT_writerestart_ensLL

!BOP
! 
! !ROUTINE: LVT_readrestart_ensLL
! 
! !INTERFACE:
  subroutine LVT_readrestart_ensLL(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for ensLL metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ensll%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for ensLL'
       stop
    end if
    
  end subroutine LVT_readrestart_ensLL

end module LVT_ensLLMod
