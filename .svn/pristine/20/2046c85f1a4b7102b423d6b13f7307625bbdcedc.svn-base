!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_RMSEMod
! \label(LVT_RMSEMod)
!
! !INTERFACE:
module LVT_RMSEMod
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
! !DESCRIPTION: 
!  This module handles the RMSE computations by comparing 
!  the LIS output to the specified observations. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP


!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initRMSE
  public :: LVT_diagnoseRMSE
  public :: LVT_computeRMSE 
  public :: LVT_writeMetric_RMSE
  public :: LVT_resetMetric_RMSE
  public :: LVT_writerestart_RMSE
  public :: LVT_readrestart_RMSE
contains
  
!BOP
! 
! !ROUTINE: LVT_initRMSE
! \label{LVT_initRMSE}
!
! !INTERFACE: 
  subroutine LVT_initRMSE(model,obs,stats,metric)
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
    if(metric%selectOpt.eq.1) then 
       allocate(stats%rmse%value_total(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))
       stats%rmse%value_total = 0.0
       allocate(stats%rmse%count_value_total(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       
       stats%rmse%count_value_total = 0
       allocate(stats%rmse%value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%rmse%value_ci = LVT_rc%udef

       allocate(stats%rmse%value_stdev_total(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))
       stats%rmse%value_stdev_total = 0.0
       allocate(stats%rmse%value_stdev_total_sxx(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))
       stats%rmse%value_stdev_total_sxx = 0.0
       allocate(stats%rmse%count_value_stdev_total(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       
       stats%rmse%count_value_stdev_total = 0
    endif
    
    if(metric%timeopt.eq.1) then 
       allocate(stats%rmse%value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%rmse%value_ts = 0.0
       allocate(stats%rmse%count_value_ts(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       
       stats%rmse%count_value_ts = 0 

       allocate(stats%rmse%tavg_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       allocate(stats%rmse%tavg_count_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%rmse%tavg_value_ts = 0.0
       stats%rmse%tavg_count_value_ts = 0.0
       
       if(metric%computeSC.eq.1) then 
          allocate(stats%rmse%value_asc(LVT_rc%ngrid, model%selectNlevs,LVT_rc%nasc))
          stats%rmse%value_asc = 0.0
          allocate(stats%rmse%count_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%rmse%count_value_asc = 0
       endif
       if(metric%computeADC.eq.1) then 
          allocate(stats%rmse%value_adc(LVT_rc%ngrid, model%selectNlevs,LVT_rc%nadc))
          stats%rmse%value_adc = 0.0
          allocate(stats%rmse%count_value_adc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nadc))
          stats%rmse%count_value_adc = 0
       endif
    endif
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1
    metric%obsData = .false. 
    metric%stdevFlag = .true. 

  end subroutine LVT_initRMSE
  
!BOP
! 
! !ROUTINE: LVT_diagnoseRMSE
! \label{LVT_diagnoseRMSE}
!
! !INTERFACE: 
  subroutine LVT_diagnoseRMSE(pass)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the RMSE calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleRMSE](\ref{diagnoseSingleRMSE})
!     updates the RMSE computation for a single variable 
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
       if(LVT_metrics%rmse%selectOpt.eq.1.or.&
            LVT_metrics%rmse%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleRMSE(obs,&
                  model, stats, &
                  LVT_metrics%rmse)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif

  end subroutine LVT_diagnoseRMSE

!BOP
! 
! !ROUTINE: diagnoseSingleRMSE
! \label{diagnoseSingleRMSE}
!
! !INTERFACE: 
  subroutine diagnoseSingleRMSE(obs, model, stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the RMSE computation (updates the running 
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
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer    :: t,k

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             if(trim(obs%units).eq.trim(model%units)) then
                if(obs%count(t,k).ne.0.and. &
                     model%count(t,k).ne.0) then      
                   if(metric%selectOpt.eq.1) then
                      stats%rmse%value_total(t,k,1) = &
                           stats%rmse%value_total(t,k,1) + &
                           (obs%value(t,k) - model%value(t,k))* & 
                           (obs%value(t,k) - model%value(t,k))
                      stats%rmse%count_value_total(t,k,1) = &
                           stats%rmse%count_value_total(t,k,1) + 1
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then 
                            stats%rmse%value_total(t,k,2) = &
                                 stats%rmse%value_total(t,k,2) + &
                                 (obs%value(t,k) - model%value(t,k))* & 
                                 (obs%value(t,k) - model%value(t,k))
                            stats%rmse%count_value_total(t,k,2) = &
                                 stats%rmse%count_value_total(t,k,2) + 1
                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then 
                            stats%rmse%value_total(t,k,3) = &
                                 stats%rmse%value_total(t,k,3) + &
                                 (obs%value(t,k) - model%value(t,k))* & 
                                 (obs%value(t,k) - model%value(t,k))
                            stats%rmse%count_value_total(t,k,3) = & 
                                 stats%rmse%count_value_total(t,k,3) + 1
                         endif
                      endif
                   endif
                   if(metric%timeOpt.eq.1) then 
                      stats%rmse%value_ts(t,k,1) = &
                           stats%rmse%value_ts(t,k,1) + &
                           (obs%value(t,k) - model%value(t,k))* & 
                           (obs%value(t,k) - model%value(t,k))
                      stats%rmse%count_value_ts(t,k,1) = &
                           stats%rmse%count_value_ts(t,k,1)+1
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then
                            stats%rmse%value_ts(t,k,2) = &
                                 stats%rmse%value_ts(t,k,2) + &
                                 (obs%value(t,k) - model%value(t,k))* & 
                                 (obs%value(t,k) - model%value(t,k))
                            stats%rmse%count_value_ts(t,k,2) = &
                                 stats%rmse%count_value_ts(t,k,2)+1
                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then 
                            stats%rmse%value_ts(t,k,3) = &
                                 stats%rmse%value_ts(t,k,3) + &
                                 (obs%value(t,k) - model%value(t,k))* & 
                                 (obs%value(t,k) - model%value(t,k))
                            stats%rmse%count_value_ts(t,k,3) = &
                                 stats%rmse%count_value_ts(t,k,3)+1
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
    
  end subroutine diagnoseSingleRMSE


!BOP
! 
! !ROUTINE: LVT_computeRMSE
! \label{LVT_computeRMSE}
!
! !INTERFACE: 
  subroutine LVT_computeRMSE(pass,alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute RMSE values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleRMSE](\ref{computeSingleRMSE})
!     updates the RMSE computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     RMSE computation has been reached
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
       if(LVT_metrics%rmse%selectOpt.eq.1.or.&
            LVT_metrics%rmse%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%rmse%timeOpt.eq.1.and.&
                  LVT_metrics%rmse%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%rmse%ftn_ts_loc(i),200,advance='no') &
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

             call computeSingleRMSE(alarm,obs,&
                  model, stats, LVT_metrics%rmse)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%rmse%timeOpt.eq.1.and.&
                  LVT_metrics%rmse%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%rmse%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif

  end subroutine LVT_ComputeRMSE

!BOP
! 
! !ROUTINE: computeSingleRMSE
! \label{computeSingleRMSE}
!
! !INTERFACE: 
  subroutine computeSingleRMSE(alarm,obs, model,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the RMSE values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     RMSE computation has been reached
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
                   if(stats%rmse%count_value_ts(t,k,l).ne.0) then 
                      stats%rmse%value_ts(t,k,l) = sqrt(stats%rmse%value_ts(t,k,l)/&
                           stats%rmse%count_value_ts(t,k,l))

                      stats%rmse%tavg_value_ts(t,k,l) = stats%rmse%tavg_value_ts(t,k,l) + &
                           stats%rmse%value_ts(t,k,l)    
                      stats%rmse%tavg_count_value_ts(t,k,l) = &
                           stats%rmse%tavg_count_value_ts(t,k,l) + 1

                   else
                      stats%rmse%value_ts(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   if(stats%rmse%count_value_ts(t,k,1).ne.0) then 
                      call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                           tind)
                      stats%rmse%value_asc(t,k,tind) = &
                           stats%rmse%value_asc(t,k,tind)+ &
                           stats%rmse%value_ts(t,k,1)
                      stats%rmse%count_value_asc(t,k,tind) = &
                           stats%rmse%count_value_asc(t,k,tind)+ 1
                   endif
                endif
                if(metric%computeADC.eq.1) then 
                   if(stats%rmse%count_value_ts(t,k,1).ne.0) then 
                      call LVT_getADCTimeIndex(tind)
                      stats%rmse%value_adc(t,k,tind) = &
                           stats%rmse%value_adc(t,k,tind)+ &
                           stats%rmse%value_ts(t,k,1)
                      stats%rmse%count_value_adc(t,k,tind) = &
                           stats%rmse%count_value_adc(t,k,tind)+ 1
                   endif                
                endif
             enddo
          enddo
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%rmse%value_ts(t,k,l).ne.LVT_rc%udef) then 
                      stats%rmse%value_stdev_total(t,k,l) = &
                           stats%rmse%value_stdev_total(t,k,l) + &
                           stats%rmse%value_ts(t,k,l)
                      stats%rmse%value_stdev_total_sxx(t,k,l) = &
                           stats%rmse%value_stdev_total_sxx(t,k,l) + &
                           stats%rmse%value_ts(t,k,l)*stats%rmse%value_ts(t,k,l)
                      stats%rmse%count_value_stdev_total(t,k,l) = &
                           stats%rmse%count_value_stdev_total(t,k,l) + 1
                   endif
                enddo
             enddo
          enddo

          if(alarm) then 
             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%rmse%tavg_count_value_ts(t,k,l).gt.0) then 
                         stats%rmse%tavg_value_ts(t,k,l) = stats%rmse%tavg_value_ts(t,k,l) / &
                              stats%rmse%tavg_count_value_ts(t,k,l) 
                      else
                         stats%rmse%tavg_value_ts(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                enddo
             enddo

             if(metric%extractTS.eq.1) then 
                call LVT_writeTSinfo(metric%ftn_ts_loc,&
                     model,&
                     LVT_rc%ngrid,&
                     stats%rmse%tavg_value_ts,&
                     stats%rmse%tavg_count_value_ts)
             endif
          endif
       endif
    endif

    if(alarm) then 
       if(metric%timeOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%rmse%tavg_count_value_ts(t,k,l).ne.0) then 
                         stats%rmse%tavg_value_ts(t,k,l) = &
                              stats%rmse%tavg_value_ts(t,k,l) /&
                              stats%rmse%tavg_count_value_ts(t,k,l)
                      endif
                   enddo
                enddo
             enddo
          endif
       endif
    endif

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%rmse%count_value_total(t,k,l).ne.0.and.&
                        stats%rmse%count_value_total(t,k,l)&
                        .gt.LVT_rc%obsCountThreshold) then 
                      
                      stats%rmse%value_total(t,k,l) = sqrt(stats%rmse%value_total(t,k,l)/&
                           stats%rmse%count_value_total(t,k,l))                      
                   else
                      stats%rmse%value_total(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1,LVT_rc%nasc
                      if(stats%rmse%count_value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then
                         stats%rmse%value_asc(t,k,l) = stats%rmse%value_asc(t,k,l)/&
                              stats%rmse%count_value_asc(t,k,l)
                      else
                         stats%rmse%value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
                
                if(metric%computeADC.eq.1) then 
                   do l=1,LVT_rc%nadc
                      if(stats%rmse%count_value_adc(t,k,l).gt.&
                           LVT_rc%ADCCountThreshold) then
                         stats%rmse%value_adc(t,k,l) = stats%rmse%value_adc(t,k,l)/&
                              stats%rmse%count_value_adc(t,k,l)
                      else
                         stats%rmse%value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
             enddo
          enddo

          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%rmse%count_value_stdev_total(t,k,l).gt.0) then 
                      stats%rmse%value_stdev_total(t,k,l) = &
                           stats%rmse%value_stdev_total(t,k,l) /&
                           stats%rmse%count_value_stdev_total(t,k,l)
                      
                      if(stats%rmse%value_stdev_total_sxx(t,k,l)/&
                           stats%rmse%count_value_stdev_total(t,k,l).gt. &
                           stats%rmse%value_stdev_total(t,k,l)**2) then 
                         stats%rmse%value_stdev_total(t,k,l) = &
                              sqrt(stats%rmse%value_stdev_total_sxx(t,k,l)/&
                              stats%rmse%count_value_stdev_total(t,k,l) - &
                              stats%rmse%value_stdev_total(t,k,l)**2)
                      else
                         stats%rmse%value_stdev_total(t,k,l) = LVT_rc%udef 
                      endif
                   else
                      stats%rmse%value_stdev_total(t,k,l) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo

          
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%rmse%value_total(:,k,l),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%rmse%value_ci(k,l))
             enddo
          enddo
       
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid,stats%rmse%value_total)      
          if(metric%computeSC.eq.1) then
             call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%rmse%value_asc,stats%rmse%count_value_asc)
          endif
          if(metric%computeADC.eq.1) then 
             call LVT_writeAvgDiurnalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%rmse%value_adc,stats%rmse%count_value_adc)
          endif
       endif
    endif


  end subroutine computeSingleRMSE

!BOP
! 
! !ROUTINE: LVT_writeMetric_RMSE
! \label{LVT_writeMetric_RMSE}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_RMSE(pass,final,vlevels,stats,obs)
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
    type(LVT_statsEntry)    :: stats
    type(LVT_metaDataEntry) :: obs
!EOP
    integer                 :: l,tind
    integer                 :: k

    if(pass.eq.LVT_metrics%rmse%npass) then 
       if(final.ne.1) then 
          if(stats%selectOpt.eq.1) then 
             do k=1,vlevels
                if(LVT_metrics%rmse%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_gridded(LVT_metrics%rmse%ftn_ts, &
                           stats%rmse%value_ts(:,k,l),stats%vid_ts(LVT_RMSEid,1),k)
                      call LVT_writevar_gridded(LVT_metrics%rmse%ftn_ts, &
                           real(stats%rmse%count_value_ts(:,k,l)),&
                           stats%vid_count_ts(LVT_RMSEid,1),k)
                   enddo
                endif
             enddo
          endif
       else
          if(pass.eq.1) then 
             if(LVT_metrics%rmse%selectOpt.eq.1) then
                if(stats%selectOpt.eq.1) then 
                   do k=1,vlevels
                      do l=1,LVT_rc%strat_nlevels
                         
                         call LVT_writevar_gridded(LVT_metrics%rmse%ftn_total, &
                              stats%rmse%value_total(:,k,l),&
                              stats%vid_total(LVT_RMSEid,1),k)
                         call LVT_writevar_gridded(LVT_metrics%rmse%ftn_total, &
                              real(stats%rmse%count_value_total(:,k,l)),&
                              stats%vid_count_total(LVT_RMSEid,1),k)
                         
                         call LVT_writevar_gridded(LVT_metrics%rmse%ftn_total, &
                              stats%rmse%value_stdev_total(:,k,l),&
                              stats%vid_stdev_total(LVT_RMSEid,1),k)
                         call LVT_writevar_gridded(LVT_metrics%rmse%ftn_total, &
                              real(stats%rmse%count_value_stdev_total(:,k,l)),&
                              stats%vid_count_stdev_total(LVT_RMSEid,1),k)
                      enddo

                      if(LVT_metrics%rmse%computeSC.eq.1) then 
                         do tind = 1,LVT_rc%nasc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%rmse%ftn_total,&
                                 stats%rmse%value_asc(:,k,tind),&
                                 stats%vid_sc_total(tind,LVT_RMSEid,1),k)
                         enddo
                      endif
                      
                      if(LVT_metrics%rmse%computeADC.eq.1) then 
                         do tind = 1,LVT_rc%nadc
                            call LVT_writevar_gridded(LVT_metrics%rmse%ftn_total,&
                                 stats%rmse%value_adc(:,k,tind),&
                                 stats%vid_adc_total(tind,LVT_RMSEid,1),k)
                         enddo
                      endif
                      call LVT_writeSummaryStats(&
                           LVT_metrics%rmse%ftn_summ,&
                           LVT_metrics%rmse%short_name,&
                           LVT_rc%ngrid,&
                           stats%rmse%value_total(:,k,:), &
                           stats%rmse%count_value_total(:,k,:),stats%standard_name,&
                           stats%rmse%value_ci(k,:))
                   enddo
                endif
             endif
          endif
       end if
    endif
  end subroutine LVT_writeMetric_RMSE

!BOP
! 
! !ROUTINE: LVT_resetMetric_RMSE
! \label(LVT_resetMetric_RMSE)
!
! !INTERFACE:
  subroutine LVT_resetMetric_RMSE(alarm)
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

    integer                :: i,k,l
    type(LVT_metadataEntry), pointer :: model
    type(LVT_statsEntry)   , pointer :: stats


    call LVT_getDataStream1Ptr(model)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       
       if(stats%selectOpt.eq.1) then 
          do k=1,model%selectNlevs
             if(LVT_metrics%rmse%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   stats%rmse%value_ts(:,k,l) = 0.0
                   stats%rmse%count_value_ts(:,k,l)=0 
                   if(alarm) then 
                      stats%rmse%tavg_value_ts(:,k,l) = 0.0
                      stats%rmse%tavg_count_value_ts(:,k,l)=0 
                   endif
                enddo
             endif
             
          enddo
          
       endif
       model => model%next
       stats => stats%next

    enddo
    

  end subroutine LVT_resetMetric_RMSE

!BOP
! 
! !ROUTINE: LVT_writerestart_RMSE
! 
! !INTERFACE:
  subroutine LVT_writerestart_RMSE(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for RMSE metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP

    integer              :: k,l
    type(LVT_metaDataEntry), pointer :: model
    type(LVT_metaDataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats

    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       
       if(LVT_metrics%rmse%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   
                   call LVT_writevar_restart(ftn,&
                        stats%rmse%value_total(:,k,l))
                   call LVT_writevar_restart(ftn,&
                        stats%rmse%count_value_total(:,k,l))
                enddo
             enddo
             if(LVT_metrics%rmse%computeSC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc
                      call LVT_writevar_restart(ftn,&
                           stats%rmse%value_asc(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%rmse%count_value_asc(:,k,l))
                   enddo
                enddo
             endif
             if(LVT_metrics%rmse%computeADC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nadc
                      call LVT_writevar_restart(ftn,&
                           stats%rmse%value_adc(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%rmse%count_value_adc(:,k,l))
                   enddo
                enddo
             endif
          end if
       end if

       model => model%next
       obs   => obs%next
       stats => stats%next

    end do
             
  end subroutine LVT_writerestart_RMSE

!BOP
! 
! !ROUTINE: LVT_readrestart_RMSE
! 
! !INTERFACE:
  subroutine LVT_readrestart_RMSE(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for RMSE metric computations
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


    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       
       if(LVT_metrics%rmse%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   
                   call LVT_readvar_restart(ftn,&
                        stats%rmse%value_total(:,k,l))
                   call LVT_readvar_restart(ftn,&
                        stats%rmse%count_value_total(:,k,l))
                enddo
             enddo
             if(LVT_metrics%rmse%computeSC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc
                      call LVT_readvar_restart(ftn,&
                           stats%rmse%value_asc(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%rmse%count_value_asc(:,k,l))
                   enddo
                enddo
             endif
             if(LVT_metrics%rmse%computeADC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nadc
                      call LVT_readvar_restart(ftn,&
                           stats%rmse%value_adc(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%rmse%count_value_adc(:,k,l))
                   enddo
                enddo
             endif
          end if
       endif

       model => model%next
       obs   => obs%next
       stats => stats%next

    enddo

             
  end subroutine LVT_readrestart_RMSE

end module LVT_RMSEMod
