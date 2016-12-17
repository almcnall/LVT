!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_ubRMSEMod
! \label(LVT_ubRMSEMod)
!
! !INTERFACE:
module LVT_ubRMSEMod
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
!  This module handles the unbiased RMSE (ubRMSE) computations by comparing 
!  the LIS output to the specified observations. 
!  
!   ubRMSE = sqrt(RMSE**2 - bias**2)
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
  public :: LVT_initubRMSE
  public :: LVT_diagnoseubRMSE
  public :: LVT_computeubRMSE 
  public :: LVT_writeMetric_ubRMSE
  public :: LVT_resetMetric_ubRMSE
  public :: LVT_writerestart_ubRMSE
  public :: LVT_readrestart_ubRMSE
  
contains
  
!BOP
! 
! !ROUTINE: LVT_initubRMSE
! \label{LVT_initubRMSE}
!
! !INTERFACE: 
  subroutine LVT_initubRMSE(model,obs,stats,metric)
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
       allocate(stats%ubrmse%value_total(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))
       allocate(stats%ubrmse%value_b2_total(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))
       allocate(stats%ubrmse%value_b1_total(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))

       stats%ubrmse%value_total = 0.0
       stats%ubrmse%value_b2_total = 0.0
       stats%ubrmse%value_b1_total = 0.0

       allocate(stats%ubrmse%count_value_total(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       
       stats%ubrmse%count_value_total = 0

       allocate(stats%ubrmse%value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%ubrmse%value_ci = LVT_rc%udef
    endif
    
    if(metric%timeopt.eq.1) then 
       allocate(stats%ubrmse%value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       allocate(stats%ubrmse%value_b1_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       allocate(stats%ubrmse%value_b2_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ubrmse%value_ts = 0.0
       stats%ubrmse%value_b1_ts = 0.0
       stats%ubrmse%value_b2_ts = 0.0
       allocate(stats%ubrmse%count_value_ts(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       
       stats%ubrmse%count_value_ts = 0 

       allocate(stats%ubrmse%tavg_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ubrmse%tavg_value_ts = 0.0
       allocate(stats%ubrmse%tavg_count_value_ts(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       
       stats%ubrmse%tavg_count_value_ts = 0 

       if(metric%computeSC.eq.1) then 
          allocate(stats%ubrmse%value_asc(LVT_rc%ngrid, model%selectNlevs,LVT_rc%nasc))
          stats%ubrmse%value_asc = 0.0
          allocate(stats%ubrmse%count_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%ubrmse%count_value_asc = 0
       endif
       if(metric%computeADC.eq.1) then 
          allocate(stats%ubrmse%value_adc(LVT_rc%ngrid, model%selectNlevs,LVT_rc%nadc))
          stats%ubrmse%value_adc = 0.0
          allocate(stats%ubrmse%count_value_adc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nadc))
          stats%ubrmse%count_value_adc = 0
       endif
    endif
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1
    metric%obsData = .false. 

  end subroutine LVT_initubRMSE
  
!BOP
! 
! !ROUTINE: LVT_diagnoseubRMSE
! \label{LVT_diagnoseubRMSE}
!
! !INTERFACE: 
  subroutine LVT_diagnoseubRMSE(pass)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the ubRMSE calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleubRMSE](\ref{diagnoseSingleubRMSE})
!     updates the ubRMSE computation for a single variable 
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
       if(LVT_metrics%ubrmse%selectOpt.eq.1.or.&
            LVT_metrics%ubrmse%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleubRMSE(obs,model,stats,&
                  LVT_metrics%ubrmse)

             model => model%next
             obs => obs%next
             stats => stats%next
             
          enddo
       endif
    endif
  end subroutine LVT_diagnoseubRMSE

!BOP
! 
! !ROUTINE: diagnoseSingleubRMSE
! \label{diagnoseSingleubRMSE}
!
! !INTERFACE: 
  subroutine diagnoseSingleubRMSE(obs, model, stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the ubRMSE computation (updates the running 
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
             if(trim(obs%units).eq.trim(model%units)) then 
                if(obs%count(t,k).ne.0.and. &
                     model%count(t,k).ne.0) then      
                   if(metric%selectOpt.eq.1) then
                      stats%ubrmse%value_b2_total(t,k,1) = &
                           stats%ubrmse%value_b2_total(t,k,1) + &
                           (obs%value(t,k) - model%value(t,k))* & 
                           (obs%value(t,k) - model%value(t,k))
                       stats%ubrmse%value_b1_total(t,k,1) = &
                           stats%ubrmse%value_b1_total(t,k,1) + &
                           (obs%value(t,k) - model%value(t,k))
                      stats%ubrmse%count_value_total(t,k,1) = &
                           stats%ubrmse%count_value_total(t,k,1) + 1
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then 
                            stats%ubrmse%value_b2_total(t,k,2) = &
                                 stats%ubrmse%value_b2_total(t,k,2) + &
                                 (obs%value(t,k) - model%value(t,k))* & 
                                 (obs%value(t,k) - model%value(t,k))
                            stats%ubrmse%value_b1_total(t,k,2) = &
                                 stats%ubrmse%value_b1_total(t,k,2) + &
                                 (obs%value(t,k) - model%value(t,k))
                            stats%ubrmse%count_value_total(t,k,2) = &
                                 stats%ubrmse%count_value_total(t,k,2) + 1
                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then 
                            stats%ubrmse%value_b2_total(t,k,3) = &
                                 stats%ubrmse%value_b2_total(t,k,3) + &
                                 (obs%value(t,k) - model%value(t,k))* & 
                                 (obs%value(t,k) - model%value(t,k))
                            stats%ubrmse%value_b1_total(t,k,3) = &
                                 stats%ubrmse%value_b1_total(t,k,3) + &
                                 (obs%value(t,k) - model%value(t,k))
                            stats%ubrmse%count_value_total(t,k,3) = & 
                                 stats%ubrmse%count_value_total(t,k,3) + 1
                         endif
                      endif
                   endif
                   if(metric%timeOpt.eq.1) then 
                      stats%ubrmse%value_b2_ts(t,k,1) = stats%ubrmse%value_b2_ts(t,k,1) + &
                           (obs%value(t,k) - model%value(t,k))* & 
                           (obs%value(t,k) - model%value(t,k))
                      stats%ubrmse%value_b1_ts(t,k,1) = stats%ubrmse%value_b1_ts(t,k,1) + &
                           (obs%value(t,k) - model%value(t,k))
                      stats%ubrmse%count_value_ts(t,k,1) = &
                           stats%ubrmse%count_value_ts(t,k,1)+1
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then
                            stats%ubrmse%value_b2_ts(t,k,2) = &
                                 stats%ubrmse%value_b2_ts(t,k,2) + &
                                 (obs%value(t,k) - model%value(t,k))* & 
                                 (obs%value(t,k) - model%value(t,k))
                            stats%ubrmse%value_b1_ts(t,k,2) = &
                                 stats%ubrmse%value_b1_ts(t,k,2) + &
                                 (obs%value(t,k) - model%value(t,k))
                            stats%ubrmse%count_value_ts(t,k,2) = &
                                 stats%ubrmse%count_value_ts(t,k,2)+1
                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then 
                            stats%ubrmse%value_b2_ts(t,k,3) = &
                                 stats%ubrmse%value_b2_ts(t,k,3) + &
                                 (obs%value(t,k) - model%value(t,k))* & 
                                 (obs%value(t,k) - model%value(t,k))
                            stats%ubrmse%value_b1_ts(t,k,3) = &
                                 stats%ubrmse%value_b1_ts(t,k,3) + &
                                 (obs%value(t,k) - model%value(t,k))
                            stats%ubrmse%count_value_ts(t,k,3) = &
                                 stats%ubrmse%count_value_ts(t,k,3)+1
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
    
  end subroutine diagnoseSingleubRMSE


!BOP
! 
! !ROUTINE: LVT_computeubRMSE
! \label{LVT_computeubRMSE}
!
! !INTERFACE: 
  subroutine LVT_computeubRMSE(pass,alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute ubRMSE values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleubRMSE](\ref{computeSingleubRMSE})
!     updates the ubRMSE computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     ubRMSE computation has been reached
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none

    integer               :: pass
    logical     :: alarm
    integer     :: i
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%ubrmse%selectOpt.eq.1.or.&
            LVT_metrics%ubrmse%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%ubrmse%timeOpt.eq.1.and.&
                  LVT_metrics%ubrmse%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ubrmse%ftn_ts_loc(i),200,advance='no') &
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
             
             call computeSingleubRMSE(alarm,&
                  obs,model,stats,LVT_metrics%ubrmse)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%ubrmse%timeOpt.eq.1.and.&
                  LVT_metrics%ubrmse%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ubrmse%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif

  end subroutine LVT_ComputeubRMSE

!BOP
! 
! !ROUTINE: computeSingleubRMSE
! \label{computeSingleubRMSE}
!
! !INTERFACE: 
  subroutine computeSingleubRMSE(alarm,obs, model,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the ubRMSE values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     ubRMSE computation has been reached
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

    real     :: rmse_val, bias_val
    integer  :: t,l,k,tind

    if(metric%timeOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%ubrmse%count_value_ts(t,k,l).ne.0) then 
                      rmse_val = sqrt(stats%ubrmse%value_b2_ts(t,k,l)/&
                           stats%ubrmse%count_value_ts(t,k,l))
                      bias_val = stats%ubrmse%value_b1_ts(t,k,l)/&
                           stats%ubrmse%count_value_ts(t,k,l)
                      stats%ubrmse%value_ts(t,k,l) = &
                           sqrt(rmse_val**2-bias_val**2)

                      stats%ubrmse%tavg_value_ts(t,k,l) = & 
                           stats%ubrmse%tavg_value_ts(t,k,l) + & 
                           stats%ubrmse%value_ts(t,k,l)

                      stats%ubrmse%tavg_count_value_ts(t,k,l) = & 
                           stats%ubrmse%tavg_count_value_ts(t,k,l) + 1
                   else
                      stats%ubrmse%value_ts(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   if(stats%ubrmse%count_value_ts(t,k,1).ne.0) then 
                      call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                           tind)
                      stats%ubrmse%value_asc(t,k,tind) = &
                           stats%ubrmse%value_asc(t,k,tind)+ &
                           stats%ubrmse%value_ts(t,k,1)
                      stats%ubrmse%count_value_asc(t,k,tind) = &
                           stats%ubrmse%count_value_asc(t,k,tind)+ 1
                   endif
                endif
                if(metric%computeADC.eq.1) then 
                   if(stats%ubrmse%count_value_ts(t,k,1).ne.0) then 
                      call LVT_getADCTimeIndex(tind)
                      stats%ubrmse%value_adc(t,k,tind) = &
                           stats%ubrmse%value_adc(t,k,tind)+ &
                           stats%ubrmse%value_ts(t,k,1)
                      stats%ubrmse%count_value_adc(t,k,tind) = &
                           stats%ubrmse%count_value_adc(t,k,tind)+ 1
                   endif                
                endif
             enddo
          enddo
          if(alarm) then 
             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%ubrmse%tavg_count_value_ts(t,k,l).gt.0) then 
                         stats%ubrmse%tavg_value_ts(t,k,l) = & 
                              stats%ubrmse%tavg_value_ts(t,k,l)/ & 
                              stats%ubrmse%tavg_count_value_ts(t,k,l)
                      endif
                   enddo
                enddo
             enddo
             if(metric%extractTS.eq.1) then 
                call LVT_writeTSinfo(metric%ftn_ts_loc,&
                     model,&
                     LVT_rc%ngrid,&
                     stats%ubrmse%tavg_value_ts,&
                     stats%ubrmse%tavg_count_value_ts)
             endif
          endif
       endif
    endif


    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%ubrmse%count_value_total(t,k,l).ne.0.and.&
                        stats%ubrmse%count_value_total(t,k,l)&
                        .gt.LVT_rc%obsCountThreshold) then 
                      rmse_val = sqrt(stats%ubrmse%value_b2_total(t,k,l)/&
                           stats%ubrmse%count_value_total(t,k,l))
                      bias_val = stats%ubrmse%value_b1_total(t,k,l)/&
                           stats%ubrmse%count_value_total(t,k,l)
                      stats%ubrmse%value_total(t,k,l) = &
                           sqrt(rmse_val**2-bias_val**2)
                   else
                      stats%ubrmse%value_total(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1,LVT_rc%nasc
                      if(stats%ubrmse%count_value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then
                         stats%ubrmse%value_asc(t,k,l) = &
                              stats%ubrmse%value_asc(t,k,l)/&
                              stats%ubrmse%count_value_asc(t,k,l)
                      else
                         stats%ubrmse%value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
                
                if(metric%computeADC.eq.1) then 
                   do l=1,LVT_rc%nadc
                      if(stats%ubrmse%count_value_adc(t,k,l).gt.&
                           LVT_rc%ADCCountThreshold) then
                         stats%ubrmse%value_adc(t,k,l) = &
                         stats%ubrmse%value_adc(t,k,l)/&
                              stats%ubrmse%count_value_adc(t,k,l)
                      else
                         stats%ubrmse%value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
             enddo
          enddo
          
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%ubrmse%value_total(:,k,l),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%ubrmse%value_ci(k,l))
             enddo
          enddo
       
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid,stats%ubrmse%value_total)      
          if(metric%computeSC.eq.1) then
             call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%ubrmse%value_asc,stats%ubrmse%count_value_asc)
          endif
          if(metric%computeADC.eq.1) then 
             call LVT_writeAvgDiurnalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%ubrmse%value_adc,stats%ubrmse%count_value_adc)
          endif
       endif       
    endif


  end subroutine computeSingleubRMSE

!BOP
! 
! !ROUTINE: LVT_writeMetric_ubRMSE
! \label{LVT_writeMetric_ubRMSE}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_ubRMSE(pass,final,vlevels,stats,obs)
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

    if(pass.eq.LVT_metrics%ubrmse%npass) then 
       if(final.ne.1) then 
          if(stats%selectOpt.eq.1) then 
             do k=1,vlevels
                if(LVT_metrics%ubrmse%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_gridded(LVT_metrics%ubrmse%ftn_ts, &
                           stats%ubrmse%value_ts(:,k,l),stats%vid_ts(LVT_ubRMSEid,1),1)
                      call LVT_writevar_gridded(LVT_metrics%ubrmse%ftn_ts, &
                           real(stats%ubrmse%count_value_ts(:,k,l)),&
                           stats%vid_count_ts(LVT_ubRMSEid,1),1)
                   enddo
                endif
             enddo
          endif
       else
          if(pass.eq.1) then 
             if(LVT_metrics%ubrmse%selectOpt.eq.1) then
                if(stats%selectOpt.eq.1) then 
                   do k=1,vlevels
                      do l=1,LVT_rc%strat_nlevels
                         
                         call LVT_writevar_gridded(LVT_metrics%ubrmse%ftn_total, &
                              stats%ubrmse%value_total(:,k,l),stats%vid_total(LVT_ubRMSEid,1))
                         call LVT_writevar_gridded(LVT_metrics%ubrmse%ftn_total, &
                              real(stats%ubrmse%count_value_total(:,k,l)),&
                              stats%vid_count_total(LVT_ubRMSEid,1))
                         
                      enddo
                      if(LVT_metrics%ubrmse%computeSC.eq.1) then 
                         do tind = 1,LVT_rc%nasc
                            call LVT_writevar_gridded(LVT_metrics%ubrmse%ftn_total,&
                                 stats%ubrmse%value_asc(:,k,tind),&
                                 stats%vid_sc_total(tind,LVT_ubRMSEid,1))
                         enddo
                      endif
                      
                      if(LVT_metrics%ubrmse%computeADC.eq.1) then 
                         do tind = 1,LVT_rc%nadc
                            call LVT_writevar_gridded(LVT_metrics%ubrmse%ftn_total,&
                                 stats%ubrmse%value_adc(:,k,tind),&
                                 stats%vid_adc_total(tind,LVT_ubRMSEid,1))
                         enddo
                      endif
                      call LVT_writeSummaryStats(&
                           LVT_metrics%ubrmse%ftn_summ,&
                           LVT_metrics%ubrmse%short_name,&
                           LVT_rc%ngrid,&
                           stats%ubrmse%value_total(:,k,:), &
                           stats%ubrmse%count_value_total(:,k,:),stats%standard_name,&
                           stats%ubrmse%value_ci(k,:))
                   enddo
                endif
             endif
          endif
       end if
    endif
  end subroutine LVT_writeMetric_ubRMSE

!BOP
! 
! !ROUTINE: LVT_resetMetric_ubRMSE
! \label(LVT_resetMetric_ubRMSE)
!
! !INTERFACE:
  subroutine LVT_resetMetric_ubRMSE(alarm)
! 
! !INPUT PARAMETERS: 
    logical            :: alarm
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
             if(LVT_metrics%ubrmse%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   stats%ubrmse%value_ts(:,k,l) = 0.0
                   stats%ubrmse%value_b1_ts(:,k,l) = 0.0
                   stats%ubrmse%value_b2_ts(:,k,l) = 0.0
                   stats%ubrmse%count_value_ts(:,k,l)=0 
                   if(alarm) then 
                      stats%ubrmse%tavg_value_ts(:,k,l) = 0.0
                      stats%ubrmse%tavg_count_value_ts(:,k,l)=0 
                   endif
                enddo
             endif
          enddo
       endif
       model => model%next
       stats => stats%next
       
    enddo
    
  end subroutine LVT_resetMetric_ubRMSE


!BOP
! 
! !ROUTINE: LVT_writerestart_ubRMSE
! 
! !INTERFACE:
  subroutine LVT_writerestart_ubRMSE(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for ubRMSE metric computations
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

       if(LVT_metrics%ubRMSE%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 

             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   call LVT_writevar_restart(ftn,&
                        stats%ubrmse%value_total(:,k,l))
                   call LVT_writevar_restart(ftn,&
                        stats%ubrmse%value_b2_total(:,k,l))
                   call LVT_writevar_restart(ftn,&
                        stats%ubrmse%value_b1_total(:,k,l))
                   call LVT_writevar_restart(ftn,&
                        stats%ubrmse%count_value_total(:,k,l))
                enddo
             enddo
          endif
          
          if(LVT_metrics%ubRMSE%computeSC.eq.1) then
             do k=1,model%selectNlevs
                do l=1,LVT_rc%nasc          
                   call LVT_writevar_restart(ftn,&
                        stats%ubrmse%value_asc(:,k,l))
                   call LVT_writevar_restart(ftn,&
                        stats%ubrmse%count_value_asc(:,k,l))
                enddo
             enddo
          endif
          if(LVT_metrics%ubRMSE%computeADC.eq.1) then
             do k=1,model%selectNlevs
                do l=1,LVT_rc%nadc           
                   call LVT_writevar_restart(ftn,&
                        stats%ubrmse%value_adc(:,k,l))
                   call LVT_writevar_restart(ftn,&
                        stats%ubrmse%count_value_adc(:,k,l))
                enddo
             enddo
          endif
       endif

       model => model%next
       obs   => obs%next
       stats => stats%next

    enddo
  end subroutine LVT_writerestart_ubRMSE

!BOP
! 
! !ROUTINE: LVT_readrestart_ubRMSE
! 
! !INTERFACE:
  subroutine LVT_readrestart_ubRMSE(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for ubRMSE metric computations
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
       if(LVT_metrics%ubRMSE%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 

             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   call LVT_readvar_restart(ftn,&
                        stats%ubrmse%value_total(:,k,l))
                   call LVT_readvar_restart(ftn,&
                        stats%ubrmse%value_b2_total(:,k,l))
                   call LVT_readvar_restart(ftn,&
                        stats%ubrmse%value_b1_total(:,k,l))
                   call LVT_readvar_restart(ftn,&
                        stats%ubrmse%count_value_total(:,k,l))
                enddo
             enddo
          endif
          
          if(LVT_metrics%ubRMSE%computeSC.eq.1) then
             do k=1,model%selectNlevs
                do l=1,LVT_rc%nasc          
                   call LVT_readvar_restart(ftn,&
                        stats%ubrmse%value_asc(:,k,l))
                   call LVT_readvar_restart(ftn,&
                        stats%ubrmse%count_value_asc(:,k,l))
                enddo
             enddo
          endif
          if(LVT_metrics%ubRMSE%computeADC.eq.1) then
             do k=1,model%selectNlevs
                do l=1,LVT_rc%nadc           
                   call LVT_readvar_restart(ftn,&
                        stats%ubrmse%value_adc(:,k,l))
                   call LVT_readvar_restart(ftn,&
                        stats%ubrmse%count_value_adc(:,k,l))
                enddo
             enddo
          endif
       endif

       model => model%next
       obs   => obs%next
       stats => stats%next
    enddo
  end subroutine LVT_readrestart_ubRMSE

end module LVT_ubRMSEMod
