!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_AnomalyCorrMod
! \label(LVT_AnomalyCorrMod)
!
! !INTERFACE:
module LVT_AnomalyCorrMod
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
!   This module handles the Anomaly correlation (pearson correlation
!   coefficient)  computations by comparing 
!   the LIS output to the specified observations. The anomaly 
!   correlation coefficients are computed by subtracting the monthly
!   mean climatology of each dataset from the corresponding raw data. 
!
!  Time series, seasonal cycles, average diurnal cycle options are 
!  not supported for Anomaly Correlation computation
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
  public :: LVT_initAnomalyCorr
  public :: LVT_diagnoseAnomalyCorr
  public :: LVT_computeAnomalyCorr 
  public :: LVT_writeMetric_AnomalyCorr
  public :: LVT_resetMetric_AnomalyCorr
  public :: LVT_writerestart_AnomalyCorr
  public :: LVT_readrestart_AnomalyCorr

!  integer :: minObsCountForClimo = 8
  integer :: minObsCountForClimo = 1
contains
  
!BOP
! 
! !ROUTINE: LVT_initAnomalyCorr
! \label{LVT_initAnomalyCorr}
!
! !INTERFACE: 
  subroutine LVT_initAnomalyCorr(model,obs,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes the data structures required for 
!  anomaly correlation computations. 
!  
!  The arguments are:  
!  \begin{description}
!   \item[model]  model variable object
!   \item[obs]    observation object
!   \item[stats]  object to hold the updated statistics
!   \item[metric] object representing the anomaly correlation metric
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
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!EOP
    if(metric%selectOpt.ge.1) then 
       allocate(stats%acorr%model_value_climo(LVT_rc%ngrid,model%selectNlevs, &
            12, LVT_rc%strat_nlevels))
       allocate(stats%acorr%obs_value_climo(LVT_rc%ngrid,model%selectNlevs, &
            12,LVT_rc%strat_nlevels))
       allocate(stats%acorr%count_model_value_climo(LVT_rc%ngrid,model%selectNlevs,&
            12,LVT_rc%strat_nlevels))
       allocate(stats%acorr%count_obs_value_climo(LVT_rc%ngrid,model%selectNlevs,&
            12,LVT_rc%strat_nlevels))
       
       allocate(stats%acorr%sxy_a(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
       allocate(stats%acorr%sxx_a(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
       allocate(stats%acorr%syy_a(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
       allocate(stats%acorr%sx_a(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
       allocate(stats%acorr%sy_a(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
       allocate(stats%acorr%count_value(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))
       allocate(stats%acorr%rval_a(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))          
       
       stats%acorr%model_value_climo = 0 
       stats%acorr%obs_value_climo   = 0 
       stats%acorr%count_model_value_climo = 0 
       stats%acorr%count_obs_value_climo   = 0 
       
       stats%acorr%sxx_a = 0 
       stats%acorr%syy_a = 0 
       stats%acorr%sxy_a = 0 
       stats%acorr%sx_a = 0 
       stats%acorr%sy_a = 0 
       stats%acorr%count_value = 0 
       stats%acorr%rval_a = 0 
       
       allocate(stats%acorr%rval_a_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%acorr%rval_a_ci = LVT_rc%udef

       if(metric%timeOpt.eq.1) then
          
          allocate(stats%acorr%sxy_ts_a(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%acorr%sxx_ts_a(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%acorr%syy_ts_a(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%acorr%sx_ts_a(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%acorr%sy_ts_a(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%acorr%count_value_ts(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%strat_nlevels))
          allocate(stats%acorr%rval_ts_a(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))          
          
          stats%acorr%sxx_ts_a = 0 
          stats%acorr%syy_ts_a = 0 
          stats%acorr%sxy_ts_a = 0 
          stats%acorr%sx_ts_a = 0 
          stats%acorr%sy_ts_a = 0 
          stats%acorr%count_value_ts = 0 
          stats%acorr%rval_ts_a = 0 

          allocate(stats%acorr%tavg_value_ts(LVT_rc%ngrid, &
               model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%acorr%tavg_count_value_ts(LVT_rc%ngrid, &
               model%selectNlevs,LVT_rc%strat_nlevels))
          stats%acorr%tavg_value_ts = 0.0
          stats%acorr%tavg_count_value_ts=0 

          if(metric%computeSC.eq.1) then 
             allocate(stats%acorr%value_asc(LVT_rc%ngrid, model%selectNlevs, LVT_rc%nasc))
             stats%acorr%value_asc = 0.0
             allocate(stats%acorr%count_value_asc(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%nasc))
             stats%acorr%count_value_asc = 0
          endif

          if(metric%computeADC.eq.1) then 
             allocate(stats%acorr%value_adc(LVT_rc%ngrid, model%selectNlevs, LVT_rc%nadc))
             stats%acorr%value_adc = 0.0
             allocate(stats%acorr%count_value_adc(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%nadc))
             stats%acorr%count_value_adc = 0
          endif

       endif
    endif
    
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 2
    metric%obsData = .false. 
    metric%stdevFlag = .false. 

  end subroutine LVT_initAnomalyCorr
  
!BOP
! 
! !ROUTINE: LVT_diagnoseAnomalyCorr
! \label{LVT_diagnoseAnomalyCorr}
!
! !INTERFACE: 
  subroutine LVT_diagnoseAnomalyCorr(pass)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the anomaly correlation
!   calculation for desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleClimatology](\ref{diagnoseSingleClimatology})
!     updates the monthly climatology computation for a single variable 
!    \item[diagnoseSingleAnomalyCorr](\ref{diagnoseSingleAnomalyCorr})
!     updates the anomaly correlation computation for a single variable 
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
       if(LVT_metrics%acorr%selectOpt.ge.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleClimatology(&
                  obs, model, stats)

             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo
       endif
    elseif(pass.eq.2) then 
       if(LVT_metrics%acorr%selectOpt.ge.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call diagnoseSingleAnomalyCorr(&
                  obs, model, stats, &
                  LVT_metrics%acorr)

             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseAnomalyCorr

!BOP
! 
! !ROUTINE: diagnoseSingleClimatology
! \label{diagnoseSingleClimatology}
!
! !INTERFACE: 
  subroutine diagnoseSingleClimatology(obs, model, stats)
! 
! !USES: 
    use LVT_coreMod,    only : LVT_rc, LVT_domain
    use LVT_histDataMod 
    use LVT_statsDataMod
    use LVT_logMod,     only : LVT_logunit, LVT_endrun

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the monthly climatology computation of the 
!   specified variable. 
!
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
!BOP
!
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry) :: stats
! 
!EOP
    integer    :: t,k

    if(stats%selectOpt.ge.1.and.obs%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             if(model%count(t,k).ne.0.and.obs%count(t,k).ne.0) then 
                stats%acorr%model_value_climo(t,k,LVT_rc%mo,1) = &
                     stats%acorr%model_value_climo(t,k,LVT_rc%mo,1) + &
                     model%value(t,k)
                stats%acorr%count_model_value_climo(t,k,LVT_rc%mo,1) = &
                     stats%acorr%count_model_value_climo(t,k,LVT_rc%mo,1) + 1
                
             endif
             if(obs%count(t,k).ne.0) then  
                stats%acorr%obs_value_climo(t,k,LVT_rc%mo,1) = &
                     stats%acorr%obs_value_climo(t,k,LVT_rc%mo,1) + &
                     obs%value(t,k)
                stats%acorr%count_obs_value_climo(t,k,LVT_rc%mo,1) = &
                     stats%acorr%count_obs_value_climo(t,k,LVT_rc%mo,1) + 1
                
             endif
             if(LVT_rc%strat_nlevels.gt.1) then 
                if(LVT_stats%strat_var(t,k).gt.&
                     LVT_rc%strat_var_threshold) then 
                   if(model%count(t,k).ne.0.and.obs%count(t,k).ne.0) then 
                      stats%acorr%model_value_climo(t,k,LVT_rc%mo,2) = &
                           stats%acorr%model_value_climo(t,k,LVT_rc%mo,2) + &
                           model%value(t,k)
                      stats%acorr%count_model_value_climo(t,k,LVT_rc%mo,2) = &
                           stats%acorr%count_model_value_climo(t,k,LVT_rc%mo,2) + 1
                   endif
                   if(obs%count(t,k).ne.0) then  
                      stats%acorr%obs_value_climo(t,k,LVT_rc%mo,2) = &
                           stats%acorr%obs_value_climo(t,k,LVT_rc%mo,2) + &
                           obs%value(t,k)
                      stats%acorr%count_obs_value_climo(t,k,LVT_rc%mo,2) = &
                           stats%acorr%count_obs_value_climo(t,k,LVT_rc%mo,2) + 1
                   endif
                elseif(LVT_stats%strat_var(t,k).le.&
                     LVT_rc%strat_var_threshold) then 
                   if(model%count(t,k).ne.0.and.obs%count(t,k).ne.0) then 
                      stats%acorr%model_value_climo(t,k,LVT_rc%mo,3) = &
                           stats%acorr%model_value_climo(t,k,LVT_rc%mo,3) + &
                           model%value(t,k)
                      stats%acorr%count_model_value_climo(t,k,LVT_rc%mo,3) = &
                           stats%acorr%count_model_value_climo(t,k,LVT_rc%mo,3) + 1
                   endif
                   if(obs%count(t,k).ne.0) then  
                      stats%acorr%obs_value_climo(t,k,LVT_rc%mo,3) = &
                           stats%acorr%obs_value_climo(t,k,LVT_rc%mo,3) + &
                           obs%value(t,k)
                      stats%acorr%count_obs_value_climo(t,k,LVT_rc%mo,3) = &
                           stats%acorr%count_obs_value_climo(t,k,LVT_rc%mo,3) + 1
                   endif
                endif
             endif
          enddo
       enddo
    endif

  end subroutine diagnoseSingleClimatology

!BOP
! 
! !ROUTINE: diagnoseSingleAnomalyCorr
! \label{diagnoseSingleAnomalyCorr}
!
! !INTERFACE: 
  subroutine diagnoseSingleAnomalyCorr(obs, model, stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the anomaly correlation computation (updates 
!  the running sum calculations of different terms). 
!
!  The arguments are: 
!  \begin{description}
!   \item[obs] observation object
!   \item[model] model variable object
!   \item[stats] object to hold the updated statistics
!   \item[metric] object representing the anomaly correlation metric
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

    if(stats%selectOpt.ge.1.and.obs%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             if(obs%count(t,k).ne.0.and. &
                  model%count(t,k).ne.0) then      
                if(metric%selectOpt.ge.1.and.&
                     stats%acorr%model_value_climo(t,k,LVT_rc%mo,1).ne.LVT_rc%udef.and.&
                     stats%acorr%obs_value_climo(t,k,LVT_rc%mo,1).ne.LVT_rc%udef) then 
                   stats%acorr%sxy_a(t,k,1) = stats%acorr%sxy_a(t,k,1) + &
                        (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,1))*&
                        (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,1))
                   stats%acorr%sx_a(t,k,1) = stats%acorr%sx_a(t,k,1) +&
                        (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,1))
                   stats%acorr%sy_a(t,k,1) = stats%acorr%sy_a(t,k,1) +&
                        (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,1))
                   stats%acorr%sxx_a(t,k,1) = stats%acorr%sxx_a(t,k,1) +&
                        (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,1))*&
                        (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,1))
                   stats%acorr%syy_a(t,k,1) = stats%acorr%syy_a(t,k,1) +&
                        (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,1))*&
                        (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,1))
                   stats%acorr%count_value(t,k,1) = stats%acorr%count_value(t,k,1) + 1
                   if(LVT_rc%strat_nlevels.gt.1) then 
                      if(LVT_stats%strat_var(t,k).gt.&
                           LVT_rc%strat_var_threshold) then 
                         stats%acorr%sxy_a(t,k,2) = stats%acorr%sxy_a(t,k,2) + &
                              (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,2))*&
                              (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,2))
                         stats%acorr%sx_a(t,k,2) = stats%acorr%sx_a(t,k,2) +&
                              (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,2))
                         stats%acorr%sy_a(t,k,2) = stats%acorr%sy_a(t,k,2) +&
                              (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,2))
                         stats%acorr%sxx_a(t,k,2) = stats%acorr%sxx_a(t,k,2) +&
                              (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,2))*&
                              (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,2))
                         stats%acorr%syy_a(t,k,2) = stats%acorr%syy_a(t,k,2) +&
                              (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,2))*&
                              (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,2))
                         stats%acorr%count_value(t,k,2) = stats%acorr%count_value(t,k,2) + 1
                      elseif(LVT_stats%strat_var(t,k).le.&
                           LVT_rc%strat_var_threshold) then 
                         stats%acorr%sxy_a(t,k,3) = stats%acorr%sxy_a(t,k,3) + &
                              (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,3))*&
                              (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,3))
                         stats%acorr%sx_a(t,k,3) = stats%acorr%sx_a(t,k,3) +&
                              (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,3))
                         stats%acorr%sy_a(t,k,3) = stats%acorr%sy_a(t,k,3) +&
                              (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,3))
                         stats%acorr%sxx_a(t,k,3) = stats%acorr%sxx_a(t,k,3) +&
                              (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,3))*&
                              (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,3))
                         stats%acorr%syy_a(t,k,3) = stats%acorr%syy_a(t,k,3) +&
                              (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,3))*&
                              (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,3))
                         stats%acorr%count_value(t,k,3) = stats%acorr%count_value(t,k,3) + 1
                         
                      endif

                      if(metric%timeOpt.eq.1) then 
                         
                         stats%acorr%sxy_ts_a(t,k,1) = stats%acorr%sxy_ts_a(t,k,1) + &
                              (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,1))*&
                              (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,1))
                         stats%acorr%sx_ts_a(t,k,1) = stats%acorr%sx_ts_a(t,k,1) +&
                              (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,1))
                         stats%acorr%sy_ts_a(t,k,1) = stats%acorr%sy_ts_a(t,k,1) +&
                              (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,1))
                         stats%acorr%sxx_ts_a(t,k,1) = stats%acorr%sxx_ts_a(t,k,1) +&
                              (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,1))*&
                              (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,1))
                         stats%acorr%syy_ts_a(t,k,1) = stats%acorr%syy_ts_a(t,k,1) +&
                              (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,1))*&
                              (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,1))
                         stats%acorr%count_value_ts(t,k,1) = stats%acorr%count_value_ts(t,k,1) + 1
                         if(LVT_rc%strat_nlevels.gt.1) then 
                            if(LVT_stats%strat_var(t,k).gt.&
                                 LVT_rc%strat_var_threshold) then 
                               stats%acorr%sxy_ts_a(t,k,2) = stats%acorr%sxy_ts_a(t,k,2) + &
                                    (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,2))*&
                                    (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,2))
                               stats%acorr%sx_ts_a(t,k,2) = stats%acorr%sx_ts_a(t,k,2) +&
                                    (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,2))
                               stats%acorr%sy_ts_a(t,k,2) = stats%acorr%sy_ts_a(t,k,2) +&
                                    (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,2))
                               stats%acorr%sxx_ts_a(t,k,2) = stats%acorr%sxx_ts_a(t,k,2) +&
                                    (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,2))*&
                                    (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,2))
                               stats%acorr%syy_ts_a(t,k,2) = stats%acorr%syy_ts_a(t,k,2) +&
                                    (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,2))*&
                                    (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,2))
                               stats%acorr%count_value_ts(t,k,2) = stats%acorr%count_value_ts(t,k,2) + 1
                            elseif(LVT_stats%strat_var(t,k).le.&
                                 LVT_rc%strat_var_threshold) then 
                               stats%acorr%sxy_ts_a(t,k,3) = stats%acorr%sxy_ts_a(t,k,3) + &
                                    (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,3))*&
                                    (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,3))
                               stats%acorr%sx_ts_a(t,k,3) = stats%acorr%sx_ts_a(t,k,3) +&
                                    (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,3))
                               stats%acorr%sy_ts_a(t,k,3) = stats%acorr%sy_ts_a(t,k,3) +&
                                    (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,3))
                               stats%acorr%sxx_ts_a(t,k,3) = stats%acorr%sxx_ts_a(t,k,3) +&
                                    (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,3))*&
                                    (model%value(t,k)-stats%acorr%model_value_climo(t,k,LVT_rc%mo,3))
                               stats%acorr%syy_ts_a(t,k,3) = stats%acorr%syy_ts_a(t,k,3) +&
                                    (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,3))*&
                                    (obs%value(t,k)-stats%acorr%obs_value_climo(t,k,LVT_rc%mo,3))
                               stats%acorr%count_value_ts(t,k,3) = stats%acorr%count_value_ts(t,k,3) + 1
                               
                            endif
                         endif
                      endif

                   endif
                endif
             endif
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleAnomalyCorr


!BOP
! 
! !ROUTINE: LVT_computeAnomalyCorr
! \label{LVT_computeAnomalyCorr}
!
! !INTERFACE: 
  subroutine LVT_computeAnomalyCorr(pass,alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute anomaly correlation values 
!   for the desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleClimatology](\ref{computeSingleClimatology})
!     updates the monthly climatology computation for a single variable 
!    \item[computeSingleAnomalyCorr](\ref{computeSingleAnomalyCorr})
!     updates the anomaly correlation computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[pass]
!     current pass number over the data points
!    \item[alarm]
!     flag indicating if the temporal output interval has been reached
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS:
    implicit none

    integer     :: pass
    logical     :: alarm
    integer     :: i
! 
!EOP
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(LVT_metrics%acorr%selectOpt.ge.1) then 
       if(pass.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call computeSingleClimatology(&
                  obs, model, stats)
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo

       elseif(pass.eq.2) then 
          
          if(LVT_metrics%acorr%selectOpt.eq.1.or.&
               LVT_metrics%acorr%timeOpt.eq.1) then 
             if(alarm) then 
                if(LVT_metrics%acorr%timeOpt.eq.1.and.&
                     LVT_metrics%acorr%extractTS.eq.1) then 
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%acorr%ftn_ts_loc(i),200,advance='no') &
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
                call computeSingleAnomalyCorr(alarm,&
                     obs, model, stats, &
                     LVT_metrics%acorr)
                model => model%next
                obs => obs%next
                stats => stats%next
             enddo

             if(alarm) then 
                if(LVT_metrics%acorr%timeOpt.eq.1.and.&
                     LVT_metrics%acorr%extractTS.eq.1) then 
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%acorr%ftn_ts_loc(i),fmt='(a1)') ''
                   enddo
                endif
             endif
          endif
       endif
    endif

  end subroutine LVT_ComputeAnomalyCorr

!BOP
! 
! !ROUTINE: computeSingleClimatology
! \label{computeSingleClimatology}
!
! !INTERFACE: 
  subroutine computeSingleClimatology(obs, model, stats)
! 
! !USES: 
    use LVT_coreMod,    only : LVT_rc, LVT_domain
    use LVT_histDataMod 
    use LVT_statsDataMod
    use LVT_logMod,     only : LVT_logunit

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the monthly climatolgy for a specified variable
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
!BOP
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry) :: stats
!EOP
    integer    :: t,k,m,l

    if(LVT_rc%endtime.eq.1) then 
       write(LVT_logunit,*) 'Computing monthly climatology '
       if(stats%selectOpt.ge.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do m=1,12
                   do l=1, LVT_rc%strat_nlevels
                      if(stats%acorr%count_model_value_climo(t,k,m,l) & 
                           .gt. minObsCountForClimo ) then 
                         stats%acorr%model_value_climo(t,k,m,l) = &
                              stats%acorr%model_value_climo(t,k,m,l) /&
                              stats%acorr%count_model_value_climo(t,k,m,l)
                      else
                         stats%acorr%model_value_climo(t,k,m,l) = LVT_rc%udef
                      endif
                      if(stats%acorr%count_obs_value_climo(t,k,m,l) &  
                         .gt.  minObsCountForClimo )then 
                         stats%acorr%obs_value_climo(t,k,m,l) = &
                              stats%acorr%obs_value_climo(t,k,m,l)/&
                              stats%acorr%count_obs_value_climo(t,k,m,l)
                      else
                         stats%acorr%obs_value_climo(t,k,m,l) = LVT_rc%udef
                      endif
                   enddo
                enddo
             enddo
          enddo
       endif

    endif
    
  end subroutine computeSingleClimatology

!BOP
! 
! !ROUTINE: computeSingleAnomalyCorr
! \label{computeSingleAnomalyCorr}
!
! !INTERFACE: 
  subroutine computeSingleAnomalyCorr(alarm,obs, model,stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the AnomalyCorr values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[alarm]
!     flag indicating if the temporal output interval has been reached
!   \item[obs] observation object
!   \item[model] model variable object
!   \item[stats] object to hold the updated statistics
!   \item[metric] object representing the anomaly correlation metric
!  \end{description}
!
!   The methods invoked are: 
!   \begin{description}
!    \item[LVT\_computeCI](\ref{LVT_computeCI})
!    computes the confidence intervals on the anomaly correlation 
!    coefficients
!    \item[LVT\_writeDataBasedStrat](\ref{LVT_writeDataBasedStrat})
!    stratifies the anomaly correlation coefficient values based on 
!    the specified external data
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS: 
    logical                 :: alarm
    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

! 
!EOP
    integer  :: t,l,k,tind
    real     :: numer
    real     :: denom_sq
    real     :: denom


    if(metric%timeOpt.eq.1) then 
       
       if(alarm) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%acorr%count_value_ts(t,k,l).ne.0.and.&
                        stats%acorr%count_value_ts(t,k,l) &
                        .gt.LVT_rc%obsCountThreshold) then 
                      
                      numer = (float(stats%acorr%count_value_ts(t,k,l))* &
                           stats%acorr%sxy_ts_a(t,k,l) - &
                           stats%acorr%sx_ts_a(t,k,l)*stats%acorr%sy_ts_a(t,k,l))
                      denom_sq =  (float(stats%acorr%count_value_ts(t,k,l))* &
                           stats%acorr%sxx_ts_a(t,k,l)-&
                           stats%acorr%sx_ts_a(t,k,l)**2)* &
                           (float(stats%acorr%count_value_ts(t,k,l))*&
                           stats%acorr%syy_ts_a(t,k,l)-&
                           stats%acorr%sy_ts_a(t,k,l)**2)
                      if(denom_sq.gt.0) then 
                         denom = sqrt(denom_sq)
                      else
                         denom = 0.0
                      endif

                      if(denom.ne.0) then 
                         stats%acorr%rval_ts_a(t,k,l) = (numer/denom)
                         
                         stats%acorr%tavg_value_ts(t,k,l) = & 
                              stats%acorr%tavg_value_ts(t,k,l) + & 
                              stats%acorr%rval_ts_a(t,k,l)
                         
                         stats%acorr%tavg_count_value_ts(t,k,l) = & 
                              stats%acorr%tavg_count_value_ts(t,k,l) + 1
                      else
                         stats%acorr%rval_ts_a(t,k,l) = LVT_rc%udef
                      endif
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   if(stats%acorr%rval_ts_a(t,k,l).ne.LVT_rc%udef) then 
                      call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                           tind)
                      stats%acorr%value_asc(t,k,tind) = & 
                           stats%acorr%value_asc(t,k,tind) + & 
                           stats%acorr%rval_ts_a(t,k,l)
                      stats%acorr%count_value_asc(t,k,tind) = & 
                           stats%acorr%count_value_asc(t,k,tind) + 1
                   endif
                endif
                
                if(metric%computeADC.eq.1) then 
                   if(stats%acorr%rval_ts_a(t,k,l).ne.LVT_rc%udef) then 
                      call LVT_getADCTimeIndex(tind)
                      stats%acorr%value_adc(t,k,tind) = & 
                           stats%acorr%value_adc(t,k,tind) + & 
                           stats%acorr%rval_ts_a(t,k,l)
                      stats%acorr%count_value_adc(t,k,tind) = & 
                           stats%acorr%count_value_adc(t,k,tind) + 1
                   endif
                endif
             enddo
          enddo
          
          if(metric%extractTS.eq.1) then 
             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%acorr%tavg_count_value_ts(t,k,l).gt.0) then 
                         stats%acorr%tavg_value_ts(t,k,l) = &
                              stats%acorr%tavg_value_ts(t,k,l)/&
                              stats%acorr%tavg_count_value_ts(t,k,l)
                      else
                         stats%acorr%tavg_value_ts(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                enddo
             enddo

             call LVT_writeTSinfo(metric%ftn_ts_loc,&
                  model,&
                  LVT_rc%ngrid,&
                  stats%acorr%tavg_value_ts,&
                  stats%acorr%tavg_count_value_ts)                
          endif
       endif
    endif

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.ge.1) then 
       if(stats%selectOpt.ge.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%acorr%count_value(t,k,l).ne.0.and.&
                        stats%acorr%count_value(t,k,l) &
                        .gt.LVT_rc%obsCountThreshold) then 
                      
                      numer = (float(stats%acorr%count_value(t,k,l))* &
                           stats%acorr%sxy_a(t,k,l) - &
                           stats%acorr%sx_a(t,k,l)*stats%acorr%sy_a(t,k,l))
                      denom_sq =  (float(stats%acorr%count_value(t,k,l))* &
                           stats%acorr%sxx_a(t,k,l)-&
                           stats%acorr%sx_a(t,k,l)**2)* &
                           (float(stats%acorr%count_value(t,k,l))*&
                           stats%acorr%syy_a(t,k,l)-&
                           stats%acorr%sy_a(t,k,l)**2)
                      if(denom_sq.gt.0) then 
                         denom = sqrt(denom_sq)
                      else
                         denom = 0.0
                      endif

                      if(denom.ne.0) then 
                         stats%acorr%rval_a(t,k,l) = (numer/denom)
                      else
                         stats%acorr%rval_a(t,k,l) = LVT_rc%udef
                         stats%acorr%count_value(t,k,l) = LVT_rc%udef
                      endif
                   else
                      stats%acorr%rval_a(t,k,l) = LVT_rc%udef
                      stats%acorr%count_value(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1,LVT_rc%nasc
                      if(stats%acorr%count_value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then 
                         stats%acorr%value_asc(t,k,l) = &
                              stats%acorr%value_asc(t,k,l)/&
                              stats%acorr%count_value_asc(t,k,l) 
                      else
                         stats%acorr%value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
                if(metric%computeADC.eq.1) then 
                   do l=1,LVT_rc%nadc
                      if(stats%acorr%count_value_adc(t,k,l).gt.&
                           LVT_rc%ADCCountThreshold) then 
                         stats%acorr%value_adc(t,k,l) = &
                              stats%acorr%value_adc(t,k,l)/&
                              stats%acorr%count_value_adc(t,k,l) 
                      else
                         stats%acorr%value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
                
             enddo
          enddo
          
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%acorr%rval_a(:,k,l),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%acorr%rval_a_ci(k,l))
             enddo
          enddo
       
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid,stats%acorr%rval_a)      

          if(metric%computeSC.eq.1) then 
             call LVT_writeSeasonalCycleInfo(model,obs,stats, metric, &
                  LVT_rc%ngrid, stats%acorr%value_asc, stats%acorr%count_value_asc)
          endif
          if(metric%computeADC.eq.1) then 

             call LVT_writeAvgDiurnalCycleInfo(model, obs, stats, metric, & 
                  LVT_rc%ngrid, stats%acorr%value_adc, stats%acorr%count_value_adc)
          endif
       endif
    endif

  end subroutine computeSingleAnomalyCorr

!BOP
! 
! !ROUTINE: LVT_writeMetric_AnomalyCorr
! \label{LVT_writeMetric_AnomalyCorr}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_AnomalyCorr(pass,final,vlevels,stats,obs)
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
!   This subroutine writes the computed anomaly correlation values to an 
!   external file. 
! 
!  \begin{description}
!   \item[pass]  current pass number over the data points
!   \item[final] boolean value indicating if computation is for the final
!                timestep (=1) or not (=0)
!   \item[vlevels] number of vertical levels of the variable being written 
!   \item[stats] object to hold the updated statistics
!  \end{description}
!   The methods invoked are: 
!   \begin{description}
!    \item[LVT\_writevar\_gridded](\ref{LVT_writevar_gridded})
!    writes the variable to a gridded output file. 
!    \item[LVT\_writeSummaryStats](\ref{LVT_writeSummaryStats})
!    writes the domain averaged summary statistics to a text file
!   \end{description}
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
    integer                 :: l
    integer                 :: k
    integer                 :: tind

    if(pass.eq.LVT_metrics%acorr%npass) then 
       if(LVT_metrics%acorr%selectOpt.ge.1) then
          if(stats%selectOpt.ge.1) then 
             
             if(final.ne.1) then
                
                call LVT_writevar_gridded(LVT_metrics%acorr%ftn_ts, &
                     stats%acorr%tavg_value_ts(:,k,l),stats%vid_ts(LVT_Acorrid,1))
                call LVT_writevar_gridded(LVT_metrics%acorr%ftn_ts, &
                     real(stats%acorr%tavg_count_value_ts(:,k,l)),&
                     stats%vid_count_ts(LVT_Acorrid,1))                       
             else
                do k=1,vlevels
                   do l=1,LVT_rc%strat_nlevels
                      
                      call LVT_writevar_gridded(LVT_metrics%acorr%ftn_total, &
                           stats%acorr%rval_a(:,k,l),stats%vid_total(LVT_Acorrid,1),k)
                      call LVT_writevar_gridded(LVT_metrics%acorr%ftn_total, &
                           real(stats%acorr%count_value(:,k,l)),&
                           stats%vid_count_total(LVT_Acorrid,1),k)                
                   enddo

                   if(LVT_metrics%acorr%computeSC.eq.1) then 
                      do tind = 1, LVT_rc%nasc
                         call LVT_writevar_gridded(&
                              LVT_metrics%acorr%ftn_total,&
                              stats%acorr%value_asc(:,k,tind), &
                              stats%vid_sc_total(tind,LVT_Acorrid,1),k)
                      enddo
                   endif
                   
                   if(LVT_metrics%acorr%computeADC.eq.1) then 
                      do tind = 1,LVT_rc%nadc
                         call LVT_writevar_gridded(LVT_metrics%acorr%ftn_total,&
                              stats%acorr%value_adc(:,k,tind),&
                              stats%vid_adc_total(tind,LVT_Acorrid,1),k)
                      enddo
                   endif
                   call LVT_writeSummaryStats(&
                        LVT_metrics%acorr%ftn_summ,&
                        LVT_metrics%acorr%short_name,&
                        LVT_rc%ngrid,&
                        stats%acorr%rval_a(:,k,:), &
                        stats%acorr%count_value(:,k,:),stats%standard_name,&
                        stats%acorr%rval_a_ci(k,:))
                enddo
             endif
          endif
       endif
    endif

  end subroutine LVT_writeMetric_AnomalyCorr

!BOP
! 
! !ROUTINE: LVT_resetMetric_AnomalyCorr
! \label{LVT_resetMetric_AnomalyCorr}
!
! !INTERFACE:
  subroutine LVT_resetMetric_AnomalyCorr(alarm)
!
! !INPUT PARAMETERS: 
    logical           :: alarm
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This subroutine resets the relevant variables during the anomaly 
!  correlation computation between each temporal output times
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
             if(LVT_metrics%acorr%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   if(alarm) then 
                      stats%acorr%sxy_ts_a(:,k,l) = 0.0
                      stats%acorr%sxx_ts_a(:,k,l) = 0.0
                      stats%acorr%syy_ts_a(:,k,l) = 0.0
                      stats%acorr%sx_ts_a(:,k,l) = 0.0
                      stats%acorr%sy_ts_a(:,k,l) = 0.0
                      stats%acorr%rval_ts_a(:,k,l) = 0.0
                      stats%acorr%count_value_ts(:,k,l)=0 
                      
                      stats%acorr%tavg_value_ts(:,k,l) = 0.0
                      stats%acorr%tavg_count_value_ts(:,k,l)=0 
                   endif

                enddo
             endif
             
          enddo
          
       endif
       model => model%next
       stats => stats%next

    enddo

  end subroutine LVT_resetMetric_AnomalyCorr

!BOP
! 
! !ROUTINE: LVT_writerestart_AnomalyCorr
! 
! !INTERFACE:
  subroutine LVT_writerestart_AnomalyCorr(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for AnomalyCorr metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    integer              :: k,m,l
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(LVT_metrics%acorr%selectOpt.ge.1) then 
       
       call LVT_getDataStream1Ptr(model)
       call LVT_getDataStream2Ptr(obs)
       call LVT_getstatsEntryPtr(stats)
       
       do while(associated(model))
          
          if(stats%selectOpt.ge.1.and.obs%selectNlevs.ge.1) then 
             do k=1,model%selectNlevs
                do m=1,12
                   do l=1,LVT_rc%strat_nlevels     
                      call LVT_writevar_restart(ftn,&
                           stats%acorr%model_value_climo(:,k,m,l))
                      call LVT_writevar_restart(ftn,&
                           stats%acorr%obs_value_climo(:,k,m,l))
                      call LVT_writevar_restart(ftn,&
                           stats%acorr%count_model_value_climo(:,k,m,l))
                      call LVT_writevar_restart(ftn,&
                           stats%acorr%count_obs_value_climo(:,k,m,l))
                   enddo
                enddo
             enddo
             if(pass.eq.2) then
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels      
                      call LVT_writevar_restart(ftn,stats%acorr%sxy_a(:,k,l))
                      call LVT_writevar_restart(ftn,stats%acorr%sxx_a(:,k,l))
                      call LVT_writevar_restart(ftn,stats%acorr%syy_a(:,k,l))
                      call LVT_writevar_restart(ftn,stats%acorr%sx_a(:,k,l))
                      call LVT_writevar_restart(ftn,stats%acorr%sy_a(:,k,l))
                      call LVT_writevar_restart(ftn,stats%acorr%count_value(:,k,l))
                      call LVT_writevar_restart(ftn,stats%acorr%rval_a(:,k,l))
                   enddo
                enddo
                if(LVT_metrics%acorr%computeSC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nasc
                         call LVT_writevar_restart(ftn,&
                              stats%acorr%value_asc(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%acorr%count_value_asc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_metrics%acorr%computeADC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nadc
                         call LVT_writevar_restart(ftn,&
                              stats%acorr%value_adc(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%acorr%count_value_adc(:,k,l))
                      enddo
                   enddo
                endif
                
             endif
          endif

          model => model%next
          obs => obs%next
          stats => stats%next
       enddo       
    endif
  end subroutine LVT_writerestart_AnomalyCorr

!BOP
! 
! !ROUTINE: LVT_readrestart_AnomalyCorr
! 
! !INTERFACE:
  subroutine LVT_readrestart_AnomalyCorr(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for AnomalyCorr metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    integer              :: k,m,l
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(LVT_metrics%acorr%selectOpt.ge.1) then 
       
       call LVT_getDataStream1Ptr(model)
       call LVT_getDataStream2Ptr(obs)
       call LVT_getstatsEntryPtr(stats)
       
       do while(associated(model)) 
          if(stats%selectOpt.ge.1.and.obs%selectNlevs.ge.1) then 
             do k=1,model%selectNlevs
                do m=1,12
                   do l=1,LVT_rc%strat_nlevels     
                      call LVT_readvar_restart(ftn,&
                           stats%acorr%model_value_climo(:,k,m,l))
                      call LVT_readvar_restart(ftn,&
                           stats%acorr%obs_value_climo(:,k,m,l))
                      call LVT_readvar_restart(ftn,&
                           stats%acorr%count_model_value_climo(:,k,m,l))
                      call LVT_readvar_restart(ftn,&
                           stats%acorr%count_obs_value_climo(:,k,m,l))
                   enddo
                enddo
             enddo
             if(LVT_rc%curr_pass.eq.2) then
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels      
                      call LVT_readvar_restart(ftn,stats%acorr%sxy_a(:,k,l))
                      call LVT_readvar_restart(ftn,stats%acorr%sxx_a(:,k,l))
                      call LVT_readvar_restart(ftn,stats%acorr%syy_a(:,k,l))
                      call LVT_readvar_restart(ftn,stats%acorr%sx_a(:,k,l))
                      call LVT_readvar_restart(ftn,stats%acorr%sy_a(:,k,l))
                      call LVT_readvar_restart(ftn,stats%acorr%count_value(:,k,l))
                      call LVT_readvar_restart(ftn,stats%acorr%rval_a(:,k,l))
                   enddo
                enddo
                if(LVT_metrics%acorr%computeSC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nasc
                         call LVT_readvar_restart(ftn,&
                              stats%acorr%value_asc(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%acorr%count_value_asc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_metrics%acorr%computeADC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nadc
                         call LVT_readvar_restart(ftn,&
                              stats%acorr%value_adc(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%acorr%count_value_adc(:,k,l))
                      enddo
                   enddo
                endif
             endif
          endif
          model => model%next
          obs => obs%next
          stats => stats%next
       enddo
    endif
  end subroutine LVT_readrestart_AnomalyCorr

end module LVT_AnomalyCorrMod
