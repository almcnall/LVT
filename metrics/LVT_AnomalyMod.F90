!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_AnomalyMod
! \label(LVT_AnomalyMod)
!
! !INTERFACE:
module LVT_AnomalyMod
! 
! !USES:   
  use LVT_coreMod
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
!  This module handles the computations required to compute mean values
!  of desired variables from the LIS output
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initAnomaly
  public :: LVT_diagnoseAnomaly
  public :: LVT_computeAnomaly
  public :: LVT_writeMetric_Anomaly
  public :: LVT_resetMetric_Anomaly
  public :: LVT_writerestart_Anomaly
  public :: LVT_readrestart_Anomaly
!EOP
  
  private

contains
  subroutine LVT_initAnomaly(model, obs, stats,metric)
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

       allocate(stats%anomaly%model_value_climo(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%anomalyTlength, LVT_rc%strat_nlevels))
       allocate(stats%anomaly%obs_value_climo(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%anomalyTlength, LVT_rc%strat_nlevels))
       
       allocate(stats%anomaly%count_model_value_climo(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%anomalyTlength, LVT_rc%strat_nlevels))
       allocate(stats%anomaly%count_obs_value_climo(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%anomalyTlength, LVT_rc%strat_nlevels))

       stats%anomaly%model_value_climo = 0 
       stats%anomaly%obs_value_climo = 0 
       stats%anomaly%count_model_value_climo = 0 
       stats%anomaly%count_obs_value_climo = 0 

       allocate(stats%anomaly%model_value_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%anomaly%model_value_total = 0.0
       allocate(stats%anomaly%count_model_value_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%anomaly%count_model_value_total = 0
       allocate(stats%anomaly%model_value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%anomaly%model_value_ci = LVT_rc%udef
       
       if(metric%computeSC.eq.1) then 
          allocate(stats%anomaly%model_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%anomaly%model_value_asc = 0.0
          allocate(stats%anomaly%count_model_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%anomaly%count_model_value_asc = 0
       endif
       if(LVT_rc%obssource(2).ne."none") then 
          if(obs%selectNlevs.ge.1) then 
             allocate(stats%anomaly%obs_value_total(LVT_rc%ngrid, obs%vlevels, &
                  LVT_rc%strat_nlevels))
             stats%anomaly%obs_value_total = 0.0
             allocate(stats%anomaly%count_obs_value_total(LVT_rc%ngrid, obs%vlevels, &
                  LVT_rc%strat_nlevels))
             stats%anomaly%count_obs_value_total = 0
             allocate(stats%anomaly%obs_value_ci(obs%vlevels,LVT_rc%strat_nlevels))
             stats%anomaly%obs_value_ci = LVT_rc%udef
             
             if(metric%computeSC.eq.1) then 
                allocate(stats%anomaly%obs_value_asc(LVT_rc%ngrid, obs%vlevels,&
                     LVT_rc%nasc))
                stats%anomaly%obs_value_asc = 0.0
                allocate(stats%anomaly%count_obs_value_asc(LVT_rc%ngrid, obs%vlevels,&
                     LVT_rc%nasc))
                stats%anomaly%count_obs_value_asc = 0
             endif
          endif
       endif
    endif
    if(metric%timeOpt.eq.1) then 
       allocate(stats%anomaly%model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%anomaly%model_value_ts = 0.0
       allocate(stats%anomaly%count_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%anomaly%count_model_value_ts = 0

       allocate(stats%anomaly%tavg_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%anomaly%tavg_model_value_ts = 0.0
       allocate(stats%anomaly%tavg_count_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%anomaly%tavg_count_model_value_ts = 0
       
       if(LVT_rc%obssource(2).ne."none") then 
          if(obs%selectNlevs.ge.1) then 
             allocate(stats%anomaly%obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%anomaly%obs_value_ts = 0.0
             allocate(stats%anomaly%count_obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%anomaly%count_obs_value_ts = 0

             allocate(stats%anomaly%tavg_obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%anomaly%tavg_obs_value_ts = 0.0
             allocate(stats%anomaly%tavg_count_obs_value_ts(LVT_rc%ngrid, &
                  model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%anomaly%tavg_count_obs_value_ts = 0

          endif
       endif
    endif

!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 2   
    if(LVT_rc%obssource(2).ne."none") then 
       metric%obsData = .true. 
    else
       metric%obsData = .false. 
    endif
    metric%stdevFlag = .false.

  end subroutine LVT_initAnomaly

!BOP
! 
! !ROUTINE: LVT_diagnoseAnomaly
! \label{LVT_diagnoseAnomaly}
!
! !INTERFACE: 
  subroutine LVT_diagnoseAnomaly(pass)
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
!   calculating the anomaly of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelAnomaly](\ref{diagnoseSingleModelAnomaly})
!     updates the anomaly computation for a single variable 
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

    
    if(pass.eq.1) then 
       if(LVT_metrics%anomaly%selectOpt.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleClimatology(&
                  obs, model, stats)

             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
       endif
    elseif(pass.eq.2) then 
       if(LVT_metrics%anomaly%selectOpt.eq.1.or.&
            LVT_metrics%anomaly%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleAnomaly(model,obs,stats,&
                  LVT_metrics%anomaly)
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseAnomaly

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
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry) :: stats
!EOP
    integer    :: t,k,tind

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
       if(LVT_rc%anomalyTlength.eq.12) then 
          tind = LVT_rc%mo
       else
          tind = 1
       endif
       
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             if(trim(obs%units).eq.trim(model%units)) then 
                if(model%count(t,k).ne.0) then 
                   stats%anomaly%model_value_climo(t,k,tind,1) = &
                        stats%anomaly%model_value_climo(t,k,tind,1) + &
                        model%value(t,k)
                   stats%anomaly%count_model_value_climo(t,k,tind,1) = &
                        stats%anomaly%count_model_value_climo(t,k,tind,1) + 1
                endif
                if(obs%count(t,k).ne.0) then  
                   stats%anomaly%obs_value_climo(t,k,tind,1) = &
                        stats%anomaly%obs_value_climo(t,k,tind,1) + &
                        obs%value(t,k)
                   stats%anomaly%count_obs_value_climo(t,k,tind,1) = &
                        stats%anomaly%count_obs_value_climo(t,k,tind,1) + 1
                endif
                if(LVT_rc%strat_nlevels.gt.1) then 
                   if(LVT_stats%strat_var(t,k).gt.&
                        LVT_rc%strat_var_threshold) then 
                      if(model%count(t,k).ne.0.0) then 
                         stats%anomaly%model_value_climo(t,k,tind,2) = &
                              stats%anomaly%model_value_climo(t,k,tind,2) + &
                              model%value(t,k)
                         stats%anomaly%count_model_value_climo(t,k,tind,2) = &
                              stats%anomaly%count_model_value_climo(t,k,tind,2) + 1
                      endif
                      if(obs%count(t,k).ne.0) then  
                         stats%anomaly%obs_value_climo(t,k,tind,2) = &
                              stats%anomaly%obs_value_climo(t,k,tind,2) + &
                              obs%value(t,k)
                         stats%anomaly%count_obs_value_climo(t,k,tind,2) = &
                              stats%anomaly%count_obs_value_climo(t,k,tind,2) + 1
                      endif
                   elseif(LVT_stats%strat_var(t,k).le.&
                        LVT_rc%strat_var_threshold) then 
                      if(model%count(t,k).ne.0.0) then 
                         stats%anomaly%model_value_climo(t,k,tind,3) = &
                              stats%anomaly%model_value_climo(t,k,tind,3) + &
                              model%value(t,k)
                         stats%anomaly%count_model_value_climo(t,k,tind,3) = &
                              stats%anomaly%count_model_value_climo(t,k,tind,3) + 1
                      endif
                      if(obs%count(t,k).ne.0) then  
                         stats%anomaly%obs_value_climo(t,k,tind,3) = &
                              stats%anomaly%obs_value_climo(t,k,tind,3) + &
                              obs%value(t,k)
                         stats%anomaly%count_obs_value_climo(t,k,tind,3) = &
                              stats%anomaly%count_obs_value_climo(t,k,tind,3) + 1
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
    
  end subroutine diagnoseSingleClimatology

!BOP
! 
! !ROUTINE: diagnoseSingleAnomaly
! \label{diagnoseSingleAnomaly}
!
! !INTERFACE: 
  subroutine diagnoseSingleAnomaly(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the anomaly computation of the 
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
    integer    :: t,k,tind, tval

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then 
       if(LVT_rc%anomalyTlength.eq.12) then 
          tval = LVT_rc%mo
       else
          tval = 1
       endif       
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             if(model%count(t,k).gt.0) then
                if(metric%selectOpt.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      stats%anomaly%model_value_total(t,k,1) = &
                           stats%anomaly%model_value_total(t,k,1) + &
                           (model%value(t,k) - &
                           stats%anomaly%model_value_climo(t,k,tval,1))
                      stats%anomaly%count_model_value_total(t,k,1) = &
                           stats%anomaly%count_model_value_total(t,k,1) + 1
                   endif
                endif
                
                if(metric%timeOpt.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      stats%anomaly%model_value_ts(t,k,1) = &
                           stats%anomaly%model_value_ts(t,k,1)+&
                           (model%value(t,k) - & 
                           stats%anomaly%model_value_climo(t,k,tval,1))
                      stats%anomaly%count_model_value_ts(t,k,1) = & 
                           stats%anomaly%count_model_value_ts(t,k,1)+1
                   endif
                endif
                if(metric%computeSC.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                           tind)
                      stats%anomaly%model_value_asc(t,k,tind) = &
                           stats%anomaly%model_value_asc(t,k,tind)+&
                           (model%value(t,k) - & 
                           stats%anomaly%model_value_climo(t,k,tval,1))
                      stats%anomaly%count_model_value_asc(t,k,tind) = &
                           stats%anomaly%count_model_value_asc(t,k,tind) + 1
                   endif
                endif
                if(LVT_rc%strat_nlevels.gt.1) then 
                   if(LVT_stats%strat_var(t,k).gt.&
                        LVT_rc%strat_var_threshold) then
                      if(metric%selectOpt.eq.1) then
                         if(model%value(t,k).ne.LVT_rc%udef) then  
                            stats%anomaly%model_value_total(t,k,2) = &
                                 stats%anomaly%model_value_total(t,k,2) + &
                                 (model%value(t,k) - & 
                                 stats%anomaly%model_value_climo(t,k,tval,2))
                            stats%anomaly%count_model_value_total(t,k,2) = &
                                 stats%anomaly%count_model_value_total(t,k,2) + 1
                         endif
                      endif
                      if(metric%timeOpt.eq.1) then 
                         if(model%value(t,k).ne.LVT_rc%udef) then 
                            stats%anomaly%model_value_ts(t,k,2) = &
                                 stats%anomaly%model_value_ts(t,k,2)+&
                                 (model%value(t,k)- & 
                                 stats%anomaly%model_value_climo(t,k,tval,2))
                            stats%anomaly%count_model_value_ts(t,k,2) = & 
                                 stats%anomaly%count_model_value_ts(t,k,2)+1
                         endif
                      endif
                   elseif(LVT_stats%strat_var(t,k).le.&
                        LVT_rc%strat_var_threshold) then
                      if(metric%selectOpt.eq.1) then 
                         if(model%value(t,k).ne.LVT_rc%udef) then 
                            stats%anomaly%model_value_total(t,k,3) = &
                                 stats%anomaly%model_value_total(t,k,3) + &
                                 (model%value(t,k) - & 
                                 stats%anomaly%model_value_climo(t,k,tval,3))
                            stats%anomaly%count_model_value_total(t,k,3) = &
                                 stats%anomaly%count_model_value_total(t,k,3) + 1
                         endif
                      endif
                      if(metric%timeOpt.eq.1) then 
                         if(model%value(t,k).ne.LVT_rc%udef) then 
                            stats%anomaly%model_value_ts(t,k,3) = &
                                 stats%anomaly%model_value_ts(t,k,3)+&
                                 (model%value(t,k) - & 
                                 stats%anomaly%model_value_climo(t,k,tval,3))
                            stats%anomaly%count_model_value_ts(t,k,3) = & 
                                 stats%anomaly%count_model_value_ts(t,k,3)+1
                         endif
                      endif
                   endif
                endif
             endif

             if(LVT_rc%obssource(2).ne."none") then
                if(obs%selectNlevs.ge.1) then
                   if(obs%count(t,k).gt.0) then 
                      if(metric%selectOpt.eq.1.and.&
                           obs%value(t,k).ne.LVT_rc%udef) then 
                         stats%anomaly%obs_value_total(t,k,1) = &
                              stats%anomaly%obs_value_total(t,k,1) + &
                              (obs%value(t,k) - & 
                              stats%anomaly%obs_value_climo(t,k,tval,1))
                         stats%anomaly%count_obs_value_total(t,k,1) = &
                              stats%anomaly%count_obs_value_total(t,k,1) + 1
                      endif
                      
                      if(metric%timeOpt.eq.1.and.obs%value(t,k).ne.LVT_rc%udef) then 
                         stats%anomaly%obs_value_ts(t,k,1) = &
                              stats%anomaly%obs_value_ts(t,k,1)+&
                              (obs%value(t,k) - & 
                              stats%anomaly%obs_value_climo(t,k,tval,1))
                         stats%anomaly%count_obs_value_ts(t,k,1) = & 
                              stats%anomaly%count_obs_value_ts(t,k,1)+1
                      endif
                      if(metric%computeSC.eq.1.and.obs%value(t,k).ne.LVT_rc%udef) then 
                         call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,tind)
                         stats%anomaly%obs_value_asc(t,k,tind) = &
                              stats%anomaly%obs_value_asc(t,k,tind)+&
                              (obs%value(t,k) - & 
                              stats%anomaly%obs_value_climo(t,k,tval,1))
                         stats%anomaly%count_obs_value_asc(t,k,tind) = &
                              stats%anomaly%count_obs_value_asc(t,k,tind) + 1
                      endif
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then
                            if(metric%selectOpt.eq.1 &
                                 .and.obs%value(t,k).ne.LVT_rc%udef) then 
                               stats%anomaly%obs_value_total(t,k,2) = &
                                    stats%anomaly%obs_value_total(t,k,2) + &
                                    (obs%value(t,k) - & 
                                    stats%anomaly%obs_value_climo(t,k,tval,2))
                               stats%anomaly%count_obs_value_total(t,k,2) = &
                                    stats%anomaly%count_obs_value_total(t,k,2) + 1
                            endif
                            
                            if(metric%timeOpt.eq.1.and.&
                                 obs%value(t,k).ne.LVT_rc%udef) then 
                               stats%anomaly%obs_value_ts(t,k,2) = &
                                    stats%anomaly%obs_value_ts(t,k,2)+&
                                    (obs%value(t,k) - & 
                                    stats%anomaly%obs_value_climo(t,k,tval,2))
                               stats%anomaly%count_obs_value_ts(t,k,2) = & 
                                    stats%anomaly%count_obs_value_ts(t,k,2)+1
                            endif
                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then
                            if(metric%selectOpt.eq.1.and.&
                                 obs%value(t,k).ne.LVT_rc%udef) then 
                               stats%anomaly%obs_value_total(t,k,3) = &
                                    stats%anomaly%obs_value_total(t,k,3) + &
                                    (obs%value(t,k) - & 
                                    stats%anomaly%obs_value_climo(t,k,tval,3))
                               stats%anomaly%count_obs_value_total(t,k,3) = &
                                    stats%anomaly%count_obs_value_total(t,k,3) + 1
                            endif
                            
                            if(metric%timeOpt.eq.1.and.&
                                 obs%value(t,k).ne.LVT_rc%udef) then 
                               stats%anomaly%obs_value_ts(t,k,3) = &
                                    stats%anomaly%obs_value_ts(t,k,3)+&
                                    (obs%value(t,k) - & 
                                    stats%anomaly%obs_value_climo(t,k,tval,3))
                               stats%anomaly%count_obs_value_ts(t,k,3) = & 
                                    stats%anomaly%count_obs_value_ts(t,k,3)+1
                            endif
                            
                         endif
                      endif
                   endif
                endif
             endif
          enddo
       enddo
    endif
  end subroutine diagnoseSingleAnomaly

!BOP
! 
! !ROUTINE: LVT_computeAnomaly
! \label{LVT_computeAnomaly}
!
! !INTERFACE: 
  subroutine LVT_computeAnomaly(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the anomaly values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelAnomaly](\ref{computeSingleModelAnomaly})
!     computes the anomaly values for a single variable
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
    integer     :: i 
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats
    
    if(pass.eq.1) then 
       if(LVT_metrics%anomaly%selectOpt.eq.1.or.&
            LVT_metrics%anomaly%timeOpt.eq.1) then

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

       endif
    elseif(pass.eq.2) then 
       if(LVT_metrics%anomaly%selectOpt.eq.1.or.&
            LVT_metrics%anomaly%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%anomaly%timeOpt.eq.1.and.&
                  LVT_metrics%anomaly%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%anomaly%ftn_ts_loc(i),200,advance='no') &
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
             
             call computeSingleAnomaly(alarm,model,obs,stats,&
                  LVT_metrics%anomaly)
             
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%anomaly%timeOpt.eq.1.and.&
                  LVT_metrics%anomaly%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%anomaly%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_computeAnomaly
  

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
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do m=1,LVT_rc%anomalyTlength
                   do l=1, LVT_rc%strat_nlevels
                      if(stats%anomaly%count_model_value_climo(t,k,m,l).ne.0) then 
                         stats%anomaly%model_value_climo(t,k,m,l) = &
                              stats%anomaly%model_value_climo(t,k,m,l) /&
                              stats%anomaly%count_model_value_climo(t,k,m,l)
                      endif
                      if(stats%anomaly%count_obs_value_climo(t,k,m,l).ne.0) then  
                         stats%anomaly%obs_value_climo(t,k,m,l) = &
                              stats%anomaly%obs_value_climo(t,k,m,l)/&
                              stats%anomaly%count_obs_value_climo(t,k,m,l)
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
! !ROUTINE: computeSingleAnomaly
! \label{computeSingleAnomaly}
!
! !INTERFACE: 
  subroutine computeSingleAnomaly(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the anomaly values
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
    real     :: diff_field(LVT_rc%ngrid)

    if(metric%timeOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%anomaly%count_model_value_ts(t,k,l).gt.0) then 
                      stats%anomaly%model_value_ts(t,k,l) =&
                           stats%anomaly%model_value_ts(t,k,l)/&
                           stats%anomaly%count_model_value_ts(t,k,l) 

                      stats%anomaly%tavg_model_value_ts(t,k,l) =&
                           stats%anomaly%tavg_model_value_ts(t,k,l) + &
                           stats%anomaly%model_value_ts(t,k,l)
                      stats%anomaly%tavg_count_model_value_ts(t,k,l) =&
                           stats%anomaly%tavg_count_model_value_ts(t,k,l) + 1
                   else
                      stats%anomaly%model_value_ts(t,k,l) = LVT_rc%udef
                   endif
                   if(LVT_rc%obssource(2).ne."none") then 
                      if(obs%selectNlevs.ge.1) then 
                         if(stats%anomaly%count_obs_value_ts(t,k,l).gt.0) then 
                            stats%anomaly%obs_value_ts(t,k,l) = &
                                 stats%anomaly%obs_value_ts(t,k,l)/&
                                 stats%anomaly%count_obs_value_ts(t,k,l)

                            stats%anomaly%tavg_obs_value_ts(t,k,l) = &
                                 stats%anomaly%tavg_obs_value_ts(t,k,l) + &
                                 stats%anomaly%obs_value_ts(t,k,l)
                            stats%anomaly%tavg_count_obs_value_ts(t,k,l) = &
                                 stats%anomaly%tavg_count_obs_value_ts(t,k,l)+ 1
                         else
                            stats%anomaly%obs_value_ts(t,k,l) = LVT_rc%udef
                         endif
                      endif
                   endif
                enddo
             enddo
          enddo
          
          if(alarm) then 

             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%anomaly%tavg_count_model_value_ts(t,k,l).gt.0) then 
                         stats%anomaly%tavg_model_value_ts(t,k,l) =&
                              stats%anomaly%tavg_model_value_ts(t,k,l)/&
                              stats%anomaly%tavg_count_model_value_ts(t,k,l) 
                      else
                         stats%anomaly%tavg_model_value_ts(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                enddo
             enddo
             if(LVT_rc%obssource(2).ne."none") then 
                do t=1,LVT_rc%ngrid
                   do k=1,obs%selectNlevs
                      do l=1,LVT_rc%strat_nlevels
                         if(stats%anomaly%tavg_count_obs_value_ts(t,k,l).gt.0) then 
                            stats%anomaly%tavg_obs_value_ts(t,k,l) =&
                                 stats%anomaly%tavg_obs_value_ts(t,k,l)/&
                                 stats%anomaly%tavg_count_obs_value_ts(t,k,l) 
                         else
                            stats%anomaly%tavg_obs_value_ts(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   enddo
                enddo
             endif

             if(metric%extractTS.eq.1) then 
                if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                   call LVT_writeTSinfo(metric%ftn_ts_loc,&
                        model,&
                        LVT_rc%ngrid,&
                        stats%anomaly%tavg_model_value_ts,&
                        stats%anomaly%tavg_count_model_value_ts,&
                        LVT_rc%ngrid,&
                        stats%anomaly%tavg_obs_value_ts,&
                        stats%anomaly%tavg_count_obs_value_ts)
                else
                   call LVT_writeTSinfo(metric%ftn_ts_loc,&
                        model,&
                        LVT_rc%ngrid,&
                        stats%anomaly%tavg_model_value_ts,&
                        stats%anomaly%tavg_count_model_value_ts)
                endif
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
                   if(stats%anomaly%count_model_value_total(t,k,l).gt.&
                        LVT_rc%obsCountThreshold) then 
                      stats%anomaly%model_value_total(t,k,l) = &
                           stats%anomaly%model_value_total(t,k,l)/&
                              stats%anomaly%count_model_value_total(t,k,l)           
                   else
                      stats%anomaly%model_value_total(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1, LVT_rc%nasc
                      if(stats%anomaly%count_model_value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then 
                         stats%anomaly%model_value_asc(t,k,l) = &
                              stats%anomaly%model_value_asc(t,k,l)/&
                              stats%anomaly%count_model_value_asc(t,k,l)           
                      else
                         stats%anomaly%model_value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif

                if(LVT_rc%obssource(2).ne."none") then 
                   if(obs%selectNlevs.ge.1) then
                      do l=1,LVT_rc%strat_nlevels                      
                         if(stats%anomaly%count_obs_value_total(t,k,l).gt.&
                              LVT_rc%obsCountThreshold) then 
                            stats%anomaly%obs_value_total(t,k,l) = &
                                 stats%anomaly%obs_value_total(t,k,l)/&
                                 stats%anomaly%count_obs_value_total(t,k,l)       
                         else
                            stats%anomaly%obs_value_total(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                      if(metric%computeSC.eq.1) then
                         do l=1, LVT_rc%nasc 
                            if(stats%anomaly%count_obs_value_asc(t,k,l).gt.&
                                 LVT_rc%SCCountThreshold) then 
                               stats%anomaly%obs_value_asc(t,k,l) = &
                                    stats%anomaly%obs_value_asc(t,k,l)/&
                                    stats%anomaly%count_obs_value_asc(t,k,l)           
                            else
                               stats%anomaly%obs_value_asc(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                      endif
                      
                   endif
                endif
             enddo
          enddo
          
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%anomaly%model_value_total(:,k,l),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%anomaly%model_value_ci(k,l))
             enddo
          enddo
          if(LVT_rc%obssource(2).ne."none") then 
             do k=1,obs%vlevels
                do l=1, LVT_rc%strat_nlevels
                   call LVT_computeCI(stats%anomaly%obs_value_total(:,k,l),LVT_rc%ngrid,&
                        LVT_rc%pval_CI,stats%anomaly%obs_value_ci(k,l))
                enddo
             enddo
          endif
          if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
             call LVT_writeDataBasedStrat(model,obs,stats,metric,&
                  LVT_rc%ngrid, stats%anomaly%model_value_total, LVT_rc%ngrid,&
                  stats%anomaly%obs_value_total)
             
             if(metric%computeSC.eq.1) then 
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%anomaly%model_value_asc,stats%anomaly%count_model_value_asc,&
                     LVT_rc%ngrid,stats%anomaly%obs_value_asc,stats%anomaly%count_obs_value_asc)          
             endif
          else
             call LVT_writeDataBasedStrat(model,obs,stats,metric,&
                   LVT_rc%ngrid,stats%anomaly%model_value_total)

             if(metric%computeSC.eq.1) then 
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%anomaly%model_value_asc,stats%anomaly%count_model_value_asc)
             endif
          endif
       endif
    endif
  end subroutine computeSingleAnomaly

!BOP
! 
! !ROUTINE: LVT_writeMetric_Anomaly
! \label(LVT_writeMetric_Anomaly)
!
! !INTERFACE:
  subroutine LVT_writeMetric_Anomaly(pass,final,vlevels,stats,obs)
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

    if(pass.eq.LVT_metrics%anomaly%npass) then 
       if(final.ne.1) then
           if(stats%selectOpt.eq.1) then 
              do k=1,vlevels
                if(LVT_metrics%anomaly%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_gridded(LVT_metrics%anomaly%ftn_ts, &
                           stats%anomaly%model_value_ts(:,k,l),&
                           stats%vid_ts(LVT_Anomalyid,1),k)

                      call LVT_writevar_gridded(LVT_metrics%anomaly%ftn_ts, &
                           real(stats%anomaly%count_model_value_ts(:,k,l)),&
                           stats%vid_count_ts(LVT_Anomalyid,1),k)
                   enddo
                   if(LVT_rc%obssource(2).ne."none") then 
                      do l=1,LVT_rc%strat_nlevels
                         call LVT_writevar_gridded(LVT_metrics%anomaly%ftn_ts, &
                              stats%anomaly%obs_value_ts(:,k,l),&
                              stats%vid_ts(LVT_Anomalyid,2),k)
                         call LVT_writevar_gridded(LVT_metrics%anomaly%ftn_ts, &
                              real(stats%anomaly%count_obs_value_ts(:,k,l)),&
                              stats%vid_count_ts(LVT_Anomalyid,2),k)
                      enddo
                   endif
                endif
             enddo
          endif
       else
          if(pass.eq.1) then 
             if(stats%selectOpt.eq.1) then 
                do k=1,vlevels
                   if(LVT_metrics%anomaly%selectOpt.eq.1) then 
                      do l=1,LVT_rc%strat_nlevels                      
                         call LVT_writevar_gridded(LVT_metrics%anomaly%ftn_total, &
                              stats%anomaly%model_value_total(:,k,l),&
                              stats%vid_total(LVT_Anomalyid,1),k)
                         call LVT_writevar_gridded(LVT_metrics%anomaly%ftn_total, &
                              real(stats%anomaly%count_model_value_total(:,k,l)),&
                              stats%vid_count_total(LVT_Anomalyid,1),k)
                         
                      enddo
                      if(LVT_rc%obssource(2).ne."none".and.&
                           obs%selectNlevs.ge.1) then 
                         do l=1,LVT_rc%strat_nlevels                      
                            call LVT_writevar_gridded(LVT_metrics%anomaly%ftn_total, &
                                 stats%anomaly%obs_value_total(:,k,l),&
                                 stats%vid_total(LVT_Anomalyid,2),k)
                            call LVT_writevar_gridded(LVT_metrics%anomaly%ftn_total, &
                                 real(stats%anomaly%count_obs_value_total(:,k,l)),&
                                 stats%vid_count_total(LVT_Anomalyid,2),k)                         
                         enddo
                      endif
                   
                      if(LVT_metrics%anomaly%computeSC.eq.1) then 
                         do tind = 1,LVT_rc%nasc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%anomaly%ftn_total,&
                                 stats%anomaly%model_value_asc(:,k,tind),&
                                 stats%vid_sc_total(tind,LVT_Anomalyid,1),k)
                         enddo
                         if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                            do tind = 1,LVT_rc%nasc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%anomaly%ftn_total,&
                                    stats%anomaly%obs_value_asc(:,k,tind),&
                                    stats%vid_sc_total(tind,LVT_Anomalyid,2),k)
                            enddo
                         endif
                      endif
                      call LVT_writeSummaryStats(&
                           LVT_metrics%anomaly%ftn_summ,&
                           LVT_metrics%anomaly%short_name,&
                           LVT_rc%ngrid,&
                           stats%anomaly%model_value_total(:,k,:), &
                           stats%anomaly%count_model_value_total(:,k,:),&
                           stats%standard_name,&
                           stats%anomaly%model_value_ci(k,:))
                      if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                         call LVT_writeSummaryStats(&
                              LVT_metrics%anomaly%ftn_summ,&
                              LVT_metrics%anomaly%short_name,&
                              LVT_rc%ngrid,&
                              stats%anomaly%obs_value_total(:,k,:), &
                              stats%anomaly%count_obs_value_total(:,k,:),&
                              "OBS_"//trim(stats%standard_name),&
                              stats%anomaly%obs_value_ci(k,:))
                      endif
                   endif
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_writeMetric_Anomaly

!BOP
! 
! !ROUTINE: LVT_resetMetric_Anomaly
! \label(LVT_resetMetric_Anomaly)
!
! !INTERFACE:
  subroutine LVT_resetMetric_Anomaly(alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
    logical         :: alarm
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
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats


    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream1Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(stats%selectOpt.eq.1) then 
          do k=1,model%selectNlevs
             if(LVT_metrics%anomaly%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   stats%anomaly%model_value_ts(:,k,l) = 0.0
                   stats%anomaly%count_model_value_ts(:,k,l)=0 
                   if(alarm) then 
                      stats%anomaly%tavg_model_value_ts(:,k,l) = 0.0
                      stats%anomaly%tavg_count_model_value_ts(:,k,l)=0 
                   endif
                enddo
                if(LVT_rc%obssource(2).ne."none".and.&
                     obs%selectNlevs.ge.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      stats%anomaly%obs_value_ts(:,k,l) = 0.0
                      stats%anomaly%count_obs_value_ts(:,k,l)=0 
                      if(alarm) then 
                         stats%anomaly%tavg_obs_value_ts(:,k,l) = 0.0
                         stats%anomaly%tavg_count_obs_value_ts(:,k,l)=0 
                      endif
                   enddo
                endif
             endif
          enddo
       endif
       
       model => model%next
       obs => obs%next
       stats => stats%next
    enddo

  end subroutine LVT_resetMetric_Anomaly

!BOP
! 
! !ROUTINE: LVT_writerestart_Anomaly
! 
! !INTERFACE:
  subroutine LVT_writerestart_Anomaly(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for Anomaly metric computations
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
       if(LVT_metrics%anomaly%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   call LVT_writevar_restart(ftn,&
                        stats%anomaly%model_value_total(:,k,l))
                   call LVT_writevar_restart(ftn,&
                        stats%anomaly%count_model_value_total(:,k,l))
                enddo
             enddo
             if(LVT_metrics%anomaly%computeSC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc         
                      call LVT_writevar_restart(ftn,&
                           stats%anomaly%model_value_asc(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%anomaly%count_model_value_asc(:,k,l))
                   enddo
                enddo
             endif
             
             if(LVT_rc%obssource(2).ne."none") then 
                if(obs%selectNlevs.ge.1) then 
                   do k=1,obs%vlevels
                      do l=1,LVT_rc%strat_nlevels
                         call LVT_writevar_restart(ftn,&
                              stats%anomaly%obs_value_total(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%anomaly%count_obs_value_total(:,k,l))
                      enddo
                   enddo
                   
                   if(LVT_metrics%anomaly%computeSC.eq.1) then 
                      do k=1,obs%vlevels
                         do l=1,LVT_rc%nasc
                            call LVT_writevar_restart(ftn,&
                                 stats%anomaly%obs_value_asc(:,k,l))
                            call LVT_writevar_restart(ftn,&
                                 stats%anomaly%count_obs_value_asc(:,k,l))
                         enddo
                      enddo
                   endif
                endif
             endif
          endif
       endif
       model => model%next
       obs   => obs%next
       stats => stats%next
    end do
  end subroutine LVT_writerestart_Anomaly


!BOP
! 
! !ROUTINE: LVT_readrestart_Anomaly
! 
! !INTERFACE:
  subroutine LVT_readrestart_Anomaly(ftn)
! !USES: 
! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for Anomaly metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP

    integer              :: k,l,index
    type(LVT_metaDataEntry), pointer :: model
    type(LVT_metaDataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats

    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(LVT_metrics%anomaly%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 

             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   call LVT_readvar_restart(ftn,&
                        stats%anomaly%model_value_total(:,k,l))
                   call LVT_readvar_restart(ftn,&
                        stats%anomaly%count_model_value_total(:,k,l))
                enddo
             enddo
             if(LVT_metrics%anomaly%computeSC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc         
                      call LVT_readvar_restart(ftn,&
                           stats%anomaly%model_value_asc(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%anomaly%count_model_value_asc(:,k,l))
                   enddo
                enddo
             endif
             
             if(LVT_rc%obssource(2).ne."none") then 
                if(obs%selectNlevs.ge.1) then 
                   do k=1,obs%vlevels
                      do l=1,LVT_rc%strat_nlevels
                         call LVT_readvar_restart(ftn,&
                              stats%anomaly%obs_value_total(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%anomaly%count_obs_value_total(:,k,l))
                      enddo
                   enddo
                   
                   if(LVT_metrics%anomaly%computeSC.eq.1) then 
                      do k=1,obs%vlevels
                         do l=1,LVT_rc%nasc
                            call LVT_readvar_restart(ftn,&
                                 stats%anomaly%obs_value_asc(:,k,l))
                            call LVT_readvar_restart(ftn,&
                                 stats%anomaly%count_obs_value_asc(:,k,l))
                         enddo
                      enddo
                   endif
                endif
             endif
          endif
       endif
       model => model%next
       obs   => obs%next
       stats => stats%next
    end do
  end subroutine LVT_readrestart_Anomaly


end module LVT_AnomalyMod
