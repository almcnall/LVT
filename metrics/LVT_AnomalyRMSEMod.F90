!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_AnomalyRMSEMod
! \label(LVT_AnomalyRMSEMod)
!
! !INTERFACE:
module LVT_AnomalyRMSEMod
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
!   the LIS output to the specified observations. 
!
! Time series, seasonal cycles, average diurnal cycle options are 
! not supported for Anomaly RMSEelation computation
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
  public :: LVT_initAnomalyRMSE
  public :: LVT_diagnoseAnomalyRMSE
  public :: LVT_computeAnomalyRMSE 
  public :: LVT_writeMetric_AnomalyRMSE
  public :: LVT_resetMetric_AnomalyRMSE
  public :: LVT_writerestart_AnomalyRMSE
  public :: LVT_readrestart_AnomalyRMSE
contains
  
!BOP
! 
! !ROUTINE: LVT_initAnomalyRMSE
! \label{LVT_initAnomalyRMSE}
!
! !INTERFACE: 
  subroutine LVT_initAnomalyRMSE(model,obs,stats,metric)
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
       allocate(stats%armse%model_value_climo(LVT_rc%ngrid,&
            model%selectNlevs, 12, LVT_rc%strat_nlevels))
       allocate(stats%armse%obs_value_climo(LVT_rc%ngrid,&
            model%selectNlevs, 12,LVT_rc%strat_nlevels))
       allocate(stats%armse%count_model_value_climo(LVT_rc%ngrid,&
            model%selectNlevs, 12,LVT_rc%strat_nlevels))
       allocate(stats%armse%count_obs_value_climo(LVT_rc%ngrid,&
            model%selectNlevs, 12,LVT_rc%strat_nlevels))
              
       stats%armse%model_value_climo = 0 
       stats%armse%obs_value_climo   = 0 
       stats%armse%count_model_value_climo = 0 
       stats%armse%count_obs_value_climo   = 0 
       
       allocate(stats%armse%value_total(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))
       stats%armse%value_total = 0.0
       allocate(stats%armse%count_value_total(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       
       stats%armse%count_value_total = 0
       allocate(stats%armse%value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%armse%value_ci = LVT_rc%udef

       if(metric%timeopt.eq.1) then 
          allocate(stats%armse%value_ts(LVT_rc%ngrid, model%selectNlevs, &
               LVT_rc%strat_nlevels))
          stats%armse%value_ts = 0.0
          allocate(stats%armse%count_value_ts(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%strat_nlevels))       
          stats%armse%count_value_ts = 0 

          allocate(stats%armse%tavg_value_ts(LVT_rc%ngrid, model%selectNlevs, &
               LVT_rc%strat_nlevels))
          stats%armse%tavg_value_ts = 0.0
          allocate(stats%armse%tavg_count_value_ts(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%strat_nlevels))       
          stats%armse%tavg_count_value_ts = 0 

          if(metric%computeSC.eq.1) then 
             allocate(stats%armse%value_asc(LVT_rc%ngrid, model%selectNlevs,LVT_rc%nasc))
             stats%armse%value_asc = 0.0
             allocate(stats%armse%count_value_asc(LVT_rc%ngrid, model%selectNlevs,&
                  LVT_rc%nasc))
             stats%armse%count_value_asc = 0
          endif
          if(metric%computeADC.eq.1) then 
             allocate(stats%armse%value_adc(LVT_rc%ngrid, model%selectNlevs,LVT_rc%nadc))
             stats%armse%value_adc = 0.0
             allocate(stats%armse%count_value_adc(LVT_rc%ngrid, model%selectNlevs,&
                  LVT_rc%nadc))
             stats%armse%count_value_adc = 0
          endif
       endif
    endif
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 2
    metric%obsData = .false. 
    metric%stdevFlag = .false. 
  end subroutine LVT_initAnomalyRMSE
  
!BOP
! 
! !ROUTINE: LVT_diagnoseAnomalyRMSE
! \label{LVT_diagnoseAnomalyRMSE}
!
! !INTERFACE: 
  subroutine LVT_diagnoseAnomalyRMSE(pass)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the AnomalyRMSE calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleAnomalyRMSE](\ref{diagnoseSingleAnomalyRMSE})
!     updates the AnomalyRMSE computation for a single variable 
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
       if(LVT_metrics%armse%selectOpt.eq.1) then 
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
       if(LVT_metrics%armse%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call diagnoseSingleAnomalyRMSE(&
                  obs, model, stats, &
                  LVT_metrics%armse)

             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
       endif
    endif
  end subroutine LVT_diagnoseAnomalyRMSE

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
    integer    :: t,k

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             if(trim(obs%units).eq.trim(model%units)) then 
                if(model%count(t,k).ne.0.and.obs%count(t,k).ne.0) then 
                   stats%armse%model_value_climo(t,k,LVT_rc%mo,1) = &
                        stats%armse%model_value_climo(t,k,LVT_rc%mo,1) + &
                        model%value(t,k)
                   stats%armse%count_model_value_climo(t,k,LVT_rc%mo,1) = &
                        stats%armse%count_model_value_climo(t,k,LVT_rc%mo,1) + 1
                endif
                if(obs%count(t,k).ne.0) then  
                   stats%armse%obs_value_climo(t,k,LVT_rc%mo,1) = &
                        stats%armse%obs_value_climo(t,k,LVT_rc%mo,1) + &
                        obs%value(t,k)
                   stats%armse%count_obs_value_climo(t,k,LVT_rc%mo,1) = &
                        stats%armse%count_obs_value_climo(t,k,LVT_rc%mo,1) + 1
                endif
                if(LVT_rc%strat_nlevels.gt.1) then 
                   if(LVT_stats%strat_var(t,k).gt.&
                        LVT_rc%strat_var_threshold) then 
                      if(model%count(t,k).ne.0.and.obs%count(t,k).ne.0) then 
                         stats%armse%model_value_climo(t,k,LVT_rc%mo,2) = &
                              stats%armse%model_value_climo(t,k,LVT_rc%mo,2) + &
                              model%value(t,k)
                         stats%armse%count_model_value_climo(t,k,LVT_rc%mo,2) = &
                              stats%armse%count_model_value_climo(t,k,LVT_rc%mo,2) + 1
                      endif
                      if(obs%count(t,k).ne.0) then  
                         stats%armse%obs_value_climo(t,k,LVT_rc%mo,2) = &
                              stats%armse%obs_value_climo(t,k,LVT_rc%mo,2) + &
                              obs%value(t,k)
                         stats%armse%count_obs_value_climo(t,k,LVT_rc%mo,2) = &
                              stats%armse%count_obs_value_climo(t,k,LVT_rc%mo,2) + 1
                      endif
                   elseif(LVT_stats%strat_var(t,k).le.&
                        LVT_rc%strat_var_threshold) then 
                      if(model%count(t,k).ne.0.and.obs%count(t,k).ne.0) then 
                         stats%armse%model_value_climo(t,k,LVT_rc%mo,3) = &
                              stats%armse%model_value_climo(t,k,LVT_rc%mo,3) + &
                              model%value(t,k)
                         stats%armse%count_model_value_climo(t,k,LVT_rc%mo,3) = &
                              stats%armse%count_model_value_climo(t,k,LVT_rc%mo,3) + 1
                      endif
                      if(obs%count(t,k).ne.0) then  
                         stats%armse%obs_value_climo(t,k,LVT_rc%mo,3) = &
                              stats%armse%obs_value_climo(t,k,LVT_rc%mo,3) + &
                              obs%value(t,k)
                         stats%armse%count_obs_value_climo(t,k,LVT_rc%mo,3) = &
                              stats%armse%count_obs_value_climo(t,k,LVT_rc%mo,3) + 1
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
! !ROUTINE: diagnoseSingleAnomalyRMSE
! \label{diagnoseSingleAnomalyRMSE}
!
! !INTERFACE: 
  subroutine diagnoseSingleAnomalyRMSE(obs, model, stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the AnomalyRMSE computation (updates the running 
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
                      stats%armse%value_total(t,k,1) = stats%armse%value_total(t,k,1) + &
                           ((obs%value(t,k)-stats%armse%obs_value_climo(t,k,LVT_rc%mo,1)) &
                           -(model%value(t,k)-stats%armse%model_value_climo(t,k,LVT_rc%mo,1)))**2 
                      stats%armse%count_value_total(t,k,1) = &
                           stats%armse%count_value_total(t,k,1) + 1
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then 
                            stats%armse%value_total(t,k,2) = stats%armse%value_total(t,k,2) + &
                                 ((obs%value(t,k)-stats%armse%obs_value_climo(t,k,LVT_rc%mo,2)) &
                                 -(model%value(t,k)-stats%armse%model_value_climo(t,k,LVT_rc%mo,2)))**2 
                            stats%armse%count_value_total(t,k,2) = &
                                 stats%armse%count_value_total(t,k,2) + 1
                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then 
                            stats%armse%value_total(t,k,3) = stats%armse%value_total(t,k,3) + &
                                 ((obs%value(t,k)-stats%armse%obs_value_climo(t,k,LVT_rc%mo,3)) &
                                 -(model%value(t,k)-stats%armse%model_value_climo(t,k,LVT_rc%mo,3)))**2 
                            stats%armse%count_value_total(t,k,3) = & 
                                 stats%armse%count_value_total(t,k,3) + 1
                         endif
                      endif
                   endif
                   if(metric%timeOpt.eq.1) then 
                      stats%armse%value_ts(t,k,1) = stats%armse%value_ts(t,k,1) + &
                           ((obs%value(t,k)-stats%armse%obs_value_climo(t,k,LVT_rc%mo,1)) &
                           -(model%value(t,k)-stats%armse%model_value_climo(t,k,LVT_rc%mo,1)))**2 
                      stats%armse%count_value_ts(t,k,1) = stats%armse%count_value_ts(t,k,1)+1
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then
                            stats%armse%value_ts(t,k,2) = stats%armse%value_ts(t,k,2) + &
                                 ((obs%value(t,k)-stats%armse%obs_value_climo(t,k,LVT_rc%mo,2)) &
                                 -(model%value(t,k)-stats%armse%model_value_climo(t,k,LVT_rc%mo,2)))**2 
                            stats%armse%count_value_ts(t,k,2) = stats%armse%count_value_ts(t,k,2)+1
                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then 
                            stats%armse%value_ts(t,k,3) = stats%armse%value_ts(t,k,3) + &
                                 ((obs%value(t,k)-stats%armse%obs_value_climo(t,k,LVT_rc%mo,3)) &
                                 -(model%value(t,k)-stats%armse%model_value_climo(t,k,LVT_rc%mo,3)))**2 
                            stats%armse%count_value_ts(t,k,3) = stats%armse%count_value_ts(t,k,3)+1
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
    
  end subroutine diagnoseSingleAnomalyRMSE


!BOP
! 
! !ROUTINE: LVT_computeAnomalyRMSE
! \label{LVT_computeAnomalyRMSE}
!
! !INTERFACE: 
  subroutine LVT_computeAnomalyRMSE(pass,alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute AnomalyRMSE values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleAnomalyRMSE](\ref{computeSingleAnomalyRMSE})
!     updates the AnomalyRMSE computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     AnomalyRMSE computation has been reached
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    implicit none

    integer               :: i,pass

    logical     :: alarm
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats


    if(LVT_metrics%armse%selectOpt.eq.1.or.&
         LVT_metrics%armse%timeOpt.eq.1) then        
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
          if(alarm) then 
             if(LVT_metrics%armse%timeOpt.eq.1.and.&
                  LVT_metrics%armse%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%armse%ftn_ts_loc(i),200,advance='no') &
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
             call computeSingleAnomalyRMSE(alarm,&
                  obs, model, stats, &
                  LVT_metrics%armse)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%armse%timeOpt.eq.1.and.&
                  LVT_metrics%armse%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%armse%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif

  end subroutine LVT_ComputeAnomalyRMSE

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
                do m=1,12
                   do l=1, LVT_rc%strat_nlevels
                      if(stats%armse%count_model_value_climo(t,k,m,l).ne.0) then 
                         stats%armse%model_value_climo(t,k,m,l) = &
                              stats%armse%model_value_climo(t,k,m,l) /&
                              stats%armse%count_model_value_climo(t,k,m,l)
                      endif
                      if(stats%armse%count_obs_value_climo(t,k,m,l).ne.0) then  
                         stats%armse%obs_value_climo(t,k,m,l) = &
                              stats%armse%obs_value_climo(t,k,m,l)/&
                              stats%armse%count_obs_value_climo(t,k,m,l)
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
! !ROUTINE: computeSingleAnomalyRMSE
! \label{computeSingleAnomalyRMSE}
!
! !INTERFACE: 
  subroutine computeSingleAnomalyRMSE(alarm,obs, model,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the AnomalyRMSE values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     AnomalyRMSE computation has been reached
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
                   if(stats%armse%count_value_ts(t,k,l).ne.0) then 
                      stats%armse%value_ts(t,k,l) = sqrt(stats%armse%value_ts(t,k,l)/&
                           stats%armse%count_value_ts(t,k,l))

                      stats%armse%tavg_value_ts(t,k,l) = &
                           stats%armse%tavg_value_ts(t,k,l) + & 
                           stats%armse%value_ts(t,k,l)
                      stats%armse%tavg_count_value_ts(t,k,l) = & 
                           stats%armse%tavg_count_value_ts(t,k,l) + 1
                   else
                      stats%armse%value_ts(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   if(stats%armse%count_value_ts(t,k,1).ne.0) then 
                      call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                           tind)
                      stats%armse%value_asc(t,k,tind) = &
                           stats%armse%value_asc(t,k,tind)+ &
                           stats%armse%value_ts(t,k,1)
                      stats%armse%count_value_asc(t,k,tind) = &
                           stats%armse%count_value_asc(t,k,tind)+ 1
                   endif
                endif
                if(metric%computeADC.eq.1) then 
                   if(stats%armse%count_value_ts(t,k,1).ne.0) then 
                      call LVT_getADCTimeIndex(tind)
                      stats%armse%value_adc(t,k,tind) = &
                           stats%armse%value_adc(t,k,tind)+ &
                           stats%armse%value_ts(t,k,1)
                      stats%armse%count_value_adc(t,k,tind) = &
                           stats%armse%count_value_adc(t,k,tind)+ 1
                   endif                
                endif
             enddo
          enddo
          if(alarm) then 
             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%armse%tavg_count_value_ts(t,k,l).ne.0) then 
                         
                         stats%armse%tavg_value_ts(t,k,l) = &
                              stats%armse%tavg_value_ts(t,k,l) /& 
                              stats%armse%tavg_count_value_ts(t,k,l) 
                      else
                         stats%armse%tavg_value_ts(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                enddo
             enddo
             
             
             if(metric%extractTS.eq.1) then 
                call LVT_writeTSinfo(metric%ftn_ts_loc,&
                     model,&
                     LVT_rc%ngrid,&
                     stats%armse%tavg_value_ts,&
                     stats%armse%tavg_count_value_ts)
             endif
          endif
       endif
    endif


    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%armse%count_value_total(t,k,l).ne.0.and.&
                        stats%armse%count_value_total(t,k,l)&
                        .gt.LVT_rc%obsCountThreshold) then 
                      
                      stats%armse%value_total(t,k,l) = sqrt(stats%armse%value_total(t,k,l)/&
                           stats%armse%count_value_total(t,k,l))                      
                   else
                      stats%armse%value_total(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1,LVT_rc%nasc
                      if(stats%armse%count_value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then
                         stats%armse%value_asc(t,k,l) = sqrt(stats%armse%value_asc(t,k,l)/&
                              stats%armse%count_value_asc(t,k,l))
                      else
                         stats%armse%value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
                
                if(metric%computeADC.eq.1) then 
                   do l=1,LVT_rc%nadc
                      if(stats%armse%count_value_adc(t,k,l).gt.&
                           LVT_rc%ADCCountThreshold) then
                         stats%armse%value_adc(t,k,l) = sqrt(stats%armse%value_adc(t,k,l)/&
                              stats%armse%count_value_adc(t,k,l))
                      else
                         stats%armse%value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
             enddo
          enddo
          
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%armse%value_total(:,k,l),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%armse%value_ci(k,l))
             enddo
          enddo
       
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid,stats%armse%value_total)      
          if(metric%computeSC.eq.1) then
             call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%armse%value_asc,stats%armse%count_value_asc)
          endif
          if(metric%computeADC.eq.1) then 
             call LVT_writeAvgDiurnalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%armse%value_adc,stats%armse%count_value_adc)
          endif
       endif       
    endif

  end subroutine computeSingleAnomalyRMSE

!BOP
! 
! !ROUTINE: LVT_writeMetric_AnomalyRMSE
! \label{LVT_writeMetric_AnomalyRMSE}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_AnomalyRMSE(pass,final,vlevels,stats,obs)
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

    if(pass.eq.LVT_metrics%armse%npass) then 
       if(final.ne.1) then 
          if(stats%selectOpt.eq.1) then 
             do k=1,vlevels
                if(LVT_metrics%armse%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_gridded(LVT_metrics%armse%ftn_ts, &
                           stats%armse%value_ts(:,k,l),stats%vid_ts(LVT_ARMSEid,1),1)
                      call LVT_writevar_gridded(LVT_metrics%armse%ftn_ts, &
                           real(stats%armse%count_value_ts(:,k,l)),&
                           stats%vid_count_ts(LVT_ARMSEid,1),1)
                   enddo
                endif
             enddo
          endif
       else
          if(pass.eq.2) then 
             if(LVT_metrics%armse%selectOpt.eq.1) then
                if(stats%selectOpt.eq.1) then 
                   do k=1,vlevels
                      do l=1,LVT_rc%strat_nlevels
                         
                         call LVT_writevar_gridded(LVT_metrics%armse%ftn_total, &
                              stats%armse%value_total(:,k,l),stats%vid_total(LVT_ARMSEid,1))
                         call LVT_writevar_gridded(LVT_metrics%armse%ftn_total, &
                              real(stats%armse%count_value_total(:,k,l)),&
                              stats%vid_count_total(LVT_ARMSEid,1))
                         
                      enddo
                      if(LVT_metrics%armse%computeSC.eq.1) then 
                         do tind = 1,LVT_rc%nasc
                            call LVT_writevar_gridded(LVT_metrics%armse%ftn_total,&
                                 stats%armse%value_asc(:,k,tind),&
                                 stats%vid_sc_total(tind,LVT_ARMSEid,1))
                         enddo
                      endif
                      
                      if(LVT_metrics%armse%computeADC.eq.1) then 
                         do tind = 1,LVT_rc%nadc
                            call LVT_writevar_gridded(LVT_metrics%armse%ftn_total,&
                                 stats%armse%value_adc(:,k,tind),&
                                 stats%vid_adc_total(tind,LVT_ARMSEid,1))
                         enddo
                      endif
                      call LVT_writeSummaryStats(&
                           LVT_metrics%armse%ftn_summ,&
                           LVT_metrics%armse%short_name,&
                           LVT_rc%ngrid,&
                           stats%armse%value_total(:,k,:), &
                           stats%armse%count_value_total(:,k,:),stats%standard_name,&
                           stats%armse%value_ci(k,:))
                   enddo
                endif
             endif
          endif
       end if
    endif


  end subroutine LVT_writeMetric_AnomalyRMSE

!BOP
! 
! !ROUTINE: LVT_resetMetric_AnomalyRMSE
! \label{LVT_resetMetric_AnomalyRMSE}
!
! !INTERFACE: 
  subroutine LVT_resetMetric_AnomalyRMSE(alarm)
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
    integer                 :: i,k,l
    type(LVT_metadataEntry), pointer :: model
    type(LVT_statsEntry)   , pointer :: stats

    call LVT_getDataStream1Ptr(model)
    call LVT_getstatsEntryPtr(stats)

    do while(associated(model))
       if(stats%selectOpt.eq.1) then 
          do k=1,model%selectNlevs
             if(LVT_metrics%armse%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   stats%armse%value_ts(:,k,l) = 0.0
                   stats%armse%count_value_ts(:,k,l)=0 
                   if(alarm) then
                      stats%armse%tavg_value_ts(:,k,l) = 0.0
                      stats%armse%tavg_count_value_ts(:,k,l)=0 
                   endif
                enddo
             endif
          enddo
       endif
       
       model => model%next
       stats => stats%next
    enddo

  end subroutine LVT_resetMetric_AnomalyRMSE


!BOP
! 
! !ROUTINE: LVT_writerestart_AnomalyRMSE
! 
! !INTERFACE:
  subroutine LVT_writerestart_AnomalyRMSE(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for AnomalyRMSE metric computations
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

    if(LVT_metrics%armse%selectOpt.eq.1) then 
       call LVT_getDataStream1Ptr(model)
       call LVT_getDataStream2Ptr(obs)
       call LVT_getstatsEntryPtr(stats)
       
       do while(associated(model))
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do k=1,model%selectNlevs
                do m = 1, 12
                   do l=1,LVT_rc%strat_nlevels  
                      
                      call LVT_writevar_restart(ftn,&
                           stats%armse%model_value_climo(:,k,m,l))
                      call LVT_writevar_restart(ftn,&
                           stats%armse%obs_value_climo(:,k,m,l))
                      call LVT_writevar_restart(ftn,&
                           stats%armse%count_model_value_climo(:,k,m,l))
                      call LVT_writevar_restart(ftn,&
                           stats%armse%count_obs_value_climo(:,k,m,l))                    
                   enddo
                enddo
             enddo
             if(pass.eq.2) then
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_restart(ftn,&
                           stats%armse%value_total(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%armse%count_value_total(:,k,l))
                   enddo
                enddo
                if(LVT_metrics%armse%computeSC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nasc
                         call LVT_writevar_restart(ftn,&
                              stats%armse%value_asc(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%armse%count_value_asc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_metrics%armse%computeADC.eq.1) then
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nadc 
                         call LVT_writevar_restart(ftn,&
                              stats%armse%value_adc(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%armse%count_value_adc(:,k,l))
                      enddo
                   enddo
                endif

             endif
          endif
          
          model => model%next
          obs => obs%next
          stats => stats%next
       enddo
    end if
    
  end subroutine LVT_writerestart_AnomalyRMSE

!BOP
! 
! !ROUTINE: LVT_readrestart_AnomalyRMSE
! 
! !INTERFACE:
  subroutine LVT_readrestart_AnomalyRMSE(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for AnomalyRMSE metric computations
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

    if(LVT_metrics%armse%selectOpt.eq.1) then 
       
       call LVT_getDataStream1Ptr(model)
       call LVT_getDataStream2Ptr(obs)
       call LVT_getstatsEntryPtr(stats)
       
       do while(associated(model))
 
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do k=1,model%selectNlevs
                do m = 1, 12
                   do l=1,LVT_rc%strat_nlevels  
                      
                      call LVT_readvar_restart(ftn,&
                           stats%armse%model_value_climo(:,k,m,l))
                      call LVT_readvar_restart(ftn,&
                           stats%armse%obs_value_climo(:,k,m,l))
                      call LVT_readvar_restart(ftn,&
                           stats%armse%count_model_value_climo(:,k,m,l))
                      call LVT_readvar_restart(ftn,&
                           stats%armse%count_obs_value_climo(:,k,m,l))                    
                   enddo
                enddo
             enddo
             if(LVT_rc%curr_pass.eq.2) then
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_readvar_restart(ftn,&
                           stats%armse%value_total(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%armse%count_value_total(:,k,l))
                   enddo
                enddo
                if(LVT_metrics%armse%computeSC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nasc
                         call LVT_readvar_restart(ftn,&
                              stats%armse%value_asc(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%armse%count_value_asc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_metrics%armse%computeADC.eq.1) then
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nadc 
                         call LVT_readvar_restart(ftn,&
                              stats%armse%value_adc(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%armse%count_value_adc(:,k,l))
                      enddo
                   enddo
                endif
                
             endif
          endif

          model => model%next
          obs => obs%next
          stats => stats%next
       enddo
    end if
  end subroutine LVT_readrestart_AnomalyRMSE


end module LVT_AnomalyRMSEMod
