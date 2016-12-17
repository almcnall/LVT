!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_RawCorrMod
! \label(LVT_RawCorrMod)
!
! !INTERFACE:
module LVT_RawCorrMod
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
!   This module handles the raw correlation (pearson correlation
!   coefficient)  computations by comparing 
!   the LIS output to the specified observations. 
!
! Time series, seasonal cycles, average diurnal cycle options are 
! not supported for Raw Correlation computation
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
  public :: LVT_initRawCorr
  public :: LVT_diagnoseRawCorr
  public :: LVT_computeRawCorr 
  public :: LVT_writeMetric_RawCorr
  public :: LVT_resetMetric_RawCorr
  public :: LVT_writerestart_RawCorr
  public :: LVT_readrestart_RawCorr

contains
  
!BOP
! 
! !ROUTINE: LVT_initRawCorr
! \label{LVT_initRawCorr}
!
! !INTERFACE: 
  subroutine LVT_initRawCorr(model,obs,stats,metric)
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
       allocate(stats%rcorr%sxy_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
       allocate(stats%rcorr%sxx_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
       allocate(stats%rcorr%syy_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
       allocate(stats%rcorr%sx_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
       allocate(stats%rcorr%sy_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
       allocate(stats%rcorr%rval_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
       allocate(stats%rcorr%count_value(LVT_rc%ngrid, model%selectNlevs, LVT_rc%strat_nlevels))       
       stats%rcorr%sxy_r = 0 
       stats%rcorr%sxx_r = 0 
       stats%rcorr%syy_r = 0 
       stats%rcorr%sx_r = 0 
       stats%rcorr%sy_r = 0 
       stats%rcorr%rval_r = 0
       stats%rcorr%count_value = 0 

       allocate(stats%rcorr%value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%rcorr%value_ci = LVT_rc%udef

       if(metric%timeOpt.eq.1) then 
          allocate(stats%rcorr%sxy_ts_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%rcorr%sxx_ts_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%rcorr%syy_ts_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%rcorr%sx_ts_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%rcorr%sy_ts_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%rcorr%rval_ts_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%rcorr%count_value_ts(LVT_rc%ngrid, model%selectNlevs, LVT_rc%strat_nlevels))       
          stats%rcorr%sxy_ts_r = 0 
          stats%rcorr%sxx_ts_r = 0 
          stats%rcorr%syy_ts_r = 0 
          stats%rcorr%sx_ts_r = 0 
          stats%rcorr%sy_ts_r = 0 
          stats%rcorr%rval_ts_r = 0
          stats%rcorr%count_value_ts = 0 

          allocate(stats%rcorr%tavg_value_ts(LVT_rc%ngrid, &
               model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%rcorr%tavg_count_value_ts(LVT_rc%ngrid, &
               model%selectNlevs,LVT_rc%strat_nlevels))
          stats%rcorr%tavg_value_ts = 0.0
          stats%rcorr%tavg_count_value_ts=0 

          if(metric%computeSC.eq.1) then 
             allocate(stats%rcorr%value_asc(LVT_rc%ngrid, model%selectNlevs, LVT_rc%nasc))
             stats%rcorr%value_asc = 0.0
             allocate(stats%rcorr%count_value_asc(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%nasc))
             stats%rcorr%count_value_asc = 0
          endif

          if(metric%computeADC.eq.1) then 
             allocate(stats%rcorr%value_adc(LVT_rc%ngrid, model%selectNlevs, LVT_rc%nadc))
             stats%rcorr%value_adc = 0.0
             allocate(stats%rcorr%count_value_adc(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%nadc))
             stats%rcorr%count_value_adc = 0
          endif
          
       endif

    endif
    
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1
    metric%obsData = .false. 
    metric%stdevFlag = .false. 

  end subroutine LVT_initRawCorr
  
!BOP
! 
! !ROUTINE: LVT_diagnoseRawCorr
! \label{LVT_diagnoseRawCorr}
!
! !INTERFACE:  
  subroutine LVT_diagnoseRawCorr(pass)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the RawCorr calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleRawCorr](\ref{diagnoseSingleRawCorr})
!     updates the RawCorr computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer                 :: pass
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%rcorr%selectOpt.eq.1.or.&
            LVT_metrics%rcorr%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleRawCorr(obs,model,stats, &
                  LVT_metrics%rcorr)
             
             model => model%next
             obs   => obs%next
             stats => stats%next

          enddo
       endif
    endif
  end subroutine LVT_diagnoseRawCorr

!BOP
! 
! !ROUTINE: diagnoseSingleRawCorr
! \label{diagnoseSingleRawCorr}
!
! !INTERFACE: 
  subroutine diagnoseSingleRawCorr(obs, model, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the RawCorr computation (updates the running 
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
                   stats%rcorr%sxy_r(t,k,1) = stats%rcorr%sxy_r(t,k,1) + &
                        (model%value(t,k))*&
                        (obs%value(t,k))
                   stats%rcorr%sx_r(t,k,1) = stats%rcorr%sx_r(t,k,1) +&
                        (model%value(t,k))
                   stats%rcorr%sy_r(t,k,1) = stats%rcorr%sy_r(t,k,1) +&
                        (obs%value(t,k))
                   stats%rcorr%sxx_r(t,k,1) = stats%rcorr%sxx_r(t,k,1) +&
                        (model%value(t,k))*&
                        (model%value(t,k))
                   stats%rcorr%syy_r(t,k,1) = stats%rcorr%syy_r(t,k,1) +&
                        (obs%value(t,k))*&
                        (obs%value(t,k))               
                   stats%rcorr%count_value(t,k,1) = stats%rcorr%count_value(t,k,1) + 1
                   if(LVT_rc%strat_nlevels.gt.1) then 
                      if(LVT_stats%strat_var(t,k).gt.&
                           LVT_rc%strat_var_threshold) then 
                         stats%rcorr%sxy_r(t,k,2) = stats%rcorr%sxy_r(t,k,2) + &
                              (model%value(t,k))*&
                              (obs%value(t,k))
                         stats%rcorr%sx_r(t,k,2) = stats%rcorr%sx_r(t,k,2) +&
                              (model%value(t,k))
                         stats%rcorr%sy_r(t,k,2) = stats%rcorr%sy_r(t,k,2) +&
                              (obs%value(t,k))
                         stats%rcorr%sxx_r(t,k,2) = stats%rcorr%sxx_r(t,k,2) +&
                              (model%value(t,k))*&
                              (model%value(t,k))
                         stats%rcorr%syy_r(t,k,2) = stats%rcorr%syy_r(t,k,2) +&
                              (obs%value(t,k))*&
                              (obs%value(t,k))               
                         stats%rcorr%count_value(t,k,2) = stats%rcorr%count_value(t,k,2) + 1
                         
                      elseif(LVT_stats%strat_var(t,k).le.&
                           LVT_rc%strat_var_threshold) then 
                         stats%rcorr%sxy_r(t,k,3) = stats%rcorr%sxy_r(t,k,3) + &
                              (model%value(t,k))*&
                              (obs%value(t,k))
                         stats%rcorr%sx_r(t,k,3) = stats%rcorr%sx_r(t,k,3) +&
                              (model%value(t,k))
                         stats%rcorr%sy_r(t,k,3) = stats%rcorr%sy_r(t,k,3) +&
                              (obs%value(t,k))
                         stats%rcorr%sxx_r(t,k,3) = stats%rcorr%sxx_r(t,k,3) +&
                              (model%value(t,k))*&
                              (model%value(t,k))
                         stats%rcorr%syy_r(t,k,3) = stats%rcorr%syy_r(t,k,3) +&
                              (obs%value(t,k))*&
                              (obs%value(t,k))               
                         stats%rcorr%count_value(t,k,3) = stats%rcorr%count_value(t,k,3) + 1
                         
                      endif
                   endif
                   if(metric%timeOpt.eq.1) then 
                      stats%rcorr%sxy_ts_r(t,k,1) = stats%rcorr%sxy_ts_r(t,k,1) + &
                           (model%value(t,k))*&
                           (obs%value(t,k))
                      stats%rcorr%sx_ts_r(t,k,1) = stats%rcorr%sx_ts_r(t,k,1) +&
                           (model%value(t,k))
                      stats%rcorr%sy_ts_r(t,k,1) = stats%rcorr%sy_ts_r(t,k,1) +&
                           (obs%value(t,k))
                      stats%rcorr%sxx_ts_r(t,k,1) = stats%rcorr%sxx_ts_r(t,k,1) +&
                           (model%value(t,k))*&
                           (model%value(t,k))
                      stats%rcorr%syy_ts_r(t,k,1) = stats%rcorr%syy_ts_r(t,k,1) +&
                           (obs%value(t,k))*&
                           (obs%value(t,k))     
                      stats%rcorr%count_value_ts(t,k,1) = stats%rcorr%count_value_ts(t,k,1) + 1
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then 
                            stats%rcorr%sxy_ts_r(t,k,2) = stats%rcorr%sxy_ts_r(t,k,2) + &
                                 (model%value(t,k))*&
                                 (obs%value(t,k))
                            stats%rcorr%sx_ts_r(t,k,2) = stats%rcorr%sx_ts_r(t,k,2) +&
                                 (model%value(t,k))
                            stats%rcorr%sy_ts_r(t,k,2) = stats%rcorr%sy_ts_r(t,k,2) +&
                                 (obs%value(t,k))
                            stats%rcorr%sxx_ts_r(t,k,2) = stats%rcorr%sxx_ts_r(t,k,2) +&
                                 (model%value(t,k))*&
                                 (model%value(t,k))
                            stats%rcorr%syy_ts_r(t,k,2) = stats%rcorr%syy_ts_r(t,k,2) +&
                                 (obs%value(t,k))*&
                                 (obs%value(t,k))               
                            stats%rcorr%count_value_ts(t,k,2) = &
                                 stats%rcorr%count_value_ts(t,k,2) + 1
                            
                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then 
                            stats%rcorr%sxy_ts_r(t,k,3) = stats%rcorr%sxy_ts_r(t,k,3) + &
                                 (model%value(t,k))*&
                                 (obs%value(t,k))
                            stats%rcorr%sx_ts_r(t,k,3) = stats%rcorr%sx_ts_r(t,k,3) +&
                                 (model%value(t,k))
                            stats%rcorr%sy_ts_r(t,k,3) = stats%rcorr%sy_ts_r(t,k,3) +&
                                 (obs%value(t,k))
                            stats%rcorr%sxx_ts_r(t,k,3) = stats%rcorr%sxx_ts_r(t,k,3) +&
                                 (model%value(t,k))*&
                                 (model%value(t,k))
                            stats%rcorr%syy_ts_r(t,k,3) = stats%rcorr%syy_ts_r(t,k,3) +&
                                 (obs%value(t,k))*&
                                 (obs%value(t,k))               
                            stats%rcorr%count_value_ts(t,k,3) = &
                                 stats%rcorr%count_value_ts(t,k,3) + 1
                            
                         endif
                      endif
                   endif
                endif
             endif
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleRawCorr


!BOP
! 
! !ROUTINE: LVT_computeRawCorr
! \label{LVT_computeRawCorr}
!
! !INTERFACE: 
  subroutine LVT_computeRawCorr(pass,alarm)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute RawCorr values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleRawCorr](\ref{computeSingleRawCorr})
!     updates the RawCorr computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     RawCorr computation has been reached
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: pass
    logical               :: alarm
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    integer               :: i

    if(pass.eq.1) then 

       if(LVT_metrics%rcorr%selectOpt.eq.1.or.&
            LVT_metrics%rcorr%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%rcorr%timeOpt.eq.1.and.&
                  LVT_metrics%rcorr%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%rcorr%ftn_ts_loc(i),200,advance='no') &
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
             call computeSingleRawCorr(alarm,obs,model,stats,&
                  LVT_metrics%rcorr)

             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo

          if(alarm) then 
             if(LVT_metrics%rcorr%timeOpt.eq.1.and.&
                  LVT_metrics%rcorr%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%rcorr%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif

  end subroutine LVT_ComputeRawCorr

!BOP
! 
! !ROUTINE: computeSingleRawCorr
! \label{computeSingleRawCorr}
!
! !INTERFACE: 
  subroutine computeSingleRawCorr(alarm,obs, model,stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the RawCorr values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     RawCorr computation has been reached
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

    logical                 :: alarm
    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer  :: t,l,k,tind
    real     :: numer
    real     :: denom_sq
    real     :: denom


    if(metric%timeOpt.eq.1) then 
       if(alarm) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%rcorr%count_value_ts(t,k,l).ne.0.and.&
                        stats%rcorr%count_value_ts(t,k,l) &
                        .gt.LVT_rc%obsCountThreshold) then 
                      
                      numer = (float(stats%rcorr%count_value_ts(t,k,l))* &
                              stats%rcorr%sxy_ts_r(t,k,l) - &
                              stats%rcorr%sx_ts_r(t,k,l)*stats%rcorr%sy_ts_r(t,k,l))
                      denom_sq =  (float(stats%rcorr%count_value_ts(t,k,l))* &
                           stats%rcorr%sxx_ts_r(t,k,l)-&
                           stats%rcorr%sx_ts_r(t,k,l)**2)* &
                           (float(stats%rcorr%count_value_ts(t,k,l))*&
                           stats%rcorr%syy_ts_r(t,k,l)-&
                           stats%rcorr%sy_ts_r(t,k,l)**2)
                      if(denom_sq.gt.0) then 
                         denom = sqrt(denom_sq)
                      else
                         denom = 0.0
                      endif
                      
                      if(denom.ne.0) then 
                         stats%rcorr%rval_ts_r(t,k,l) = numer/denom
                         stats%rcorr%tavg_value_ts(t,k,l) = & 
                              stats%rcorr%tavg_value_ts(t,k,l) + & 
                              stats%rcorr%rval_ts_r(t,k,l)
                         stats%rcorr%tavg_count_value_ts(t,k,l) = & 
                              stats%rcorr%tavg_count_value_ts(t,k,l) + 1
                      else
                         stats%rcorr%rval_ts_r(t,k,l) = LVT_rc%udef
                      endif
                      
                   endif
                enddo

                if(metric%computeSC.eq.1) then 
                   if(stats%rcorr%rval_ts_r(t,k,l).ne.LVT_rc%udef) then 
                      call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                           tind)
                      stats%rcorr%value_asc(t,k,tind) = & 
                           stats%rcorr%value_asc(t,k,tind) + & 
                           stats%rcorr%rval_ts_r(t,k,l)
                      stats%rcorr%count_value_asc(t,k,tind) = & 
                           stats%rcorr%count_value_asc(t,k,tind) + 1
                   endif
                endif
                
                if(metric%computeADC.eq.1) then 
                   if(stats%rcorr%rval_ts_r(t,k,l).ne.LVT_rc%udef) then 
                      call LVT_getADCTimeIndex(tind)
                      stats%rcorr%value_adc(t,k,tind) = & 
                           stats%rcorr%value_adc(t,k,tind) + & 
                           stats%rcorr%rval_ts_r(t,k,l)
                      stats%rcorr%count_value_adc(t,k,tind) = & 
                           stats%rcorr%count_value_adc(t,k,tind) + 1
                   endif
                endif

             enddo
          enddo

          if(metric%extractTS.eq.1) then 
             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%rcorr%tavg_count_value_ts(t,k,l).gt.0) then 
                         stats%rcorr%tavg_value_ts(t,k,l) = &
                              stats%rcorr%tavg_value_ts(t,k,l)/&
                              stats%rcorr%tavg_count_value_ts(t,k,l)
                      else
                         stats%rcorr%tavg_value_ts(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                enddo
             enddo

             call LVT_writeTSinfo(metric%ftn_ts_loc,&
                  model,&
                  LVT_rc%ngrid,&
                  stats%rcorr%tavg_value_ts,&
                  stats%rcorr%tavg_count_value_ts)
          endif
          
       endif
    endif
    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%rcorr%count_value(t,k,l).ne.0.and.&
                        stats%rcorr%count_value(t,k,l) &
                        .gt.LVT_rc%obsCountThreshold) then 

                      numer = (float(stats%rcorr%count_value(t,k,l))* &
                           stats%rcorr%sxy_r(t,k,l) - &
                           stats%rcorr%sx_r(t,k,l)*stats%rcorr%sy_r(t,k,l))
                      denom_sq =  (float(stats%rcorr%count_value(t,k,l))* &
                           stats%rcorr%sxx_r(t,k,l)-&
                           stats%rcorr%sx_r(t,k,l)**2)* &
                           (float(stats%rcorr%count_value(t,k,l))*&
                           stats%rcorr%syy_r(t,k,l)-&
                           stats%rcorr%sy_r(t,k,l)**2)
                      if(denom_sq.gt.0) then 
                         denom = sqrt(denom_sq)
                      else
                         denom = 0.0
                      endif

                      if(denom.ne.0) then 
                         stats%rcorr%rval_r(t,k,l) = numer/denom
                      else
                         stats%rcorr%rval_r(t,k,l) = LVT_rc%udef
                         stats%rcorr%count_value(t,k,l) = LVT_rc%udef
                      endif
                   else
                      stats%rcorr%rval_r(t,k,l) = LVT_rc%udef
                      stats%rcorr%count_value(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1,LVT_rc%nasc
                      if(stats%rcorr%count_value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then 
                         stats%rcorr%value_asc(t,k,l) = &
                              stats%rcorr%value_asc(t,k,l)/&
                              stats%rcorr%count_value_asc(t,k,l) 
                      else
                         stats%rcorr%value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
                if(metric%computeADC.eq.1) then 
                   do l=1,LVT_rc%nadc
                      if(stats%rcorr%count_value_adc(t,k,l).gt.&
                           LVT_rc%ADCCountThreshold) then 
                         stats%rcorr%value_adc(t,k,l) = &
                              stats%rcorr%value_adc(t,k,l)/&
                              stats%rcorr%count_value_adc(t,k,l) 
                      else
                         stats%rcorr%value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif


             enddo
          enddo
          
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%rcorr%rval_r(:,k,l),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%rcorr%value_ci(k,l))
             enddo
          enddo
       
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid,stats%rcorr%rval_r)      

          if(metric%computeSC.eq.1) then 
             call LVT_writeSeasonalCycleInfo(model,obs,stats, metric, &
                  LVT_rc%ngrid, stats%rcorr%value_asc, stats%rcorr%count_value_asc)
          endif
          if(metric%computeADC.eq.1) then 

             call LVT_writeAvgDiurnalCycleInfo(model, obs, stats, metric, & 
                  LVT_rc%ngrid, stats%rcorr%value_adc, stats%rcorr%count_value_adc)
          endif
       endif
    endif

  end subroutine computeSingleRawCorr

!BOP
! 
! !ROUTINE: LVT_writeMetric_RawCorr
! \label{LVT_writeMetric_RawCorr}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_RawCorr(pass,final,vlevels,stats,obs)
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

    if(pass.eq.LVT_metrics%rcorr%npass) then 
       if(final.ne.1) then 
          if(stats%selectOpt.eq.1) then 
             do k=1,vlevels
                do l=1,LVT_rc%strat_nlevels
                   
                   call LVT_writevar_gridded(LVT_metrics%rcorr%ftn_ts, &
                        stats%rcorr%tavg_value_ts(:,k,l),stats%vid_ts(LVT_Rcorrid,1))
                   call LVT_writevar_gridded(LVT_metrics%rcorr%ftn_ts, &
                        real(stats%rcorr%tavg_count_value_ts(:,k,l)),&
                        stats%vid_count_ts(LVT_Rcorrid,1))                
                enddo
             enddo
          endif
       else
          if(LVT_metrics%rcorr%selectOpt.eq.1) then
             if(stats%selectOpt.eq.1) then 
                do k=1,vlevels
                   do l=1,LVT_rc%strat_nlevels
                      
                      call LVT_writevar_gridded(LVT_metrics%rcorr%ftn_total, &
                           stats%rcorr%rval_r(:,k,l),stats%vid_total(LVT_Rcorrid,1))
                      call LVT_writevar_gridded(LVT_metrics%rcorr%ftn_total, &
                           real(stats%rcorr%count_value(:,k,l)),&
                           stats%vid_count_total(LVT_Rcorrid,1))                
                   enddo

                   if(LVT_metrics%rcorr%computeSC.eq.1) then 
                      do tind = 1, LVT_rc%nasc
                         call LVT_writevar_gridded(&
                              LVT_metrics%rcorr%ftn_total,&
                              stats%rcorr%value_asc(:,k,tind), &
                              stats%vid_sc_total(tind,LVT_Rcorrid,1),k)
                      enddo
                   endif
                   
                   if(LVT_metrics%rcorr%computeADC.eq.1) then 
                      do tind = 1,LVT_rc%nadc
                         call LVT_writevar_gridded(LVT_metrics%rcorr%ftn_total,&
                              stats%rcorr%value_adc(:,k,tind),&
                              stats%vid_adc_total(tind,LVT_Rcorrid,1),k)
                      enddo
                   endif

                   call LVT_writeSummaryStats(&
                        LVT_metrics%rcorr%ftn_summ,&
                        LVT_metrics%rcorr%short_name,&
                        LVT_rc%ngrid,&
                        stats%rcorr%rval_r(:,k,:), &
                        stats%rcorr%count_value(:,k,:),stats%standard_name,&
                        stats%rcorr%value_ci(k,:))
                enddo
             endif
          endif
       endif
    endif

  end subroutine LVT_writeMetric_RawCorr


!BOP
! 
! !ROUTINE: LVT_resetMetric_RawCorr
! \label(LVT_resetMetric_RawCorr)
!
! !INTERFACE:
  subroutine LVT_resetMetric_RawCorr(alarm)

! !INPUT PARAMETERS: 
    logical           :: alarm
 

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
             if(LVT_metrics%rcorr%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels

                   if(alarm) then 
                      stats%rcorr%sxy_ts_r(:,k,l) = 0.0
                      stats%rcorr%sxx_ts_r(:,k,l) = 0.0
                      stats%rcorr%syy_ts_r(:,k,l) = 0.0
                      stats%rcorr%sx_ts_r(:,k,l) = 0.0
                      stats%rcorr%sy_ts_r(:,k,l) = 0.0
                      stats%rcorr%rval_ts_r(:,k,l) = 0.0
                      stats%rcorr%count_value_ts(:,k,l)=0 
                      
                      stats%rcorr%tavg_value_ts(:,k,l) = 0.0
                      stats%rcorr%tavg_count_value_ts(:,k,l)=0 
                   endif

                enddo
             endif
             
          enddo
          
       endif
       model => model%next
       stats => stats%next

    enddo
  end subroutine LVT_resetMetric_RawCorr


!BOP
! 
! !ROUTINE: LVT_writerestart_RawCorr
! 
! !INTERFACE:
  subroutine LVT_writerestart_RawCorr(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for RawCorr metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    integer              :: k,l
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(LVT_metrics%rcorr%selectOpt.eq.1) then 
       
       call LVT_getDataStream1Ptr(model)
       call LVT_getDataStream2Ptr(obs)
       call LVT_getstatsEntryPtr(stats)

       do while(associated(model))
          
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels     
                   call LVT_writevar_restart(ftn,stats%rcorr%sxy_r(:,k,l))
                   call LVT_writevar_restart(ftn,stats%rcorr%sxx_r(:,k,l))
                   call LVT_writevar_restart(ftn,stats%rcorr%syy_r(:,k,l))
                   call LVT_writevar_restart(ftn,stats%rcorr%sx_r(:,k,l))
                   call LVT_writevar_restart(ftn,stats%rcorr%sy_r(:,k,l))
                   call LVT_writevar_restart(ftn,stats%rcorr%rval_r(:,k,l))
                   call LVT_writevar_restart(ftn,stats%rcorr%count_value(:,k,l))
                enddo
             enddo

             if(LVT_metrics%rcorr%computeSC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc
                      call LVT_writevar_restart(ftn,&
                           stats%rcorr%value_asc(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%rcorr%count_value_asc(:,k,l))
                   enddo
                enddo
             endif
             if(LVT_metrics%rcorr%computeADC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nadc
                      call LVT_writevar_restart(ftn,&
                           stats%rcorr%value_adc(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%rcorr%count_value_adc(:,k,l))
                   enddo
                enddo
             endif
             
          endif

          model => model%next
          obs   => obs%next
          stats => stats%next
          
       enddo
    end if
    
  end subroutine LVT_writerestart_RawCorr

!BOP
! 
! !ROUTINE: LVT_readrestart_RawCorr
! 
! !INTERFACE:
  subroutine LVT_readrestart_RawCorr(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for RawCorr metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats
    integer              :: k,l

    if(LVT_metrics%rcorr%selectOpt.eq.1) then 
       call LVT_getDataStream1Ptr(model)
       call LVT_getDataStream2Ptr(obs)
       call LVT_getstatsEntryPtr(stats)
       
       do while(associated(model))
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   call LVT_readvar_restart(ftn,stats%rcorr%sxy_r(:,k,l))
                   call LVT_readvar_restart(ftn,stats%rcorr%sxx_r(:,k,l))
                   call LVT_readvar_restart(ftn,stats%rcorr%syy_r(:,k,l))
                   call LVT_readvar_restart(ftn,stats%rcorr%sx_r(:,k,l))
                   call LVT_readvar_restart(ftn,stats%rcorr%sy_r(:,k,l))
                   call LVT_readvar_restart(ftn,stats%rcorr%rval_r(:,k,l))
                   call LVT_readvar_restart(ftn,stats%rcorr%count_value(:,k,l))
                enddo
             enddo

             if(LVT_metrics%rcorr%computeSC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc
                      call LVT_readvar_restart(ftn,&
                           stats%rcorr%value_asc(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%rcorr%count_value_asc(:,k,l))
                   enddo
                enddo
             endif
             if(LVT_metrics%rcorr%computeADC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nadc
                      call LVT_readvar_restart(ftn,&
                           stats%rcorr%value_adc(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%rcorr%count_value_adc(:,k,l))
                   enddo
                enddo
             endif
             
          endif
          
          model => model%next
          obs   => obs%next
          stats => stats%next

       enddo
    end if

  end subroutine LVT_readrestart_RawCorr
end module LVT_RawCorrMod
