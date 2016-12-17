!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_SumMod
! \label(LVT_SumMod)
!
! !INTERFACE:
module LVT_SumMod
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
!  This module handles the computations required to compute sum values
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
  public :: LVT_initSum
  public :: LVT_diagnoseSum
  public :: LVT_computeSum
  public :: LVT_writeMetric_Sum
  public :: LVT_resetMetric_Sum
  public :: LVT_writerestart_Sum
  public :: LVT_readrestart_Sum
!EOP
  
  private

contains
  subroutine LVT_initSum(model, obs, stats,metric)
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
       allocate(stats%sum%model_value_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%sum%model_value_total = 0.0
       allocate(stats%sum%count_model_value_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%sum%count_model_value_total = 0
       allocate(stats%sum%model_value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%sum%model_value_ci = LVT_rc%udef
       
       if(metric%computeSC.eq.1) then 
          allocate(stats%sum%model_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%sum%model_value_asc = 0.0
          allocate(stats%sum%count_model_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%sum%count_model_value_asc = 0
       endif
       if(metric%computeADC.eq.1) then 
          allocate(stats%sum%model_value_adc(LVT_rc%ngrid, obs%selectNlevs,&
               LVT_rc%nadc))
          stats%sum%model_value_adc = 0.0
          allocate(stats%sum%count_model_value_adc(LVT_rc%ngrid, obs%selectNlevs,&
               LVT_rc%nadc))
          stats%sum%count_model_value_adc = 0
       endif
       
       if(LVT_rc%obssource(2).ne."none") then 
          if(obs%selectNlevs.ge.1) then 
             allocate(stats%sum%obs_value_total(LVT_rc%ngrid, obs%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%sum%obs_value_total = 0.0
             allocate(stats%sum%count_obs_value_total(LVT_rc%ngrid, obs%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%sum%count_obs_value_total = 0
             allocate(stats%sum%obs_value_ci(obs%selectNlevs,LVT_rc%strat_nlevels))
             stats%sum%obs_value_ci = LVT_rc%udef
             
             if(metric%computeSC.eq.1) then 
                allocate(stats%sum%obs_value_asc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nasc))
                stats%sum%obs_value_asc = 0.0
                allocate(stats%sum%count_obs_value_asc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nasc))
                stats%sum%count_obs_value_asc = 0
             endif
             if(metric%computeADC.eq.1) then 
                allocate(stats%sum%obs_value_adc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nadc))
                stats%sum%obs_value_adc = 0.0
                allocate(stats%sum%count_obs_value_adc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nadc))
                stats%sum%count_obs_value_adc = 0
             endif
          endif
       endif
    endif
    if(metric%timeOpt.eq.1) then 
       allocate(stats%sum%model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%sum%model_value_ts = 0.0
       allocate(stats%sum%count_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%sum%count_model_value_ts = 0

       allocate(stats%sum%tavg_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%sum%tavg_model_value_ts = 0.0
       allocate(stats%sum%tavg_count_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%sum%tavg_count_model_value_ts = 0
       
       if(LVT_rc%obssource(2).ne."none") then 
          if(obs%selectNlevs.ge.1) then 
             allocate(stats%sum%obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%sum%obs_value_ts = 0.0
             allocate(stats%sum%count_obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%sum%count_obs_value_ts = 0

             allocate(stats%sum%tavg_obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%sum%tavg_obs_value_ts = 0.0
             allocate(stats%sum%tavg_count_obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%sum%tavg_count_obs_value_ts = 0

          endif
       endif
    endif

!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1    
    if(LVT_rc%obssource(2).ne."none") then 
       metric%obsData = .true. 
    else
       metric%obsData = .false. 
    endif
    metric%stdevFlag = .false. 

  end subroutine LVT_initSum

!BOP
! 
! !ROUTINE: LVT_diagnoseSum
! \label{LVT_diagnoseSum}
!
! !INTERFACE: 
  subroutine LVT_diagnoseSum(pass)
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
!   calculating the sum of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelSum](\ref{diagnoseSingleModelSum})
!     updates the sum computation for a single variable 
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
       if(LVT_metrics%sum%selectOpt.eq.1.or.&
            LVT_metrics%sum%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleSum(model,obs,stats,&
                  LVT_metrics%sum)
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseSum

!BOP
! 
! !ROUTINE: diagnoseSingleSum
! \label{diagnoseSingleSum}
!
! !INTERFACE: 
  subroutine diagnoseSingleSum(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the sum computation of the 
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
    integer    :: t,k,tind

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then        
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             if(model%count(t,k).gt.0) then
                if(metric%selectOpt.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      stats%sum%model_value_total(t,k,1) = &
                           stats%sum%model_value_total(t,k,1) + &
                           model%value(t,k)
                      stats%sum%count_model_value_total(t,k,1) = &
                           stats%sum%count_model_value_total(t,k,1) + 1
                   endif
                endif
                
                if(metric%timeOpt.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      stats%sum%model_value_ts(t,k,1) = stats%sum%model_value_ts(t,k,1)+&
                           model%value(t,k)
                      stats%sum%count_model_value_ts(t,k,1) = & 
                           stats%sum%count_model_value_ts(t,k,1)+1
                   endif
                endif
                if(metric%computeSC.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                           tind)
                      stats%sum%model_value_asc(t,k,tind) = &
                           stats%sum%model_value_asc(t,k,tind)+&
                           model%value(t,k)
                      stats%sum%count_model_value_asc(t,k,tind) = &
                           stats%sum%count_model_value_asc(t,k,tind) + 1
                   endif
                endif
                if(metric%computeADC.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      call LVT_getADCTimeIndex(tind)
                      stats%sum%model_value_adc(t,k,tind) = &
                           stats%sum%model_value_adc(t,k,tind)+&
                           model%value(t,k)
                      stats%sum%count_model_value_adc(t,k,tind) = &
                           stats%sum%count_model_value_adc(t,k,tind) + 1
                   endif
                endif
                if(LVT_rc%strat_nlevels.gt.1) then 
                   if(LVT_stats%strat_var(t,k).gt.&
                        LVT_rc%strat_var_threshold) then
                      if(metric%selectOpt.eq.1) then
                         if(model%value(t,k).ne.LVT_rc%udef) then  
                            stats%sum%model_value_total(t,k,2) = &
                                 stats%sum%model_value_total(t,k,2) + &
                                 model%value(t,k)
                            stats%sum%count_model_value_total(t,k,2) = &
                                 stats%sum%count_model_value_total(t,k,2) + 1
                         endif
                      endif
                      if(metric%timeOpt.eq.1) then 
                         if(model%value(t,k).ne.LVT_rc%udef) then 
                            stats%sum%model_value_ts(t,k,2) = &
                                 stats%sum%model_value_ts(t,k,2)+&
                                 model%value(t,k)
                            stats%sum%count_model_value_ts(t,k,2) = & 
                                 stats%sum%count_model_value_ts(t,k,2)+1
                         endif
                      endif
                   elseif(LVT_stats%strat_var(t,k).le.&
                        LVT_rc%strat_var_threshold) then
                      if(metric%selectOpt.eq.1) then 
                         if(model%value(t,k).ne.LVT_rc%udef) then 
                            stats%sum%model_value_total(t,k,3) = &
                                 stats%sum%model_value_total(t,k,3) + &
                                 model%value(t,k)
                            stats%sum%count_model_value_total(t,k,3) = &
                                 stats%sum%count_model_value_total(t,k,3) + 1
                         endif
                      endif
                      if(metric%timeOpt.eq.1) then 
                         if(model%value(t,k).ne.LVT_rc%udef) then 
                            stats%sum%model_value_ts(t,k,3) = &
                                 stats%sum%model_value_ts(t,k,3)+&
                                 model%value(t,k)
                            stats%sum%count_model_value_ts(t,k,3) = & 
                                 stats%sum%count_model_value_ts(t,k,3)+1
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
                         stats%sum%obs_value_total(t,k,1) = &
                              stats%sum%obs_value_total(t,k,1) + &
                              obs%value(t,k)
                         stats%sum%count_obs_value_total(t,k,1) = &
                              stats%sum%count_obs_value_total(t,k,1) + 1
                      endif
                      
                      if(metric%timeOpt.eq.1.and.obs%value(t,k).ne.LVT_rc%udef) then 
                         stats%sum%obs_value_ts(t,k,1) = stats%sum%obs_value_ts(t,k,1)+&
                              obs%value(t,k)
                         stats%sum%count_obs_value_ts(t,k,1) = & 
                              stats%sum%count_obs_value_ts(t,k,1)+1
                      endif
                      if(metric%computeSC.eq.1.and.obs%value(t,k).ne.LVT_rc%udef) then 
                         call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,tind)
                         stats%sum%obs_value_asc(t,k,tind) = &
                              stats%sum%obs_value_asc(t,k,tind)+&
                              obs%value(t,k)
                         stats%sum%count_obs_value_asc(t,k,tind) = &
                              stats%sum%count_obs_value_asc(t,k,tind) + 1
                      endif
                      if(metric%computeADC.eq.1.and.obs%value(t,k).ne.LVT_rc%udef) then 
                         call LVT_getADCTimeIndex(tind)
                         stats%sum%obs_value_adc(t,k,tind) = &
                              stats%sum%obs_value_adc(t,k,tind)+&
                              obs%value(t,k)
                         stats%sum%count_obs_value_adc(t,k,tind) = &
                              stats%sum%count_obs_value_adc(t,k,tind) + 1
                      endif
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then
                            if(metric%selectOpt.eq.1 &
                                 .and.obs%value(t,k).ne.LVT_rc%udef) then 
                               stats%sum%obs_value_total(t,k,2) = &
                                    stats%sum%obs_value_total(t,k,2) + &
                                    obs%value(t,k)
                               stats%sum%count_obs_value_total(t,k,2) = &
                                    stats%sum%count_obs_value_total(t,k,2) + 1
                            endif
                            
                            if(metric%timeOpt.eq.1.and.&
                                 obs%value(t,k).ne.LVT_rc%udef) then 
                               stats%sum%obs_value_ts(t,k,2) = &
                                    stats%sum%obs_value_ts(t,k,2)+&
                                    obs%value(t,k)
                               stats%sum%count_obs_value_ts(t,k,2) = & 
                                    stats%sum%count_obs_value_ts(t,k,2)+1
                            endif
                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then
                            if(metric%selectOpt.eq.1.and.&
                                 obs%value(t,k).ne.LVT_rc%udef) then 
                               stats%sum%obs_value_total(t,k,3) = &
                                    stats%sum%obs_value_total(t,k,3) + &
                                    obs%value(t,k)
                               stats%sum%count_obs_value_total(t,k,3) = &
                                    stats%sum%count_obs_value_total(t,k,3) + 1
                            endif
                            
                            if(metric%timeOpt.eq.1.and.&
                                 obs%value(t,k).ne.LVT_rc%udef) then 
                               stats%sum%obs_value_ts(t,k,3) = &
                                    stats%sum%obs_value_ts(t,k,3)+&
                                    obs%value(t,k)
                               stats%sum%count_obs_value_ts(t,k,3) = & 
                                    stats%sum%count_obs_value_ts(t,k,3)+1
                            endif
                            
                         endif
                      endif
                   endif
                endif
             endif
          enddo
       enddo
    endif
  end subroutine diagnoseSingleSum

!BOP
! 
! !ROUTINE: LVT_computeSum
! \label{LVT_computeSum}
!
! !INTERFACE: 
  subroutine LVT_computeSum(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the sum values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelSum](\ref{computeSingleModelSum})
!     computes the sum values for a single variable
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
       if(LVT_metrics%sum%selectOpt.eq.1.or.&
            LVT_metrics%sum%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%sum%timeOpt.eq.1.and.&
                  LVT_metrics%sum%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%sum%ftn_ts_loc(i),200,advance='no') &
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
             
             call computeSingleSum(alarm,model,obs,stats,&
                  LVT_metrics%sum)
             
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%sum%timeOpt.eq.1.and.&
                  LVT_metrics%sum%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%sum%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_computeSum
  

!BOP
! 
! !ROUTINE: computeSingleSum
! \label{computeSingleSum}
!
! !INTERFACE: 
  subroutine computeSingleSum(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the sum values
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
                   if(stats%sum%count_model_value_ts(t,k,l).gt.0) then 
                      stats%sum%model_value_ts(t,k,l) = stats%sum%model_value_ts(t,k,l)
                      
                      stats%sum%tavg_model_value_ts(t,k,l) = & 
                           stats%sum%tavg_model_value_ts(t,k,l) + & 
                           stats%sum%model_value_ts(t,k,l)
                      stats%sum%tavg_count_model_value_ts(t,k,l) = & 
                           stats%sum%tavg_count_model_value_ts(t,k,l) + 1 
                   else
                      stats%sum%model_value_ts(t,k,l) = LVT_rc%udef
                   endif
                   if(LVT_rc%obssource(2).ne."none") then 
                      if(obs%selectNlevs.ge.1) then 
                         if(stats%sum%count_obs_value_ts(t,k,l).gt.0) then 
                            stats%sum%obs_value_ts(t,k,l) = &
                                 stats%sum%obs_value_ts(t,k,l)

                            stats%sum%tavg_obs_value_ts(t,k,l) = & 
                                 stats%sum%tavg_obs_value_ts(t,k,l) + & 
                                 stats%sum%obs_value_ts(t,k,l)
                            stats%sum%tavg_count_obs_value_ts(t,k,l) = & 
                                 stats%sum%tavg_count_obs_value_ts(t,k,l) + 1 
                         else
                            stats%sum%obs_value_ts(t,k,l) = LVT_rc%udef
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
                      if(stats%sum%tavg_count_model_value_ts(t,k,l).gt.0) then 
                         stats%sum%tavg_model_value_ts(t,k,l) = & 
                              stats%sum%tavg_model_value_ts(t,k,l) / & 
                              stats%sum%tavg_count_model_value_ts(t,k,l) 
                      endif
                   enddo
                enddo
             enddo

             if(metric%extractTS.eq.1) then 
                if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                   call LVT_writeTSinfo(metric%ftn_ts_loc,&
                        model,&
                        LVT_rc%ngrid,&
                        stats%sum%tavg_model_value_ts,&
                        stats%sum%tavg_count_model_value_ts,&
                        LVT_rc%ngrid,&
                        stats%sum%tavg_obs_value_ts,&
                        stats%sum%tavg_count_obs_value_ts)
                else
                   call LVT_writeTSinfo(metric%ftn_ts_loc,&
                        model,&
                        LVT_rc%ngrid,&
                        stats%sum%tavg_model_value_ts,&
                        stats%sum%tavg_count_model_value_ts)
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
                   if(stats%sum%count_model_value_total(t,k,l).gt.&
                        LVT_rc%obsCountThreshold) then 
                      stats%sum%model_value_total(t,k,l) = &
                           stats%sum%model_value_total(t,k,l)
                   else
                      stats%sum%model_value_total(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1, LVT_rc%nasc
                      if(stats%sum%count_model_value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then 
                         stats%sum%model_value_asc(t,k,l) = &
                              stats%sum%model_value_asc(t,k,l)
                      else
                         stats%sum%model_value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif

                if(metric%computeADC.eq.1) then 
                   do l=1, LVT_rc%nadc
                      if(stats%sum%count_model_value_adc(t,k,l).gt.&
                           LVT_rc%ADCCountThreshold) then 
                         stats%sum%model_value_adc(t,k,l) = &
                              stats%sum%model_value_adc(t,k,l)
                      else
                         stats%sum%model_value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
                if(LVT_rc%obssource(2).ne."none") then 
                   if(obs%selectNlevs.ge.1) then
                      do l=1,LVT_rc%strat_nlevels                      
                         if(stats%sum%count_obs_value_total(t,k,l).gt.&
                              LVT_rc%obsCountThreshold) then 
                            stats%sum%obs_value_total(t,k,l) = &
                                 stats%sum%obs_value_total(t,k,l)
                         else
                            stats%sum%obs_value_total(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                      if(metric%computeSC.eq.1) then
                         do l=1, LVT_rc%nasc 
                            if(stats%sum%count_obs_value_asc(t,k,l).gt.&
                                 LVT_rc%SCCountThreshold) then 
                               stats%sum%obs_value_asc(t,k,l) = &
                                    stats%sum%obs_value_asc(t,k,l)
                            else
                               stats%sum%obs_value_asc(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                      endif
                      
                      if(metric%computeADC.eq.1) then
                         do l=1, LVT_rc%nadc  
                            if(stats%sum%count_obs_value_adc(t,k,l).gt.&
                                 LVT_rc%ADCCountThreshold) then 
                               stats%sum%obs_value_adc(t,k,l) = &
                                    stats%sum%obs_value_adc(t,k,l)
                            else
                               stats%sum%obs_value_adc(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                      endif
                   endif
                endif
             enddo
          enddo
          
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%sum%model_value_total(:,k,l),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%sum%model_value_ci(k,l))
             enddo
          enddo
          if(LVT_rc%obssource(2).ne."none") then 
             do k=1,obs%selectNlevs
                do l=1, LVT_rc%strat_nlevels
                   call LVT_computeCI(stats%sum%obs_value_total(:,k,l),LVT_rc%ngrid,&
                        LVT_rc%pval_CI,stats%sum%obs_value_ci(k,l))
                enddo
             enddo
          endif

          if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
             call LVT_writeDataBasedStrat(model,obs,stats,metric,&
                  LVT_rc%ngrid, stats%sum%model_value_total, LVT_rc%ngrid,stats%sum%obs_value_total)
             
             if(metric%computeSC.eq.1) then 
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%sum%model_value_asc,stats%sum%count_model_value_asc,&
                     LVT_rc%ngrid,stats%sum%obs_value_asc,stats%sum%count_obs_value_asc)          
             endif
             if(metric%computeADC.eq.1) then 
                call LVT_writeAvgDiurnalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%sum%model_value_adc,stats%sum%count_model_value_adc,&
                     LVT_rc%ngrid,stats%sum%obs_value_adc,stats%sum%count_obs_value_adc)     
             endif
          else
             call LVT_writeDataBasedStrat(model,obs,stats,metric,&
                   LVT_rc%ngrid,stats%sum%model_value_total)

             if(metric%computeSC.eq.1) then 
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%sum%model_value_asc,stats%sum%count_model_value_asc)
             endif
             if(metric%computeADC.eq.1) then 
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%sum%model_value_adc,stats%sum%count_model_value_adc)
             endif
          endif
       endif
    endif
  end subroutine computeSingleSum

!BOP
! 
! !ROUTINE: LVT_writeMetric_Sum
! \label(LVT_writeMetric_Sum)
!
! !INTERFACE:
  subroutine LVT_writeMetric_Sum(pass,final,vlevels,stats,obs)
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

    if(pass.eq.LVT_metrics%sum%npass) then 
       if(final.ne.1) then
           if(stats%selectOpt.eq.1) then 
              do k=1,vlevels
                if(LVT_metrics%sum%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_gridded(LVT_metrics%sum%ftn_ts, &
                           stats%sum%model_value_ts(:,k,l),&
                           stats%vid_ts(LVT_Sumid,1),k)

                      call LVT_writevar_gridded(LVT_metrics%sum%ftn_ts, &
                           real(stats%sum%count_model_value_ts(:,k,l)),&
                           stats%vid_count_ts(LVT_Sumid,1),k)
                   enddo
                   if(LVT_rc%obssource(2).ne."none") then 
                      do l=1,LVT_rc%strat_nlevels
                         call LVT_writevar_gridded(LVT_metrics%sum%ftn_ts, &
                              stats%sum%obs_value_ts(:,k,l),&
                              stats%vid_ts(LVT_Sumid,2),k)
                         call LVT_writevar_gridded(LVT_metrics%sum%ftn_ts, &
                              real(stats%sum%count_obs_value_ts(:,k,l)),&
                              stats%vid_count_ts(LVT_Sumid,2),k)
                      enddo
                   endif
                endif
             enddo
          endif
       else
          if(pass.eq.1) then 
             if(stats%selectOpt.eq.1) then 
                do k=1,vlevels
                   if(LVT_metrics%sum%selectOpt.eq.1) then 
                      do l=1,LVT_rc%strat_nlevels                      
                         call LVT_writevar_gridded(LVT_metrics%sum%ftn_total, &
                              stats%sum%model_value_total(:,k,l),&
                              stats%vid_total(LVT_Sumid,1),k)
                         call LVT_writevar_gridded(LVT_metrics%sum%ftn_total, &
                              real(stats%sum%count_model_value_total(:,k,l)),&
                              stats%vid_count_total(LVT_Sumid,1),k)
                         
                      enddo
                      if(LVT_rc%obssource(2).ne."none".and.&
                           obs%selectNlevs.ge.1) then 
                         do l=1,LVT_rc%strat_nlevels                      
                            call LVT_writevar_gridded(LVT_metrics%sum%ftn_total, &
                                 stats%sum%obs_value_total(:,k,l),&
                                 stats%vid_total(LVT_Sumid,2),k)
                            call LVT_writevar_gridded(LVT_metrics%sum%ftn_total, &
                                 real(stats%sum%count_obs_value_total(:,k,l)),&
                                 stats%vid_count_total(LVT_Sumid,2),k)                         
                         enddo
                      endif
                   
                      if(LVT_metrics%sum%computeSC.eq.1) then 
                         do tind = 1,LVT_rc%nasc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%sum%ftn_total,&
                                 stats%sum%model_value_asc(:,k,tind),&
                                 stats%vid_sc_total(tind,LVT_Sumid,1),k)
                         enddo
                         if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                            do tind = 1,LVT_rc%nasc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%sum%ftn_total,&
                                    stats%sum%obs_value_asc(:,k,tind),&
                                    stats%vid_sc_total(tind,LVT_Sumid,2),k)
                            enddo
                         endif
                      endif
                      if(LVT_metrics%sum%computeADC.eq.1) then 
                         do tind = 1,LVT_rc%nadc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%sum%ftn_total,&
                                 stats%sum%model_value_adc(:,k,tind),&
                                 stats%vid_adc_total(tind,LVT_Sumid,1),k)
                         enddo
                         if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                            do tind = 1,LVT_rc%nadc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%sum%ftn_total,&
                                    stats%sum%obs_value_adc(:,k,tind),&
                                    stats%vid_adc_total(tind,LVT_Sumid,2),k)
                            enddo
                         endif
                      endif
                      call LVT_writeSummaryStats(&
                           LVT_metrics%sum%ftn_summ,&
                           LVT_metrics%sum%short_name,&
                           LVT_rc%ngrid,&
                           stats%sum%model_value_total(:,k,:), &
                           stats%sum%count_model_value_total(:,k,:),&
                           stats%standard_name,&
                           stats%sum%model_value_ci(k,:))
                      if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                         call LVT_writeSummaryStats(&
                              LVT_metrics%sum%ftn_summ,&
                              LVT_metrics%sum%short_name,&
                              LVT_rc%ngrid,&
                              stats%sum%obs_value_total(:,k,:), &
                              stats%sum%count_obs_value_total(:,k,:),&
                              "OBS_"//trim(stats%standard_name),&
                              stats%sum%obs_value_ci(k,:))
                      endif
                   endif
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_writeMetric_Sum

!BOP
! 
! !ROUTINE: LVT_resetMetric_Sum
! \label(LVT_resetMetric_Sum)
!
! !INTERFACE:
  subroutine LVT_resetMetric_Sum(alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
    logical                  :: alarm
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
             if(LVT_metrics%sum%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   stats%sum%model_value_ts(:,k,l) = 0.0
                   stats%sum%count_model_value_ts(:,k,l)=0 
                   if(alarm) then
                      stats%sum%tavg_model_value_ts(:,k,l) = 0.0
                      stats%sum%tavg_count_model_value_ts(:,k,l)=0 
                   endif
                enddo
                if(LVT_rc%obssource(2).ne."none".and.&
                     obs%selectNlevs.ge.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      stats%sum%obs_value_ts(:,k,l) = 0.0
                      stats%sum%count_obs_value_ts(:,k,l)=0 
                      if(alarm) then 
                         stats%sum%tavg_obs_value_ts(:,k,l) = 0.0
                         stats%sum%tavg_count_obs_value_ts(:,k,l)=0 
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

  end subroutine LVT_resetMetric_Sum

!BOP
! 
! !ROUTINE: LVT_writerestart_Sum
! 
! !INTERFACE:
  subroutine LVT_writerestart_Sum(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for Sum metric computations
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
       if(LVT_metrics%sum%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   call LVT_writevar_restart(ftn,&
                        stats%sum%model_value_total(:,k,l))
                   call LVT_writevar_restart(ftn,&
                        stats%sum%count_model_value_total(:,k,l))
                enddo
             enddo
             if(LVT_metrics%sum%computeSC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc         
                      call LVT_writevar_restart(ftn,&
                           stats%sum%model_value_asc(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%sum%count_model_value_asc(:,k,l))
                   enddo
                enddo
             endif
             if(LVT_metrics%sum%computeADC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nadc         
                      call LVT_writevar_restart(ftn,&
                           stats%sum%model_value_adc(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%sum%count_model_value_adc(:,k,l))
                   enddo
                enddo
             endif
             
             if(LVT_rc%obssource(2).ne."none") then 
                if(obs%selectNlevs.ge.1) then 
                   do k=1,obs%selectNlevs
                      do l=1,LVT_rc%strat_nlevels
                         call LVT_writevar_restart(ftn,&
                              stats%sum%obs_value_total(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%sum%count_obs_value_total(:,k,l))
                      enddo
                   enddo
                   
                   if(LVT_metrics%sum%computeSC.eq.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%nasc
                            call LVT_writevar_restart(ftn,&
                                 stats%sum%obs_value_asc(:,k,l))
                            call LVT_writevar_restart(ftn,&
                                 stats%sum%count_obs_value_asc(:,k,l))
                         enddo
                      enddo
                   endif
                   if(LVT_metrics%sum%computeADC.eq.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%nasc
                            call LVT_writevar_restart(ftn,&
                                 stats%sum%obs_value_adc(:,k,l))
                            call LVT_writevar_restart(ftn,&
                                 stats%sum%count_obs_value_adc(:,k,l))
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
  end subroutine LVT_writerestart_Sum


!BOP
! 
! !ROUTINE: LVT_readrestart_Sum
! 
! !INTERFACE:
  subroutine LVT_readrestart_Sum(ftn)
! !USES: 
! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for Sum metric computations
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
       if(LVT_metrics%sum%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 

             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   call LVT_readvar_restart(ftn,&
                        stats%sum%model_value_total(:,k,l))
                   call LVT_readvar_restart(ftn,&
                        stats%sum%count_model_value_total(:,k,l))
                enddo
             enddo
             if(LVT_metrics%sum%computeSC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc         
                      call LVT_readvar_restart(ftn,&
                           stats%sum%model_value_asc(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%sum%count_model_value_asc(:,k,l))
                   enddo
                enddo
             endif
             if(LVT_metrics%sum%computeADC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nadc         
                      call LVT_readvar_restart(ftn,&
                           stats%sum%model_value_adc(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%sum%count_model_value_adc(:,k,l))
                   enddo
                enddo
             endif
             
             if(LVT_rc%obssource(2).ne."none") then 
                if(obs%selectNlevs.ge.1) then 
                   do k=1,obs%selectNlevs
                      do l=1,LVT_rc%strat_nlevels
                         call LVT_readvar_restart(ftn,&
                              stats%sum%obs_value_total(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%sum%count_obs_value_total(:,k,l))
                      enddo
                   enddo
                   
                   if(LVT_metrics%sum%computeSC.eq.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%nasc
                            call LVT_readvar_restart(ftn,&
                                 stats%sum%obs_value_asc(:,k,l))
                            call LVT_readvar_restart(ftn,&
                                 stats%sum%count_obs_value_asc(:,k,l))
                         enddo
                      enddo
                   endif
                   if(LVT_metrics%sum%computeADC.eq.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%nasc
                            call LVT_readvar_restart(ftn,&
                                 stats%sum%obs_value_adc(:,k,l))
                            call LVT_readvar_restart(ftn,&
                                 stats%sum%count_obs_value_adc(:,k,l))
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
  end subroutine LVT_readrestart_Sum


end module LVT_SumMod
