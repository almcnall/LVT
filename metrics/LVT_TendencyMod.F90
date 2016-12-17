!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_TendencyMod
! \label(LVT_TendencyMod)
!
! !INTERFACE:
module LVT_TendencyMod
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
  public :: LVT_initTendency
  public :: LVT_diagnoseTendency
  public :: LVT_computeTendency
  public :: LVT_writeMetric_Tendency
  public :: LVT_resetMetric_Tendency
  public :: LVT_writerestart_Tendency
  public :: LVT_readrestart_Tendency
!EOP
  
  private

contains
  subroutine LVT_initTendency(model, obs, stats,metric)
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
       allocate(stats%tendency%model_value_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%tendency%model_value_total = 0.0
       allocate(stats%tendency%count_model_value_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%tendency%count_model_value_total = 0
       allocate(stats%tendency%model_value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%tendency%model_value_ci = LVT_rc%udef
       
       if(metric%computeSC.eq.1) then 
          allocate(stats%tendency%model_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%tendency%model_value_asc = 0.0
          allocate(stats%tendency%count_model_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%tendency%count_model_value_asc = 0
       endif
       if(metric%computeADC.eq.1) then 
          allocate(stats%tendency%model_value_adc(LVT_rc%ngrid, obs%selectNlevs,&
               LVT_rc%nadc))
          stats%tendency%model_value_adc = 0.0
          allocate(stats%tendency%count_model_value_adc(LVT_rc%ngrid, obs%selectNlevs,&
               LVT_rc%nadc))
          stats%tendency%count_model_value_adc = 0
       endif
       
       if(LVT_rc%obssource(2).ne."none") then 
          if(obs%selectNlevs.ge.1) then 
             allocate(stats%tendency%obs_value_total(LVT_rc%ngrid, obs%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%tendency%obs_value_total = 0.0
             allocate(stats%tendency%count_obs_value_total(LVT_rc%ngrid, obs%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%tendency%count_obs_value_total = 0
             allocate(stats%tendency%obs_value_ci(obs%selectNlevs,LVT_rc%strat_nlevels))
             stats%tendency%obs_value_ci = LVT_rc%udef
             
             if(metric%computeSC.eq.1) then 
                allocate(stats%tendency%obs_value_asc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nasc))
                stats%tendency%obs_value_asc = 0.0
                allocate(stats%tendency%count_obs_value_asc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nasc))
                stats%tendency%count_obs_value_asc = 0
             endif
             if(metric%computeADC.eq.1) then 
                allocate(stats%tendency%obs_value_adc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nadc))
                stats%tendency%obs_value_adc = 0.0
                allocate(stats%tendency%count_obs_value_adc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nadc))
                stats%tendency%count_obs_value_adc = 0
             endif
          endif
       endif
       
       allocate(stats%tendency%model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%tendency%model_value_ts = 0.0
       allocate(stats%tendency%count_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%tendency%count_model_value_ts = 0
       
       allocate(stats%tendency%tavg_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%tendency%tavg_model_value_ts = 0.0
       allocate(stats%tendency%tavg_count_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%tendency%tavg_count_model_value_ts = 0
       
       allocate(stats%tendency%model_value_cval_ts(LVT_rc%ngrid,model%selectNlevs, &
            LVT_rc%strat_nlevels))
       allocate(stats%tendency%model_value_pval_ts(LVT_rc%ngrid,model%selectNlevs, &
            LVT_rc%strat_nlevels))
       
       stats%tendency%model_value_cval_ts = LVT_rc%udef
       stats%tendency%model_value_pval_ts = LVT_rc%udef

       if(LVT_rc%obssource(2).ne."none") then 
          if(obs%selectNlevs.ge.1) then 
             allocate(stats%tendency%obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%tendency%obs_value_ts = 0.0
             allocate(stats%tendency%count_obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%tendency%count_obs_value_ts = 0
             
             allocate(stats%tendency%tavg_obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%tendency%tavg_obs_value_ts = 0.0
             allocate(stats%tendency%tavg_count_obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%tendency%tavg_count_obs_value_ts = 0

             allocate(stats%tendency%obs_value_cval_ts(LVT_rc%ngrid,obs%selectNlevs, &
                  LVT_rc%strat_nlevels))
             allocate(stats%tendency%obs_value_pval_ts(LVT_rc%ngrid,obs%selectNlevs, &
                  LVT_rc%strat_nlevels))
             
             stats%tendency%obs_value_cval_ts = LVT_rc%udef
             stats%tendency%obs_value_pval_ts = LVT_rc%udef
             
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

  end subroutine LVT_initTendency

!BOP
! 
! !ROUTINE: LVT_diagnoseTendency
! \label{LVT_diagnoseTendency}
!
! !INTERFACE: 
  subroutine LVT_diagnoseTendency(pass)
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
!   calculating the tendency of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelTendency](\ref{diagnoseSingleModelTendency})
!     updates the tendency computation for a single variable 
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
       if(LVT_metrics%tendency%selectOpt.eq.1.or.&
            LVT_metrics%tendency%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleTendency(model,obs,stats,&
                  LVT_metrics%tendency)
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseTendency

!BOP
! 
! !ROUTINE: diagnoseSingleTendency
! \label{diagnoseSingleTendency}
!
! !INTERFACE: 
  subroutine diagnoseSingleTendency(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the tendency computation of the 
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
             stats%tendency%model_value_pval_ts(t,k,1) = &
                  stats%tendency%model_value_cval_ts(t,k,1) 
             stats%tendency%model_value_cval_ts(t,k,1) = &
                  model%value(t,k)
             
             if(stats%tendency%model_value_cval_ts(t,k,1).ne.LVT_rc%udef.and.&
                  stats%tendency%model_value_pval_ts(t,k,1).ne.LVT_rc%udef) then
                stats%tendency%model_value_ts(t,k,1) = &
                     (stats%tendency%model_value_cval_ts(t,k,1) - &
                     stats%tendency%model_value_pval_ts(t,k,1))
                stats%tendency%count_model_value_ts(t,k,1) = & 
                     stats%tendency%count_model_value_ts(t,k,1)+1
             endif
             if(metric%computeSC.eq.1) then 
                if(stats%tendency%model_value_ts(t,k,1).ne.LVT_rc%udef) then 
                   call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                        tind)
                   stats%tendency%model_value_asc(t,k,tind) = &
                        stats%tendency%model_value_asc(t,k,tind)+&
                        stats%tendency%model_value_ts(t,k,1)
                   stats%tendency%count_model_value_asc(t,k,tind) = &
                        stats%tendency%count_model_value_asc(t,k,tind) + 1
                endif
             endif
             if(metric%computeADC.eq.1) then 
                if(stats%tendency%model_value_ts(t,k,1).ne.LVT_rc%udef) then 
                   call LVT_getADCTimeIndex(tind)
                   stats%tendency%model_value_adc(t,k,tind) = &
                        stats%tendency%model_value_adc(t,k,tind)+&
                        stats%tendency%model_value_ts(t,k,1)
                   stats%tendency%count_model_value_adc(t,k,tind) = &
                        stats%tendency%count_model_value_adc(t,k,tind) + 1
                endif
             endif
             if(LVT_rc%strat_nlevels.gt.1) then 
                if(LVT_stats%strat_var(t,k).gt.&
                     LVT_rc%strat_var_threshold) then
                   stats%tendency%model_value_pval_ts(t,k,2) = &
                        stats%tendency%model_value_cval_ts(t,k,2) 
                   stats%tendency%model_value_cval_ts(t,k,2) = &
                        model%value(t,k)
                   
                   if(stats%tendency%model_value_cval_ts(t,k,2).ne.LVT_rc%udef.and.&
                        stats%tendency%model_value_pval_ts(t,k,2).ne.LVT_rc%udef) then
                      stats%tendency%model_value_ts(t,k,2) = &
                           (stats%tendency%model_value_cval_ts(t,k,2) - &
                           stats%tendency%model_value_pval_ts(t,k,2))
                      
                      stats%tendency%count_model_value_ts(t,k,2) = & 
                           stats%tendency%count_model_value_ts(t,k,2)+1
                   endif
                elseif(LVT_stats%strat_var(t,k).le.&
                     LVT_rc%strat_var_threshold) then
                   stats%tendency%model_value_pval_ts(t,k,3) = &
                        stats%tendency%model_value_cval_ts(t,k,3) 
                   stats%tendency%model_value_cval_ts(t,k,3) = &
                        model%value(t,k)
                   
                   if(stats%tendency%model_value_cval_ts(t,k,3).ne.LVT_rc%udef.and.&
                        stats%tendency%model_value_pval_ts(t,k,3).ne.LVT_rc%udef) then
                      stats%tendency%model_value_ts(t,k,3) = &
                           (stats%tendency%model_value_cval_ts(t,k,3) - &
                           stats%tendency%model_value_pval_ts(t,k,3))
                      
                      stats%tendency%count_model_value_ts(t,k,3) = & 
                           stats%tendency%count_model_value_ts(t,k,3)+1
                   endif
                endif
             endif
          end do
          do k=1,obs%selectNlevs
             if(LVT_rc%obssource(2).ne."none") then
                if(obs%selectNlevs.ge.1) then

                   stats%tendency%obs_value_pval_ts(t,k,1) = &
                        stats%tendency%obs_value_cval_ts(t,k,1) 
                   stats%tendency%obs_value_cval_ts(t,k,1) = &
                        obs%value(t,k)
                   
                   if(stats%tendency%obs_value_cval_ts(t,k,1).ne.LVT_rc%udef.and.&
                        stats%tendency%obs_value_pval_ts(t,k,1).ne.LVT_rc%udef) then
                      stats%tendency%obs_value_ts(t,k,1) = &
                           (stats%tendency%obs_value_cval_ts(t,k,1) - &
                           stats%tendency%obs_value_pval_ts(t,k,1))
                      
                      stats%tendency%count_obs_value_ts(t,k,1) = & 
                           stats%tendency%count_obs_value_ts(t,k,1)+1
                   endif
                   if(metric%computeSC.eq.1) then 
                      if(stats%tendency%obs_value_ts(t,k,1).ne.LVT_rc%udef) then 
                         call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,tind)
                         stats%tendency%obs_value_asc(t,k,tind) = &
                              stats%tendency%obs_value_asc(t,k,tind)+&
                              stats%tendency%obs_value_ts(t,k,1)
                         stats%tendency%count_obs_value_asc(t,k,tind) = &
                              stats%tendency%count_obs_value_asc(t,k,tind) + 1
                      endif
                   endif
                   if(metric%computeADC.eq.1.and.&
                        stats%tendency%obs_value_ts(t,k,1).ne.LVT_rc%udef) then 
                      call LVT_getADCTimeIndex(tind)
                      stats%tendency%obs_value_adc(t,k,tind) = &
                           stats%tendency%obs_value_adc(t,k,tind)+&
                           stats%tendency%obs_value_ts(t,k,1)
                      stats%tendency%count_obs_value_adc(t,k,tind) = &
                           stats%tendency%count_obs_value_adc(t,k,tind) + 1
                   endif
                   if(LVT_rc%strat_nlevels.gt.1) then 
                      if(LVT_stats%strat_var(t,k).gt.&
                           LVT_rc%strat_var_threshold) then
                         stats%tendency%obs_value_pval_ts(t,k,2) = &
                              stats%tendency%obs_value_cval_ts(t,k,2) 
                         stats%tendency%obs_value_cval_ts(t,k,2) = &
                              obs%value(t,k)
                         
                         if(stats%tendency%obs_value_cval_ts(t,k,2).ne.LVT_rc%udef.and.&
                              stats%tendency%obs_value_pval_ts(t,k,2).ne.LVT_rc%udef) then
                            stats%tendency%obs_value_ts(t,k,2) = &
                                 (stats%tendency%obs_value_cval_ts(t,k,2) - &
                                 stats%tendency%obs_value_pval_ts(t,k,2))
                            
                            stats%tendency%count_obs_value_ts(t,k,2) = & 
                                 stats%tendency%count_obs_value_ts(t,k,2)+1
                         endif
                      elseif(LVT_stats%strat_var(t,k).le.&
                           LVT_rc%strat_var_threshold) then
                         stats%tendency%obs_value_pval_ts(t,k,3) = &
                              stats%tendency%obs_value_cval_ts(t,k,3) 
                         stats%tendency%obs_value_cval_ts(t,k,3) = &
                              obs%value(t,k)
                         
                         if(stats%tendency%obs_value_cval_ts(t,k,3).ne.LVT_rc%udef.and.&
                              stats%tendency%obs_value_pval_ts(t,k,3).ne.LVT_rc%udef) then
                            stats%tendency%obs_value_ts(t,k,3) = &
                                 (stats%tendency%obs_value_cval_ts(t,k,3) - &
                                 stats%tendency%obs_value_pval_ts(t,k,3))
                            
                            stats%tendency%count_obs_value_ts(t,k,3) = & 
                                 stats%tendency%count_obs_value_ts(t,k,3)+1
                         endif
                      endif
                   endif
                end if
             end if
          end do
       end do
    end if
  end subroutine diagnoseSingleTendency

!BOP
! 
! !ROUTINE: LVT_computeTendency
! \label{LVT_computeTendency}
!
! !INTERFACE: 
  subroutine LVT_computeTendency(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the tendency values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelTendency](\ref{computeSingleModelTendency})
!     computes the tendency values for a single variable
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
       if(LVT_metrics%tendency%selectOpt.eq.1.or.&
            LVT_metrics%tendency%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%tendency%timeOpt.eq.1.and.&
                  LVT_metrics%tendency%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%tendency%ftn_ts_loc(i),200,advance='no') &
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
             
             call computeSingleTendency(alarm,model,obs,stats,&
                  LVT_metrics%tendency)
             
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%tendency%timeOpt.eq.1.and.&
                  LVT_metrics%tendency%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%tendency%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_computeTendency
  

!BOP
! 
! !ROUTINE: computeSingleTendency
! \label{computeSingleTendency}
!
! !INTERFACE: 
  subroutine computeSingleTendency(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the tendency values
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
                   if(stats%tendency%count_model_value_ts(t,k,l).gt.0) then 
                      stats%tendency%model_value_ts(t,k,l) = &
                           stats%tendency%model_value_ts(t,k,l)/&
                           stats%tendency%count_model_value_ts(t,k,l)                  
                      stats%tendency%tavg_model_value_ts(t,k,l) = &
                           stats%tendency%tavg_model_value_ts(t,k,l) + &
                           stats%tendency%model_value_ts(t,k,l)
                      stats%tendency%tavg_count_model_value_ts(t,k,l) = & 
                           stats%tendency%tavg_count_model_value_ts(t,k,l) + 1
                   else
                      stats%tendency%model_value_ts(t,k,l) = LVT_rc%udef
                   endif
                   if(LVT_rc%obssource(2).ne."none") then 
                      if(obs%selectNlevs.ge.1) then 
                         if(stats%tendency%count_obs_value_ts(t,k,l).gt.0) then 
                            stats%tendency%obs_value_ts(t,k,l) = &
                                 stats%tendency%obs_value_ts(t,k,l)/&
                                 stats%tendency%count_obs_value_ts(t,k,l)

                            stats%tendency%tavg_obs_value_ts(t,k,l) = &
                                 stats%tendency%tavg_obs_value_ts(t,k,l) + & 
                                 stats%tendency%obs_value_ts(t,k,l)
                            stats%tendency%tavg_count_obs_value_ts(t,k,l) = &
                                 stats%tendency%tavg_count_obs_value_ts(t,k,l) + 1 
                         else
                            stats%tendency%obs_value_ts(t,k,l) = LVT_rc%udef
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
                      if(stats%tendency%tavg_count_model_value_ts(t,k,l).gt.0) then 
                         stats%tendency%tavg_model_value_ts(t,k,l) = &
                              stats%tendency%tavg_model_value_ts(t,k,l)/&
                           stats%tendency%tavg_count_model_value_ts(t,k,l) 
                         
                      else
                         stats%tendency%tavg_model_value_ts(t,k,l) = LVT_rc%udef
                      endif
                      if(LVT_rc%obssource(2).ne."none") then 
                         if(obs%selectNlevs.ge.1) then 
                            if(stats%tendency%tavg_count_obs_value_ts(t,k,l).gt.0) then 
                               stats%tendency%tavg_obs_value_ts(t,k,l) = &
                                    stats%tendency%tavg_obs_value_ts(t,k,l)/&
                                    stats%tendency%tavg_count_obs_value_ts(t,k,l)
                            else
                               stats%tendency%tavg_obs_value_ts(t,k,l) = LVT_rc%udef
                            endif
                         endif
                      endif
                   end do
                enddo
             enddo

             if(metric%extractTS.eq.1) then 
                if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                   call LVT_writeTSinfo(metric%ftn_ts_loc,&
                        model,&
                        LVT_rc%ngrid,&
                        stats%tendency%tavg_model_value_ts,&
                        stats%tendency%tavg_count_model_value_ts,&
                        LVT_rc%ngrid,&
                        stats%tendency%tavg_obs_value_ts,&
                        stats%tendency%tavg_count_obs_value_ts)
                else
                   call LVT_writeTSinfo(metric%ftn_ts_loc,&
                        model,&
                        LVT_rc%ngrid,&
                        stats%tendency%tavg_model_value_ts,&
                        stats%tendency%tavg_count_model_value_ts)
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
                   if(stats%tendency%count_model_value_total(t,k,l).gt.&
                        LVT_rc%obsCountThreshold) then 
                      stats%tendency%model_value_total(t,k,l) = &
                           stats%tendency%model_value_total(t,k,l)/&
                              stats%tendency%count_model_value_total(t,k,l)           
                   else
                      stats%tendency%model_value_total(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1, LVT_rc%nasc
                      if(stats%tendency%count_model_value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then 
                         stats%tendency%model_value_asc(t,k,l) = &
                              stats%tendency%model_value_asc(t,k,l)/&
                              stats%tendency%count_model_value_asc(t,k,l)           
                      else
                         stats%tendency%model_value_asc(t,k,l) = LVT_rc%udef
                      endif                      
                   enddo
                endif

                if(metric%computeADC.eq.1) then 
                   do l=1, LVT_rc%nadc
                      if(stats%tendency%count_model_value_adc(t,k,l).gt.&
                           LVT_rc%ADCCountThreshold) then 
                         stats%tendency%model_value_adc(t,k,l) = &
                              stats%tendency%model_value_adc(t,k,l)/&
                              stats%tendency%count_model_value_adc(t,k,l)           
                      else
                         stats%tendency%model_value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
                if(LVT_rc%obssource(2).ne."none") then 
                   if(obs%selectNlevs.ge.1) then
                      do l=1,LVT_rc%strat_nlevels                      
                         if(stats%tendency%count_obs_value_total(t,k,l).gt.&
                              LVT_rc%obsCountThreshold) then 
                            stats%tendency%obs_value_total(t,k,l) = &
                                 stats%tendency%obs_value_total(t,k,l)/&
                                 stats%tendency%count_obs_value_total(t,k,l)       
                         else
                            stats%tendency%obs_value_total(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                      if(metric%computeSC.eq.1) then
                         do l=1, LVT_rc%nasc 
                            if(stats%tendency%count_obs_value_asc(t,k,l).gt.&
                                 LVT_rc%SCCountThreshold) then 
                               stats%tendency%obs_value_asc(t,k,l) = &
                                    stats%tendency%obs_value_asc(t,k,l)/&
                                    stats%tendency%count_obs_value_asc(t,k,l)           
                            else
                               stats%tendency%obs_value_asc(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                      endif
                      
                      if(metric%computeADC.eq.1) then
                         do l=1, LVT_rc%nadc  
                            if(stats%tendency%count_obs_value_adc(t,k,l).gt.&
                                 LVT_rc%ADCCountThreshold) then 
                               stats%tendency%obs_value_adc(t,k,l) = &
                                    stats%tendency%obs_value_adc(t,k,l)/&
                                    stats%tendency%count_obs_value_adc(t,k,l)           
                            else
                               stats%tendency%obs_value_adc(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                      endif
                   endif
                endif
             enddo
          enddo
          
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%tendency%model_value_total(:,k,l),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%tendency%model_value_ci(k,l))
             enddo
          enddo
          if(LVT_rc%obssource(2).ne."none") then 
             do k=1,obs%selectNlevs
                do l=1, LVT_rc%strat_nlevels
                   call LVT_computeCI(stats%tendency%obs_value_total(:,k,l),LVT_rc%ngrid,&
                        LVT_rc%pval_CI,stats%tendency%obs_value_ci(k,l))
                enddo
             enddo
          endif

          if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
             call LVT_writeDataBasedStrat(model,obs,stats,metric,&
                  LVT_rc%ngrid, stats%tendency%model_value_total, &
                  LVT_rc%ngrid,stats%tendency%obs_value_total)
             
             if(metric%computeSC.eq.1) then 
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%tendency%model_value_asc,&
                     stats%tendency%count_model_value_asc,&
                     LVT_rc%ngrid,stats%tendency%obs_value_asc,&
                     stats%tendency%count_obs_value_asc)          
             endif
             if(metric%computeADC.eq.1) then 
                call LVT_writeAvgDiurnalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%tendency%model_value_adc,&
                     stats%tendency%count_model_value_adc,&
                     LVT_rc%ngrid,stats%tendency%obs_value_adc,&
                     stats%tendency%count_obs_value_adc)     
             endif
          else
             call LVT_writeDataBasedStrat(model,obs,stats,metric,&
                   LVT_rc%ngrid,stats%tendency%model_value_total)

             if(metric%computeSC.eq.1) then 
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%tendency%model_value_asc,&
                     stats%tendency%count_model_value_asc)
             endif
             if(metric%computeADC.eq.1) then 
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%tendency%model_value_adc,&
                     stats%tendency%count_model_value_adc)
             endif
          endif
       endif
    endif
  end subroutine computeSingleTendency

!BOP
! 
! !ROUTINE: LVT_writeMetric_Tendency
! \label(LVT_writeMetric_Tendency)
!
! !INTERFACE:
  subroutine LVT_writeMetric_Tendency(pass,final,vlevels,stats,obs)
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

    if(pass.eq.LVT_metrics%tendency%npass) then 
       if(final.ne.1) then
           if(stats%selectOpt.eq.1) then 
              do k=1,vlevels
                if(LVT_metrics%tendency%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_gridded(LVT_metrics%tendency%ftn_ts, &
                           stats%tendency%tavg_model_value_ts(:,k,l),&
                           stats%vid_ts(LVT_Tendencyid,1),k)
                      
                      call LVT_writevar_gridded(LVT_metrics%tendency%ftn_ts, &
                           real(stats%tendency%tavg_count_model_value_ts(:,k,l)),&
                           stats%vid_count_ts(LVT_Tendencyid,1),k)
                   enddo

                   if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                      do l=1,LVT_rc%strat_nlevels
                         call LVT_writevar_gridded(LVT_metrics%tendency%ftn_ts, &
                              stats%tendency%tavg_obs_value_ts(:,k,l),&
                              stats%vid_ts(LVT_Tendencyid,2),k)
                         call LVT_writevar_gridded(LVT_metrics%tendency%ftn_ts, &
                              real(stats%tendency%tavg_count_obs_value_ts(:,k,l)),&
                              stats%vid_count_ts(LVT_Tendencyid,2),k)
                      enddo
                   endif
                endif
             enddo
          endif
       else
          if(pass.eq.1) then 
             if(stats%selectOpt.eq.1) then 
                do k=1,vlevels
                   if(LVT_metrics%tendency%selectOpt.eq.1) then 
                      do l=1,LVT_rc%strat_nlevels                      
                         call LVT_writevar_gridded(LVT_metrics%tendency%ftn_total, &
                              stats%tendency%model_value_total(:,k,l),&
                              stats%vid_total(LVT_Tendencyid,1),k)
                         call LVT_writevar_gridded(LVT_metrics%tendency%ftn_total, &
                              real(stats%tendency%count_model_value_total(:,k,l)),&
                              stats%vid_count_total(LVT_Tendencyid,1),k)
                         
                      enddo
                      if(LVT_rc%obssource(2).ne."none".and.&
                           obs%selectNlevs.ge.1) then 
                         do l=1,LVT_rc%strat_nlevels                      
                            call LVT_writevar_gridded(LVT_metrics%tendency%ftn_total, &
                                 stats%tendency%obs_value_total(:,k,l),&
                                 stats%vid_total(LVT_Tendencyid,2),k)
                            call LVT_writevar_gridded(LVT_metrics%tendency%ftn_total, &
                                 real(stats%tendency%count_obs_value_total(:,k,l)),&
                                 stats%vid_count_total(LVT_Tendencyid,2),k)                         
                         enddo
                      endif
                   
                      if(LVT_metrics%tendency%computeSC.eq.1) then 
                         do tind = 1,LVT_rc%nasc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%tendency%ftn_total,&
                                 stats%tendency%model_value_asc(:,k,tind),&
                                 stats%vid_sc_total(tind,LVT_Tendencyid,1),k)
                         enddo
                         if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                            do tind = 1,LVT_rc%nasc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%tendency%ftn_total,&
                                    stats%tendency%obs_value_asc(:,k,tind),&
                                    stats%vid_sc_total(tind,LVT_Tendencyid,2),k)
                            enddo
                         endif
                      endif
                      if(LVT_metrics%tendency%computeADC.eq.1) then 
                         do tind = 1,LVT_rc%nadc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%tendency%ftn_total,&
                                 stats%tendency%model_value_adc(:,k,tind),&
                                 stats%vid_adc_total(tind,LVT_Tendencyid,1),k)
                         enddo
                         if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                            do tind = 1,LVT_rc%nadc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%tendency%ftn_total,&
                                    stats%tendency%obs_value_adc(:,k,tind),&
                                    stats%vid_adc_total(tind,LVT_Tendencyid,2),k)
                            enddo
                         endif
                      endif
                      call LVT_writeSummaryStats(&
                           LVT_metrics%tendency%ftn_summ,&
                           LVT_metrics%tendency%short_name,&
                           LVT_rc%ngrid,&
                           stats%tendency%model_value_total(:,k,:), &
                           stats%tendency%count_model_value_total(:,k,:),&
                           stats%standard_name,&
                           stats%tendency%model_value_ci(k,:))
                      if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                         call LVT_writeSummaryStats(&
                              LVT_metrics%tendency%ftn_summ,&
                              LVT_metrics%tendency%short_name,&
                              LVT_rc%ngrid,&
                              stats%tendency%obs_value_total(:,k,:), &
                              stats%tendency%count_obs_value_total(:,k,:),&
                              "OBS_"//trim(stats%standard_name),&
                              stats%tendency%obs_value_ci(k,:))
                      endif
                   endif
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_writeMetric_Tendency

!BOP
! 
! !ROUTINE: LVT_resetMetric_Tendency
! \label(LVT_resetMetric_Tendency)
!
! !INTERFACE:
  subroutine LVT_resetMetric_Tendency(alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
    logical                :: alarm
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
             do l=1,LVT_rc%strat_nlevels
                stats%tendency%count_model_value_ts(:,k,l) = 0 
                if(alarm) then 
                   stats%tendency%tavg_model_value_ts(:,k,l) = 0.0
                   stats%tendency%tavg_count_model_value_ts(:,k,l)=0 
                endif
             enddo
             
             if(LVT_rc%obssource(2).ne."none".and.&
                  obs%selectNlevs.ge.1) then 
                do l=1,LVT_rc%strat_nlevels
                   stats%tendency%count_obs_value_ts(:,k,l) = 0 
                   if(alarm) then 
                      stats%tendency%tavg_obs_value_ts(:,k,l) = 0.0
                      stats%tendency%tavg_count_obs_value_ts(:,k,l)=0 
                   endif
                enddo
             endif
          enddo
       endif
       
       model => model%next
       obs => obs%next
       stats => stats%next
    enddo

  end subroutine LVT_resetMetric_Tendency

!BOP
! 
! !ROUTINE: LVT_writerestart_Tendency
! 
! !INTERFACE:
  subroutine LVT_writerestart_Tendency(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for Tendency metric computations
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
       if(LVT_metrics%tendency%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 
             if(LVT_metrics%tendency%computeSC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc         
                      call LVT_writevar_restart(ftn,&
                           stats%tendency%model_value_asc(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%tendency%count_model_value_asc(:,k,l))
                   enddo
                enddo
             endif
             if(LVT_metrics%tendency%computeADC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nadc         
                      call LVT_writevar_restart(ftn,&
                           stats%tendency%model_value_adc(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%tendency%count_model_value_adc(:,k,l))
                   enddo
                enddo
             endif
             
             if(LVT_rc%obssource(2).ne."none") then 
                if(obs%selectNlevs.ge.1) then 
                   
                   if(LVT_metrics%tendency%computeSC.eq.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%nasc
                            call LVT_writevar_restart(ftn,&
                                 stats%tendency%obs_value_asc(:,k,l))
                            call LVT_writevar_restart(ftn,&
                                 stats%tendency%count_obs_value_asc(:,k,l))
                         enddo
                      enddo
                   endif
                   if(LVT_metrics%tendency%computeADC.eq.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%nasc
                            call LVT_writevar_restart(ftn,&
                                 stats%tendency%obs_value_adc(:,k,l))
                            call LVT_writevar_restart(ftn,&
                                 stats%tendency%count_obs_value_adc(:,k,l))
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
  end subroutine LVT_writerestart_Tendency


!BOP
! 
! !ROUTINE: LVT_readrestart_Tendency
! 
! !INTERFACE:
  subroutine LVT_readrestart_Tendency(ftn)
! !USES: 
! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for Tendency metric computations
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
       if(LVT_metrics%tendency%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 

             if(LVT_metrics%tendency%computeSC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc         
                      call LVT_readvar_restart(ftn,&
                           stats%tendency%model_value_asc(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%tendency%count_model_value_asc(:,k,l))
                   enddo
                enddo
             endif
             if(LVT_metrics%tendency%computeADC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nadc         
                      call LVT_readvar_restart(ftn,&
                           stats%tendency%model_value_adc(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%tendency%count_model_value_adc(:,k,l))
                   enddo
                enddo
             endif
             
             if(LVT_rc%obssource(2).ne."none") then 
                if(obs%selectNlevs.ge.1) then 
                   
                   if(LVT_metrics%tendency%computeSC.eq.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%nasc
                            call LVT_readvar_restart(ftn,&
                                 stats%tendency%obs_value_asc(:,k,l))
                            call LVT_readvar_restart(ftn,&
                                 stats%tendency%count_obs_value_asc(:,k,l))
                         enddo
                      enddo
                   endif
                   if(LVT_metrics%tendency%computeADC.eq.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%nasc
                            call LVT_readvar_restart(ftn,&
                                 stats%tendency%obs_value_adc(:,k,l))
                            call LVT_readvar_restart(ftn,&
                                 stats%tendency%count_obs_value_adc(:,k,l))
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
  end subroutine LVT_readrestart_Tendency


end module LVT_TendencyMod
