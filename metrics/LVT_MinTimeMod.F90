!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_MintimeMod
! \label(LVT_MintimeMod)
!
! !INTERFACE:
module LVT_MintimeMod
! 
! !USES:  
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
  use LVT_TSMod
  use LVT_logMod
  use LVT_timeMgrMod
  use LVT_CIMod
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the computations required to compute minimum values
!  (temporally) of desired variables from the LIS output
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
  public :: LVT_initMintime
  public :: LVT_diagnoseMintime
  public :: LVT_computeMintime
  public :: LVT_writeMetric_Mintime
  public :: LVT_resetMetric_Mintime
  public :: LVT_writerestart_Mintime
  public :: LVT_readrestart_Mintime
!EOP
  
  type, public :: mintimedec
    character*5          :: time_option
  end type mintimedec

  type(mintimedec), save :: LVT_mintime_struc

  private

contains
  subroutine LVT_initMintime(model, obs, stats,metric)
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
! 
!EOP
    integer      :: rc
    if(metric%selectOpt.eq.1) then 

       call ESMF_ConfigGetAttribute(LVT_config,LVT_mintime_struc%time_option, &
            label="Time specification option for MinTime metric:",rc=rc)
       if(rc.ne.0) then 
          write(LVT_logunit,*) "[ERR] "
          write(LVT_logunit,*) "[ERR] Time specification option for MinTime metric: not specified"
          write(LVT_logunit,*) "[ERR] options are 'doy/da/hr/mn/ss' "
          write(LVT_logunit,*) "[ERR] "
          call LVT_endrun()
       endif

       allocate(stats%mintime%model_value_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%mintime%model_value_total = 1E10
       allocate(stats%mintime%model_value_min_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%mintime%model_value_min_total = 1E10
       allocate(stats%mintime%model_value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%mintime%model_value_ci = LVT_rc%udef
       
       if(metric%computeSC.eq.1) then 
          allocate(stats%mintime%model_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%mintime%model_value_asc = 1E10
          allocate(stats%mintime%model_value_min_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%mintime%model_value_min_asc = 1E10
       endif
       if(metric%computeADC.eq.1) then 
          allocate(stats%mintime%model_value_adc(LVT_rc%ngrid, obs%selectNlevs,&
               LVT_rc%nadc))
          stats%mintime%model_value_adc = 1E10
          allocate(stats%mintime%model_value_adc(LVT_rc%ngrid, obs%selectNlevs,&
               LVT_rc%nadc))
          stats%mintime%model_value_min_adc = 1E10
       endif
       
       if(LVT_rc%obssource(2).ne."none") then 
          if(obs%selectNlevs.ge.1) then 
             allocate(stats%mintime%obs_value_total(LVT_rc%ngrid, obs%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%mintime%obs_value_total = 1E10

             allocate(stats%mintime%obs_value_min_total(LVT_rc%ngrid, obs%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%mintime%obs_value_min_total = 1E10
             allocate(stats%mintime%obs_value_ci(obs%selectNlevs,LVT_rc%strat_nlevels))
             stats%mintime%obs_value_ci = LVT_rc%udef
             
             if(metric%computeSC.eq.1) then 
                allocate(stats%mintime%obs_value_asc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nasc))
                stats%mintime%obs_value_asc = 1E10
                allocate(stats%mintime%obs_value_min_asc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nasc))
                stats%mintime%obs_value_min_asc = 1E10
             endif
             if(metric%computeADC.eq.1) then 
                allocate(stats%mintime%obs_value_adc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nadc))
                stats%mintime%obs_value_adc = 1E10
                allocate(stats%mintime%obs_value_min_adc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nadc))
                stats%mintime%obs_value_min_adc = 1E10
             endif
          endif
       endif
    endif
    if(metric%timeOpt.eq.1) then 
       allocate(stats%mintime%model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%mintime%model_value_ts = 1E10

       allocate(stats%mintime%model_value_min_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%mintime%model_value_min_ts = 1E10

       if(LVT_rc%obssource(2).ne."none") then 
          if(obs%selectNlevs.ge.1) then 
             allocate(stats%mintime%obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%mintime%obs_value_ts = 1E10

             allocate(stats%mintime%obs_value_min_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%mintime%obs_value_min_ts = 1E10

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

  end subroutine LVT_initMintime

!BOP
! 
! !ROUTINE: LVT_diagnoseMintime
! \label{LVT_diagnoseMintime}
!
! !INTERFACE: 
  subroutine LVT_diagnoseMintime(pass)
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
!   calculating the mintime of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelMintime](\ref{diagnoseSingleModelMintime})
!     updates the mintime computation for a single variable 
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
       if(LVT_metrics%mintime%selectOpt.eq.1.or.&
            LVT_metrics%mintime%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleMintime(model,obs,stats,&
                  LVT_metrics%mintime)
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseMintime

!BOP
! 
! !ROUTINE: diagnoseSingleMintime
! \label{diagnoseSingleMintime}
!
! !INTERFACE: 
  subroutine diagnoseSingleMintime(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the mintime computation of the 
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
    integer    :: time

!----------------------------------------------------------------
!It is assumed temporal lag is not chosen in this setup
!----------------------------------------------------------------
    if(LVT_rc%tlag.gt.0) then 
       write(LVT_logunit,*) "[ERR] "
       write(LVT_logunit,*) "[ERR] Non-zero temporal lag specification is not "
       write(LVT_logunit,*) "[ERR] supported for MinTime metric"
       call LVT_endrun()
    endif
    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then        
       if(LVT_mintime_struc%time_option.eq."doy") then 
          time = LVT_rc%doy
       elseif(LVT_mintime_struc%time_option.eq."da") then 
          time = LVT_rc%da
       elseif(LVT_mintime_struc%time_option.eq."hr") then 
          time = LVT_rc%hr
       elseif(LVT_mintime_struc%time_option.eq."mn") then 
          time = LVT_rc%mn
       elseif(LVT_mintime_struc%time_option.eq."ss") then 
          time = LVT_rc%ss
       endif
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             if(model%count(t,k).gt.0) then
                if(metric%selectOpt.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      if(model%value(t,k).lt.stats%mintime%model_value_min_total(t,k,1)) then 
                         stats%mintime%model_value_min_total(t,k,1) = &
                              model%value(t,k)
                         stats%mintime%model_value_total(t,k,1) = &
                              time
                      endif
                   endif
                endif
                
                if(metric%timeOpt.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      if(model%value(t,k).lt.stats%mintime%model_value_min_ts(t,k,1)) then 
                         stats%mintime%model_value_min_ts(t,k,1) = model%value(t,k)
                         stats%mintime%model_value_ts(t,k,1) = time

                      endif
                   endif
                endif
                if(metric%computeSC.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                           tind)
                      if(model%value(t,k).lt.&
                           stats%mintime%model_value_min_asc(t,k,tind)) then
                         stats%mintime%model_value_min_asc(t,k,tind) = &
                              model%value(t,k)
                         stats%mintime%model_value_asc(t,k,tind) = &
                              time
                      endif
                   endif
                endif
                if(metric%computeADC.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      call LVT_getADCTimeIndex(tind)
                      if(model%value(t,k).lt.&
                           stats%mintime%model_value_min_adc(t,k,tind)) then 
                         stats%mintime%model_value_min_adc(t,k,tind) = & 
                              model%value(t,k)
                         stats%mintime%model_value_adc(t,k,tind) = &
                              time
                      endif
                   endif
                endif
                if(LVT_rc%strat_nlevels.gt.1) then 
                   if(LVT_stats%strat_var(t,k).gt.&
                        LVT_rc%strat_var_threshold) then
                      if(metric%selectOpt.eq.1) then
                         if(model%value(t,k).ne.LVT_rc%udef) then  
                            if(model%value(t,k).lt.&
                                 stats%mintime%model_value_min_total(t,k,2)) then
                               stats%mintime%model_value_min_total(t,k,2) = & 
                                    model%value(t,k)
                               stats%mintime%model_value_total(t,k,2) = &
                                    time
                            endif
                         endif
                      endif
                      if(metric%timeOpt.eq.1) then 
                         if(model%value(t,k).ne.LVT_rc%udef) then 
                            if(model%value(t,k).lt.&
                                 stats%mintime%model_value_min_ts(t,k,2)) then 
                               stats%mintime%model_value_min_ts(t,k,2) = & 
                                    model%value(t,k)
                               stats%mintime%model_value_ts(t,k,2) = &
                                    time
                            endif
                         endif
                      endif
                   elseif(LVT_stats%strat_var(t,k).le.&
                        LVT_rc%strat_var_threshold) then
                      if(metric%selectOpt.eq.1) then 
                         if(model%value(t,k).ne.LVT_rc%udef) then 
                            if(model%value(t,k).lt.&
                                 stats%mintime%model_value_min_total(t,k,3)) then 
                               stats%mintime%model_value_min_total(t,k,3) = & 
                                    model%value(t,k)
                               stats%mintime%model_value_total(t,k,3) = &
                                    time
                            endif
                         endif
                      endif
                      if(metric%timeOpt.eq.1) then 
                         if(model%value(t,k).ne.LVT_rc%udef) then 
                            if(model%value(t,k).lt.&
                                 stats%mintime%model_value_min_ts(t,k,3)) then 
                               stats%mintime%model_value_min_ts(t,k,3) = & 
                                    model%value(t,k)
                               stats%mintime%model_value_ts(t,k,3) = &
                                    time
                            endif
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
                         if(obs%value(t,k).lt.&
                              stats%mintime%obs_value_min_total(t,k,1)) then 
                            stats%mintime%obs_value_min_total(t,k,1) = & 
                                 obs%value(t,k)
                            stats%mintime%obs_value_total(t,k,1) = &
                                 time
                         endif
                      endif

                      if(metric%timeOpt.eq.1.and.obs%value(t,k).ne.LVT_rc%udef) then 
                         if(obs%value(t,k).lt.&
                              stats%mintime%obs_value_min_ts(t,k,1)) then 
                            stats%mintime%obs_value_min_ts(t,k,1) = obs%value(t,k)
                            stats%mintime%obs_value_ts(t,k,1) = time
                         endif
                      endif
                      if(metric%computeSC.eq.1.and.obs%value(t,k).ne.LVT_rc%udef) then 
                         call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,tind)
                         if(obs%value(t,k).lt.stats%mintime%obs_value_min_asc(t,k,tind)) then 
                            stats%mintime%obs_value_min_asc(t,k,tind) = &
                                 obs%value(t,k)
                            stats%mintime%obs_value_asc(t,k,tind) = &
                                 time
                         endif
                      endif
                      if(metric%computeADC.eq.1.and.obs%value(t,k).ne.LVT_rc%udef) then 
                         call LVT_getADCTimeIndex(tind)
                         if(obs%value(t,k).lt.stats%mintime%obs_value_min_adc(t,k,tind)) then 
                            stats%mintime%obs_value_min_adc(t,k,tind) = &
                                 obs%value(t,k)
                            stats%mintime%obs_value_adc(t,k,tind) = &
                                 time
                         endif
                      endif
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then
                            if(metric%selectOpt.eq.1 &
                                 .and.obs%value(t,k).ne.LVT_rc%udef) then 
                               if(obs%value(t,k).lt.stats%mintime%obs_value_min_total(t,k,2)) then 
                                  stats%mintime%obs_value_min_total(t,k,2) = &
                                       obs%value(t,k)
                                  stats%mintime%obs_value_total(t,k,2) = &
                                       time
                               endif
                            endif
                            
                            if(metric%timeOpt.eq.1.and.&
                                 obs%value(t,k).ne.LVT_rc%udef) then 
                               if(obs%value(t,k).lt.stats%mintime%obs_value_min_ts(t,k,2)) then 
                                  stats%mintime%obs_value_min_ts(t,k,2) = & 
                                       obs%value(t,k)
                                  stats%mintime%obs_value_ts(t,k,2) = &
                                       time
                               endif
                            endif
                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then
                            if(metric%selectOpt.eq.1.and.&
                                 obs%value(t,k).ne.LVT_rc%udef) then 
                               if(obs%value(t,k).lt.stats%mintime%obs_value_min_total(t,k,3)) then 
                                  stats%mintime%obs_value_min_total(t,k,3) = & 
                                       obs%value(t,k)
                                  stats%mintime%obs_value_total(t,k,3) = &
                                       time
                               endif
                            endif
                            
                            if(metric%timeOpt.eq.1.and.&
                                 obs%value(t,k).ne.LVT_rc%udef) then 
                               if(obs%value(t,k).lt.stats%mintime%obs_value_min_ts(t,k,3)) then 
                                  stats%mintime%obs_value_min_ts(t,k,3) = & 
                                       obs%value(t,k)
                                  stats%mintime%obs_value_ts(t,k,3) = &
                                       time
                               endif
                            endif
                         endif
                      endif
                   endif
                endif
             endif
          enddo
       enddo
    endif
  end subroutine diagnoseSingleMintime

!BOP
! 
! !ROUTINE: LVT_computeMintime
! \label{LVT_computeMintime}
!
! !INTERFACE: 
  subroutine LVT_computeMintime(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the mintime values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelMintime](\ref{computeSingleModelMintime})
!     computes the mintime values for a single variable
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
       if(LVT_metrics%mintime%selectOpt.eq.1.or.&
            LVT_metrics%mintime%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%mintime%timeOpt.eq.1.and.&
                  LVT_metrics%mintime%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%mintime%ftn_ts_loc(i),200,advance='no') &
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
             
             call computeSingleMintime(alarm,model,obs,stats,&
                  LVT_metrics%mintime)
             
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%mintime%timeOpt.eq.1.and.&
                  LVT_metrics%mintime%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%mintime%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_computeMintime
  

!BOP
! 
! !ROUTINE: computeSingleMintime
! \label{computeSingleMintime}
!
! !INTERFACE: 
  subroutine computeSingleMintime(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the mintime values
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
    integer     :: count_model_asc(LVT_rc%ngrid,model%selectNlevs,LVT_rc%nasc)
    integer     :: count_obs_asc(LVT_rc%ngrid,obs%selectNlevs,LVT_rc%nasc)
    integer     :: count_model_adc(LVT_rc%ngrid,model%selectNlevs,LVT_rc%nadc)
    integer     :: count_obs_adc(LVT_rc%ngrid,obs%selectNlevs,LVT_rc%nadc)
    
    integer    :: count_model(LVT_rc%ngrid, model%selectNlevs, LVT_rc%strat_nlevels)
    integer    :: count_obs(LVT_rc%ngrid, model%selectNlevs, LVT_rc%strat_nlevels)

    count_model = 1
    count_obs = 1
    count_model_asc = 1
    count_model_adc = 1
    count_obs_adc = 1

    if(metric%timeOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%mintime%model_value_ts(t,k,l).lt.1E10) then 
                      
                   else
                      stats%mintime%model_value_ts(t,k,l) = LVT_rc%udef
                   endif
                   if(LVT_rc%obssource(2).ne."none") then 
                      if(obs%selectNlevs.ge.1) then 
                         if(stats%mintime%obs_value_ts(t,k,l).lt.1E10) then 
                            
                         else
                            stats%mintime%obs_value_ts(t,k,l) = LVT_rc%udef
                         endif
                      endif
                   endif
                enddo
             enddo
          enddo
                   
          if(alarm) then 

             if(metric%extractTS.eq.1) then 
                if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                   call LVT_writeTSinfo(metric%ftn_ts_loc,&
                        model,&
                        LVT_rc%ngrid,&
                        stats%mintime%model_value_ts,&
                        count_model,&
                        LVT_rc%ngrid,&
                        stats%mintime%obs_value_ts,&
                        count_obs)
                else
                   call LVT_writeTSinfo(metric%ftn_ts_loc,&
                        model,&
                        LVT_rc%ngrid,&
                        stats%mintime%model_value_ts,&
                        count_model)
                endif

             end if
             
          endif
       end if
    end if

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%mintime%model_value_total(t,k,l).lt.1E10) then 
                      stats%mintime%model_value_total(t,k,l) = &
                           stats%mintime%model_value_total(t,k,l)
                   else
                      stats%mintime%model_value_total(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1, LVT_rc%nasc
                      if(stats%mintime%model_value_asc(t,k,l).lt.1E10) then 
                         stats%mintime%model_value_asc(t,k,l) = &
                              count_model_asc(t,1,l)           
                      else
                         stats%mintime%model_value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif

                if(metric%computeADC.eq.1) then 
                   do l=1, LVT_rc%nadc
                      if(stats%mintime%model_value_adc(t,k,l).lt.1E10) then 
                         stats%mintime%model_value_adc(t,k,l) = &
                              stats%mintime%model_value_adc(t,k,l)
                      else
                         stats%mintime%model_value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
                if(LVT_rc%obssource(2).ne."none") then 
                   if(obs%selectNlevs.ge.1) then
                      do l=1,LVT_rc%strat_nlevels                      
                         if(stats%mintime%obs_value_total(t,k,l).lt.1E10) then 
                            stats%mintime%obs_value_total(t,k,l) = &
                                 stats%mintime%obs_value_total(t,k,l)
                         else
                            stats%mintime%obs_value_total(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                      if(metric%computeSC.eq.1) then
                         do l=1, LVT_rc%nasc 
                            if(stats%mintime%obs_value_asc(t,k,l).lt.1E10) then 
                               stats%mintime%obs_value_asc(t,k,l) = &
                                    stats%mintime%obs_value_asc(t,k,l)
                            else
                               stats%mintime%obs_value_asc(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                      endif
                      
                      if(metric%computeADC.eq.1) then
                         do l=1, LVT_rc%nadc  
                            if(stats%mintime%obs_value_adc(t,k,l).lt.1E10) then
                               stats%mintime%obs_value_adc(t,k,l) = &
                                    stats%mintime%obs_value_adc(t,k,l)
                            else
                               stats%mintime%obs_value_adc(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                      endif
                   endif
                endif
             enddo
          enddo
          
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%mintime%model_value_total(:,k,l),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%mintime%model_value_ci(k,l))
             enddo
          enddo
          if(LVT_rc%obssource(2).ne."none") then 
             do k=1,obs%selectNlevs
                do l=1, LVT_rc%strat_nlevels
                   call LVT_computeCI(stats%mintime%obs_value_total(:,k,l),LVT_rc%ngrid,&
                        LVT_rc%pval_CI,stats%mintime%obs_value_ci(k,l))
                enddo
             enddo
          endif

          if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
             call LVT_writeDataBasedStrat(model,obs,stats,metric,&
                  LVT_rc%ngrid, stats%mintime%model_value_total, &
                  LVT_rc%ngrid,stats%mintime%obs_value_total)
             
             if(metric%computeSC.eq.1) then 
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%mintime%model_value_asc,count_model_asc,&
                     LVT_rc%ngrid,stats%mintime%obs_value_asc,count_obs_asc)          
             endif
             if(metric%computeADC.eq.1) then 
                call LVT_writeAvgDiurnalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%mintime%model_value_adc,count_model_adc,&
                     LVT_rc%ngrid,stats%mintime%obs_value_adc,count_obs_adc)     
             endif
          else
             call LVT_writeDataBasedStrat(model,obs,stats,metric,&
                   LVT_rc%ngrid,stats%mintime%model_value_total)

             if(metric%computeSC.eq.1) then 
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%mintime%model_value_asc,count_model_asc)
             endif
             if(metric%computeADC.eq.1) then 
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%mintime%model_value_adc,count_model_adc)
             endif
          endif
       endif
    endif
  end subroutine computeSingleMintime

!BOP
! 
! !ROUTINE: LVT_writeMetric_Mintime
! \label(LVT_writeMetric_Mintime)
!
! !INTERFACE:
  subroutine LVT_writeMetric_Mintime(pass,final,vlevels,stats,obs)
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
    integer    :: count_var(LVT_rc%ngrid,1,1)
    integer     :: count_obs_var(LVT_rc%ngrid,1,1)
    
    count_var = 1
    count_obs_var = 1

    if(pass.eq.LVT_metrics%mintime%npass) then 
       if(final.ne.1) then
           if(stats%selectOpt.eq.1) then 
              do k=1,vlevels
                if(LVT_metrics%mintime%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_gridded(LVT_metrics%mintime%ftn_ts, &
                           stats%mintime%model_value_ts(:,k,l),&
                           stats%vid_ts(LVT_Mintimeid,1),k)

                      call LVT_writevar_gridded(LVT_metrics%mintime%ftn_ts, &
                           real(count_var(:,1,l)),&
                           stats%vid_count_ts(LVT_Mintimeid,1),k)
                   enddo
                   if(LVT_rc%obssource(2).ne."none") then 
                      do l=1,LVT_rc%strat_nlevels
                         call LVT_writevar_gridded(LVT_metrics%mintime%ftn_ts, &
                              stats%mintime%obs_value_ts(:,k,l),&
                              stats%vid_ts(LVT_Mintimeid,2),k)
                         call LVT_writevar_gridded(LVT_metrics%mintime%ftn_ts, &
                              real(count_obs_var(:,1,l)),&
                              stats%vid_count_ts(LVT_Mintimeid,2),k)
                      enddo
                   endif
                endif
             enddo
          endif
       else
          if(pass.eq.1) then 
             if(stats%selectOpt.eq.1) then 
                do k=1,vlevels
                   if(LVT_metrics%mintime%selectOpt.eq.1) then 
                      do l=1,LVT_rc%strat_nlevels                      
                         call LVT_writevar_gridded(LVT_metrics%mintime%ftn_total, &
                              stats%mintime%model_value_total(:,k,l),&
                              stats%vid_total(LVT_Mintimeid,1),k)
                         call LVT_writevar_gridded(LVT_metrics%mintime%ftn_total, &
                              real(count_var(:,1,l)),&
                              stats%vid_count_total(LVT_Mintimeid,1),k)
                         
                      enddo
                      if(LVT_rc%obssource(2).ne."none".and.&
                           obs%selectNlevs.ge.1) then 
                         do l=1,LVT_rc%strat_nlevels                      
                            call LVT_writevar_gridded(LVT_metrics%mintime%ftn_total, &
                                 stats%mintime%obs_value_total(:,k,l),&
                                 stats%vid_total(LVT_Mintimeid,2),k)
                            call LVT_writevar_gridded(LVT_metrics%mintime%ftn_total, &
                                 real(count_obs_var(:,1,l)),&
                                 stats%vid_count_total(LVT_Mintimeid,2),k)                         
                         enddo
                      endif
                   
                      if(LVT_metrics%mintime%computeSC.eq.1) then 
                         do tind = 1,LVT_rc%nasc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%mintime%ftn_total,&
                                 stats%mintime%model_value_asc(:,k,tind),&
                                 stats%vid_sc_total(tind,LVT_Mintimeid,1),k)
                         enddo
                         if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                            do tind = 1,LVT_rc%nasc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%mintime%ftn_total,&
                                    stats%mintime%obs_value_asc(:,k,tind),&
                                    stats%vid_sc_total(tind,LVT_Mintimeid,2),k)
                            enddo
                         endif
                      endif
                      if(LVT_metrics%mintime%computeADC.eq.1) then 
                         do tind = 1,LVT_rc%nadc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%mintime%ftn_total,&
                                 stats%mintime%model_value_adc(:,k,tind),&
                                 stats%vid_adc_total(tind,LVT_Mintimeid,1),k)
                         enddo
                         if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                            do tind = 1,LVT_rc%nadc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%mintime%ftn_total,&
                                    stats%mintime%obs_value_adc(:,k,tind),&
                                    stats%vid_adc_total(tind,LVT_Mintimeid,2),k)
                            enddo
                         endif
                      endif
                      call LVT_writeSummaryStats(&
                           LVT_metrics%mintime%ftn_summ,&
                           LVT_metrics%mintime%short_name,&
                           LVT_rc%ngrid,&
                           stats%mintime%model_value_total(:,k,:), &
                           count_var(:,1,:),&
                           stats%standard_name,&
                           stats%mintime%model_value_ci(k,:))
                      if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                         call LVT_writeSummaryStats(&
                              LVT_metrics%mintime%ftn_summ,&
                              LVT_metrics%mintime%short_name,&
                              LVT_rc%ngrid,&
                              stats%mintime%obs_value_total(:,k,:), &
                              count_obs_var(:,1,:),&
                              "OBS_"//trim(stats%standard_name),&
                              stats%mintime%obs_value_ci(k,:))
                      endif
                   endif
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_writeMetric_Mintime

!BOP
! 
! !ROUTINE: LVT_resetMetric_Mintime
! \label(LVT_resetMetric_Mintime)
!
! !INTERFACE:
  subroutine LVT_resetMetric_Mintime(alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
    logical          :: alarm
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
             if(LVT_metrics%mintime%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   if(alarm) then 
                      stats%mintime%model_value_ts(:,k,l) = 1E10
                      stats%mintime%model_value_min_ts(:,k,l) = 1E10 
                   endif
                enddo
                if(LVT_rc%obssource(2).ne."none".and.&
                     obs%selectNlevs.ge.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      if(alarm) then 
                         stats%mintime%obs_value_ts(:,k,l) = 1E10
                         stats%mintime%obs_value_min_ts(:,k,l) = 1E10
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

  end subroutine LVT_resetMetric_Mintime

!BOP
! 
! !ROUTINE: LVT_writerestart_Mintime
! 
! !INTERFACE:
  subroutine LVT_writerestart_Mintime(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for Mintime metric computations
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
       if(LVT_metrics%mintime%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   call LVT_writevar_restart(ftn,&
                        stats%mintime%model_value_total(:,k,l))
                enddo
             enddo
             if(LVT_metrics%mintime%computeSC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc         
                      call LVT_writevar_restart(ftn,&
                           stats%mintime%model_value_asc(:,k,l))
                   enddo
                enddo
             endif
             if(LVT_metrics%mintime%computeADC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nadc         
                      call LVT_writevar_restart(ftn,&
                           stats%mintime%model_value_adc(:,k,l))
                   enddo
                enddo
             endif
             
             if(LVT_rc%obssource(2).ne."none") then 
                if(obs%selectNlevs.ge.1) then 
                   do k=1,obs%selectNlevs
                      do l=1,LVT_rc%strat_nlevels
                         call LVT_writevar_restart(ftn,&
                              stats%mintime%obs_value_total(:,k,l))
                      enddo
                   enddo
                   
                   if(LVT_metrics%mintime%computeSC.eq.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%nasc
                            call LVT_writevar_restart(ftn,&
                                 stats%mintime%obs_value_asc(:,k,l))
                         enddo
                      enddo
                   endif
                   if(LVT_metrics%mintime%computeADC.eq.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%nasc
                            call LVT_writevar_restart(ftn,&
                                 stats%mintime%obs_value_adc(:,k,l))
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
  end subroutine LVT_writerestart_Mintime


!BOP
! 
! !ROUTINE: LVT_readrestart_Mintime
! 
! !INTERFACE:
  subroutine LVT_readrestart_Mintime(ftn)
! !USES: 
! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for Mintime metric computations
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
       if(LVT_metrics%mintime%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 

             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   call LVT_readvar_restart(ftn,&
                        stats%mintime%model_value_total(:,k,l))
                enddo
             enddo
             if(LVT_metrics%mintime%computeSC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc         
                      call LVT_readvar_restart(ftn,&
                           stats%mintime%model_value_asc(:,k,l))
                   enddo
                enddo
             endif
             if(LVT_metrics%mintime%computeADC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nadc         
                      call LVT_readvar_restart(ftn,&
                           stats%mintime%model_value_adc(:,k,l))
                   enddo
                enddo
             endif
             
             if(LVT_rc%obssource(2).ne."none") then 
                if(obs%selectNlevs.ge.1) then 
                   do k=1,obs%selectNlevs
                      do l=1,LVT_rc%strat_nlevels
                         call LVT_readvar_restart(ftn,&
                              stats%mintime%obs_value_total(:,k,l))
                      enddo
                   enddo
                   
                   if(LVT_metrics%mintime%computeSC.eq.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%nasc
                            call LVT_readvar_restart(ftn,&
                                 stats%mintime%obs_value_asc(:,k,l))
                         enddo
                      enddo
                   endif
                   if(LVT_metrics%mintime%computeADC.eq.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%nasc
                            call LVT_readvar_restart(ftn,&
                                 stats%mintime%obs_value_adc(:,k,l))
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
  end subroutine LVT_readrestart_Mintime


end module LVT_MintimeMod
