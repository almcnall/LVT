!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_MEANMod
! \label(LVT_MEANMod)
!
! !INTERFACE:
module LVT_MEANMod
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
  public :: LVT_initMEAN
  public :: LVT_diagnoseMEAN
  public :: LVT_computeMEAN
  public :: LVT_writeMetric_MEAN
  public :: LVT_resetMetric_MEAN
  public :: LVT_writerestart_MEAN
  public :: LVT_readrestart_MEAN
!EOP
  
  private

contains
  subroutine LVT_initMEAN(model, obs, stats,metric)
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
       allocate(stats%mean%model_value_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%mean%model_value_total = 0.0
       allocate(stats%mean%count_model_value_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%mean%count_model_value_total = 0
       allocate(stats%mean%model_value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%mean%model_value_ci = LVT_rc%udef
       
       if(metric%computeSC.eq.1) then 
          allocate(stats%mean%model_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%mean%model_value_asc = 0.0
          allocate(stats%mean%count_model_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%mean%count_model_value_asc = 0
       endif
       if(metric%computeADC.eq.1) then 
          allocate(stats%mean%model_value_adc(LVT_rc%ngrid, obs%selectNlevs,&
               LVT_rc%nadc))
          stats%mean%model_value_adc = 0.0
          allocate(stats%mean%count_model_value_adc(LVT_rc%ngrid, obs%selectNlevs,&
               LVT_rc%nadc))
          stats%mean%count_model_value_adc = 0
       endif
       
       if(LVT_rc%obssource(2).ne."none") then 
          if(obs%selectNlevs.ge.1) then 
             allocate(stats%mean%obs_value_total(LVT_rc%ngrid, obs%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%mean%obs_value_total = 0.0
             allocate(stats%mean%count_obs_value_total(LVT_rc%ngrid, obs%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%mean%count_obs_value_total = 0
             allocate(stats%mean%obs_value_ci(obs%selectNlevs,LVT_rc%strat_nlevels))
             stats%mean%obs_value_ci = LVT_rc%udef
             
             if(metric%computeSC.eq.1) then 
                allocate(stats%mean%obs_value_asc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nasc))
                stats%mean%obs_value_asc = 0.0
                allocate(stats%mean%count_obs_value_asc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nasc))
                stats%mean%count_obs_value_asc = 0
             endif
             if(metric%computeADC.eq.1) then 
                allocate(stats%mean%obs_value_adc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nadc))
                stats%mean%obs_value_adc = 0.0
                allocate(stats%mean%count_obs_value_adc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nadc))
                stats%mean%count_obs_value_adc = 0
             endif
          endif
       endif
    endif
    if(metric%timeOpt.eq.1) then 
       allocate(stats%mean%model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%mean%model_value_ts = 0.0
       allocate(stats%mean%count_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%mean%count_model_value_ts = 0

       allocate(stats%mean%tavg_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%mean%tavg_model_value_ts = 0.0
       allocate(stats%mean%tavg_count_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%mean%tavg_count_model_value_ts = 0
       
       if(LVT_rc%obssource(2).ne."none") then 
          if(obs%selectNlevs.ge.1) then 
             allocate(stats%mean%obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%mean%obs_value_ts = 0.0
             allocate(stats%mean%count_obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%mean%count_obs_value_ts = 0


             allocate(stats%mean%tavg_obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%mean%tavg_obs_value_ts = 0.0
             allocate(stats%mean%tavg_count_obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%mean%tavg_count_obs_value_ts = 0

             if(obs%stdev_flag) then 
                allocate(stats%mean%obs_value_stdev_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
                stats%mean%obs_value_stdev_ts = 0.0
                allocate(stats%mean%count_obs_value_stdev_ts(LVT_rc%ngrid, model%selectNlevs, &
                     LVT_rc%strat_nlevels))
                stats%mean%count_obs_value_stdev_ts = 0

                allocate(stats%mean%tavg_obs_value_stdev_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
                stats%mean%tavg_obs_value_stdev_ts = 0.0
                allocate(stats%mean%tavg_count_obs_value_stdev_ts(LVT_rc%ngrid, model%selectNlevs, &
                     LVT_rc%strat_nlevels))
                stats%mean%tavg_count_obs_value_stdev_ts = 0
             endif
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

  end subroutine LVT_initMEAN

!BOP
! 
! !ROUTINE: LVT_diagnoseMEAN
! \label{LVT_diagnoseMEAN}
!
! !INTERFACE: 
  subroutine LVT_diagnoseMEAN(pass)
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
!   calculating the mean of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelMEAN](\ref{diagnoseSingleModelMEAN})
!     updates the mean computation for a single variable 
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
       if(LVT_metrics%mean%selectOpt.eq.1.or.&
            LVT_metrics%mean%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleMEAN(model,obs,stats,&
                  LVT_metrics%mean)
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseMEAN

!BOP
! 
! !ROUTINE: diagnoseSingleMEAN
! \label{diagnoseSingleMEAN}
!
! !INTERFACE: 
  subroutine diagnoseSingleMEAN(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the mean computation of the 
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
                      stats%mean%model_value_total(t,k,1) = &
                           stats%mean%model_value_total(t,k,1) + &
                           model%value(t,k)
                      stats%mean%count_model_value_total(t,k,1) = &
                           stats%mean%count_model_value_total(t,k,1) + 1
                   endif
                endif
                
                if(metric%timeOpt.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      stats%mean%model_value_ts(t,k,1) = stats%mean%model_value_ts(t,k,1)+&
                           model%value(t,k)
                      stats%mean%count_model_value_ts(t,k,1) = & 
                           stats%mean%count_model_value_ts(t,k,1)+1

                   endif
                endif
                if(metric%computeSC.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                           tind)
                      stats%mean%model_value_asc(t,k,tind) = &
                           stats%mean%model_value_asc(t,k,tind)+&
                           model%value(t,k)
                      stats%mean%count_model_value_asc(t,k,tind) = &
                           stats%mean%count_model_value_asc(t,k,tind) + 1
                   endif
                endif
                if(metric%computeADC.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      call LVT_getADCTimeIndex(tind)
                      stats%mean%model_value_adc(t,k,tind) = &
                           stats%mean%model_value_adc(t,k,tind)+&
                           model%value(t,k)
                      stats%mean%count_model_value_adc(t,k,tind) = &
                           stats%mean%count_model_value_adc(t,k,tind) + 1
                   endif
                endif
                if(LVT_rc%strat_nlevels.gt.1) then 
                   if(LVT_stats%strat_var(t,k).gt.&
                        LVT_rc%strat_var_threshold) then
                      if(metric%selectOpt.eq.1) then
                         if(model%value(t,k).ne.LVT_rc%udef) then  
                            stats%mean%model_value_total(t,k,2) = &
                                 stats%mean%model_value_total(t,k,2) + &
                                 model%value(t,k)
                            stats%mean%count_model_value_total(t,k,2) = &
                                 stats%mean%count_model_value_total(t,k,2) + 1
                         endif
                      endif
                      if(metric%timeOpt.eq.1) then 
                         if(model%value(t,k).ne.LVT_rc%udef) then 
                            stats%mean%model_value_ts(t,k,2) = &
                                 stats%mean%model_value_ts(t,k,2)+&
                                 model%value(t,k)
                            stats%mean%count_model_value_ts(t,k,2) = & 
                                 stats%mean%count_model_value_ts(t,k,2)+1
                         endif
                      endif
                   elseif(LVT_stats%strat_var(t,k).le.&
                        LVT_rc%strat_var_threshold) then
                      if(metric%selectOpt.eq.1) then 
                         if(model%value(t,k).ne.LVT_rc%udef) then 
                            stats%mean%model_value_total(t,k,3) = &
                                 stats%mean%model_value_total(t,k,3) + &
                                 model%value(t,k)
                            stats%mean%count_model_value_total(t,k,3) = &
                                 stats%mean%count_model_value_total(t,k,3) + 1
                         endif
                      endif
                      if(metric%timeOpt.eq.1) then 
                         if(model%value(t,k).ne.LVT_rc%udef) then 
                            stats%mean%model_value_ts(t,k,3) = &
                                 stats%mean%model_value_ts(t,k,3)+&
                                 model%value(t,k)
                            stats%mean%count_model_value_ts(t,k,3) = & 
                                 stats%mean%count_model_value_ts(t,k,3)+1
                         endif
                      endif
                   endif
                endif
             endif
          end do
          do k=1,obs%selectNlevs
             if(LVT_rc%obssource(2).ne."none") then
                if(obs%selectNlevs.ge.1) then
                   if(obs%count(t,k).gt.0) then 
                      if(metric%selectOpt.eq.1.and.&
                           obs%value(t,k).ne.LVT_rc%udef) then 
                         stats%mean%obs_value_total(t,k,1) = &
                              stats%mean%obs_value_total(t,k,1) + &
                              obs%value(t,k)
                         stats%mean%count_obs_value_total(t,k,1) = &
                              stats%mean%count_obs_value_total(t,k,1) + 1
                      endif
                      
                      if(metric%timeOpt.eq.1.and.obs%value(t,k).ne.LVT_rc%udef) then 
                         stats%mean%obs_value_ts(t,k,1) = stats%mean%obs_value_ts(t,k,1)+&
                              obs%value(t,k)
                         stats%mean%count_obs_value_ts(t,k,1) = & 
                              stats%mean%count_obs_value_ts(t,k,1)+1
                         if(obs%stdev_flag) then 
                            stats%mean%obs_value_stdev_ts(t,k,1) = &
                                 stats%mean%obs_value_stdev_ts(t,k,1)+&
                                 obs%stdev(t,k)
                            stats%mean%count_obs_value_stdev_ts(t,k,1) = & 
                                 stats%mean%count_obs_value_stdev_ts(t,k,1)+1
                         endif
                      endif
                      if(metric%computeSC.eq.1.and.obs%value(t,k).ne.LVT_rc%udef) then 
                         call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,tind)
                         stats%mean%obs_value_asc(t,k,tind) = &
                              stats%mean%obs_value_asc(t,k,tind)+&
                              obs%value(t,k)
                         stats%mean%count_obs_value_asc(t,k,tind) = &
                              stats%mean%count_obs_value_asc(t,k,tind) + 1
                      endif
                      if(metric%computeADC.eq.1.and.obs%value(t,k).ne.LVT_rc%udef) then 
                         call LVT_getADCTimeIndex(tind)
                         stats%mean%obs_value_adc(t,k,tind) = &
                              stats%mean%obs_value_adc(t,k,tind)+&
                              obs%value(t,k)
                         stats%mean%count_obs_value_adc(t,k,tind) = &
                              stats%mean%count_obs_value_adc(t,k,tind) + 1
                      endif
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then
                            if(metric%selectOpt.eq.1 &
                                 .and.obs%value(t,k).ne.LVT_rc%udef) then 
                               stats%mean%obs_value_total(t,k,2) = &
                                    stats%mean%obs_value_total(t,k,2) + &
                                    obs%value(t,k)
                               stats%mean%count_obs_value_total(t,k,2) = &
                                    stats%mean%count_obs_value_total(t,k,2) + 1
                            endif
                            
                            if(metric%timeOpt.eq.1.and.&
                                 obs%value(t,k).ne.LVT_rc%udef) then 
                               stats%mean%obs_value_ts(t,k,2) = &
                                    stats%mean%obs_value_ts(t,k,2)+&
                                    obs%value(t,k)
                               stats%mean%count_obs_value_ts(t,k,2) = & 
                                    stats%mean%count_obs_value_ts(t,k,2)+1
                               if(obs%stdev_flag) then 
                                  stats%mean%obs_value_stdev_ts(t,k,2) = &
                                       stats%mean%obs_value_stdev_ts(t,k,2)+&
                                       obs%stdev(t,k)
                                  stats%mean%count_obs_value_stdev_ts(t,k,2) = & 
                                       stats%mean%count_obs_value_stdev_ts(t,k,2)+1
                               endif
                            endif
                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then
                            if(metric%selectOpt.eq.1.and.&
                                 obs%value(t,k).ne.LVT_rc%udef) then 
                               stats%mean%obs_value_total(t,k,3) = &
                                    stats%mean%obs_value_total(t,k,3) + &
                                    obs%value(t,k)
                               stats%mean%count_obs_value_total(t,k,3) = &
                                    stats%mean%count_obs_value_total(t,k,3) + 1
                            endif
                            
                            if(metric%timeOpt.eq.1.and.&
                                 obs%value(t,k).ne.LVT_rc%udef) then 
                               stats%mean%obs_value_ts(t,k,3) = &
                                    stats%mean%obs_value_ts(t,k,3)+&
                                    obs%value(t,k)
                               stats%mean%count_obs_value_ts(t,k,3) = & 
                                    stats%mean%count_obs_value_ts(t,k,3)+1
                               if(obs%stdev_flag) then 
                                  stats%mean%obs_value_stdev_ts(t,k,3) = &
                                       stats%mean%obs_value_stdev_ts(t,k,3)+&
                                       obs%stdev(t,k)
                                  stats%mean%count_obs_value_stdev_ts(t,k,3) = & 
                                       stats%mean%count_obs_value_stdev_ts(t,k,3)+1
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

  end subroutine diagnoseSingleMEAN

!BOP
! 
! !ROUTINE: LVT_computeMEAN
! \label{LVT_computeMEAN}
!
! !INTERFACE: 
  subroutine LVT_computeMEAN(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the mean values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelMEAN](\ref{computeSingleModelMEAN})
!     computes the mean values for a single variable
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
       if(LVT_metrics%mean%selectOpt.eq.1.or.&
            LVT_metrics%mean%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%mean%timeOpt.eq.1.and.&
                  LVT_metrics%mean%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%mean%ftn_ts_loc(i),200,advance='no') &
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
             
             call computeSingleMEAN(alarm,model,obs,stats,&
                  LVT_metrics%mean)
             
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%mean%timeOpt.eq.1.and.&
                  LVT_metrics%mean%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%mean%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_computeMEAN
  

!BOP
! 
! !ROUTINE: computeSingleMEAN
! \label{computeSingleMEAN}
!
! !INTERFACE: 
  subroutine computeSingleMEAN(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the mean values
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
                   if(stats%mean%count_model_value_ts(t,k,l).gt.0) then 
                      stats%mean%model_value_ts(t,k,l) =&
                           stats%mean%model_value_ts(t,k,l)/&
                           stats%mean%count_model_value_ts(t,k,l)                   

                      stats%mean%tavg_model_value_ts(t,k,l) = &
                           stats%mean%tavg_model_value_ts(t,k,l) + &
                           stats%mean%model_value_ts(t,k,l)
                      stats%mean%tavg_count_model_value_ts(t,k,l) = & 
                           stats%mean%tavg_count_model_value_ts(t,k,l) + 1
                   else
                      stats%mean%model_value_ts(t,k,l) = LVT_rc%udef
                   endif
                   if(LVT_rc%obssource(2).ne."none") then 
                      if(obs%selectNlevs.ge.1) then 
                         if(stats%mean%count_obs_value_ts(t,k,l).gt.0) then 
                            stats%mean%obs_value_ts(t,k,l) = &
                                 stats%mean%obs_value_ts(t,k,l)/&
                                 stats%mean%count_obs_value_ts(t,k,l)

                            stats%mean%tavg_obs_value_ts(t,k,l) = &
                                 stats%mean%tavg_obs_value_ts(t,k,l) + & 
                                 stats%mean%obs_value_ts(t,k,l)
                            stats%mean%tavg_count_obs_value_ts(t,k,l) = &
                                 stats%mean%tavg_count_obs_value_ts(t,k,l) + 1 
                         else
                            stats%mean%obs_value_ts(t,k,l) = LVT_rc%udef
                         endif
                         if(obs%stdev_flag) then 
                            if(stats%mean%count_obs_value_stdev_ts(t,k,l).gt.0) then 
                               stats%mean%obs_value_stdev_ts(t,k,l) = &
                                    stats%mean%obs_value_stdev_ts(t,k,l)/&
                                    stats%mean%count_obs_value_stdev_ts(t,k,l)

                               stats%mean%tavg_obs_value_stdev_ts(t,k,l) = &
                                    stats%mean%tavg_obs_value_stdev_ts(t,k,l) + & 
                                    stats%mean%obs_value_stdev_ts(t,k,l)
                               stats%mean%tavg_count_obs_value_stdev_ts(t,k,l) = &
                                    stats%mean%tavg_count_obs_value_stdev_ts(t,k,l) + 1
                            else
                               stats%mean%obs_value_stdev_ts(t,k,l) = LVT_rc%udef
                            endif
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
                      if(stats%mean%tavg_count_model_value_ts(t,k,l).gt.0) then 
                         stats%mean%tavg_model_value_ts(t,k,l) = &
                              stats%mean%tavg_model_value_ts(t,k,l)/&
                           stats%mean%tavg_count_model_value_ts(t,k,l)                   
                      else
                         stats%mean%tavg_model_value_ts(t,k,l) = LVT_rc%udef
                      endif
                      if(LVT_rc%obssource(2).ne."none") then 
                         if(obs%selectNlevs.ge.1) then 
                            if(stats%mean%tavg_count_obs_value_ts(t,k,l).gt.0) then 
                               stats%mean%tavg_obs_value_ts(t,k,l) = &
                                    stats%mean%tavg_obs_value_ts(t,k,l)/&
                                    stats%mean%tavg_count_obs_value_ts(t,k,l)
                            else
                               stats%mean%tavg_obs_value_ts(t,k,l) = LVT_rc%udef
                            endif
                            if(obs%stdev_flag) then 
                               if(stats%mean%tavg_count_obs_value_stdev_ts(t,k,l).gt.0) then 
                                  stats%mean%tavg_obs_value_stdev_ts(t,k,l) = &
                                       stats%mean%tavg_obs_value_stdev_ts(t,k,l)/&
                                       stats%mean%tavg_count_obs_value_stdev_ts(t,k,l)
                               else
                                  stats%mean%tavg_obs_value_stdev_ts(t,k,l) = LVT_rc%udef
                               endif
                            endif
                         endif
                      endif
                   end do
                enddo
             enddo

             if(metric%extractTS.eq.1) then 
                if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                   if(obs%stdev_flag) then 
                      call LVT_writeTSinfo(metric%ftn_ts_loc,&
                           model,&
                           LVT_rc%ngrid,&
                           stats%mean%tavg_model_value_ts,&
                           stats%mean%tavg_count_model_value_ts,&
                           LVT_rc%ngrid,&
                           stats%mean%tavg_obs_value_ts,&
                           stats%mean%tavg_count_obs_value_ts, &
                           stats%mean%tavg_obs_value_stdev_ts, &
                           stats%mean%tavg_count_obs_value_stdev_ts)
                   else
                      call LVT_writeTSinfo(metric%ftn_ts_loc,&
                           model,&
                           LVT_rc%ngrid,&
                           stats%mean%tavg_model_value_ts,&
                           stats%mean%tavg_count_model_value_ts,&
                           LVT_rc%ngrid,&
                           stats%mean%tavg_obs_value_ts,&
                           stats%mean%tavg_count_obs_value_ts)
                   endif
                else
                   call LVT_writeTSinfo(metric%ftn_ts_loc,&
                        model,&
                        LVT_rc%ngrid,&
                        stats%mean%tavg_model_value_ts,&
                        stats%mean%tavg_count_model_value_ts)
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
                   if(stats%mean%count_model_value_total(t,k,l).gt.&
                        LVT_rc%obsCountThreshold) then 
                      stats%mean%model_value_total(t,k,l) = &
                           stats%mean%model_value_total(t,k,l)/&
                              stats%mean%count_model_value_total(t,k,l)           
                   else
                      stats%mean%model_value_total(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1, LVT_rc%nasc
                      if(stats%mean%count_model_value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then 
                         stats%mean%model_value_asc(t,k,l) = &
                              stats%mean%model_value_asc(t,k,l)/&
                              stats%mean%count_model_value_asc(t,k,l)           
                      else
                         stats%mean%model_value_asc(t,k,l) = LVT_rc%udef
                      endif                      
                   enddo
                endif

                if(metric%computeADC.eq.1) then 
                   do l=1, LVT_rc%nadc
                      if(stats%mean%count_model_value_adc(t,k,l).gt.&
                           LVT_rc%ADCCountThreshold) then 
                         stats%mean%model_value_adc(t,k,l) = &
                              stats%mean%model_value_adc(t,k,l)/&
                              stats%mean%count_model_value_adc(t,k,l)           
                      else
                         stats%mean%model_value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
                if(LVT_rc%obssource(2).ne."none") then 
                   if(obs%selectNlevs.ge.1) then
                      do l=1,LVT_rc%strat_nlevels                      
                         if(stats%mean%count_obs_value_total(t,k,l).gt.&
                              LVT_rc%obsCountThreshold) then 
                            stats%mean%obs_value_total(t,k,l) = &
                                 stats%mean%obs_value_total(t,k,l)/&
                                 stats%mean%count_obs_value_total(t,k,l)       
                         else
                            stats%mean%obs_value_total(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                      if(metric%computeSC.eq.1) then
                         do l=1, LVT_rc%nasc 
                            if(stats%mean%count_obs_value_asc(t,k,l).gt.&
                                 LVT_rc%SCCountThreshold) then 
                               stats%mean%obs_value_asc(t,k,l) = &
                                    stats%mean%obs_value_asc(t,k,l)/&
                                    stats%mean%count_obs_value_asc(t,k,l)           
                            else
                               stats%mean%obs_value_asc(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                      endif
                      
                      if(metric%computeADC.eq.1) then
                         do l=1, LVT_rc%nadc  
                            if(stats%mean%count_obs_value_adc(t,k,l).gt.&
                                 LVT_rc%ADCCountThreshold) then 
                               stats%mean%obs_value_adc(t,k,l) = &
                                    stats%mean%obs_value_adc(t,k,l)/&
                                    stats%mean%count_obs_value_adc(t,k,l)           
                            else
                               stats%mean%obs_value_adc(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                      endif
                   endif
                endif
             enddo
          enddo
          
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%mean%model_value_total(:,k,l),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%mean%model_value_ci(k,l))
             enddo
          enddo
          if(LVT_rc%obssource(2).ne."none") then 
             do k=1,obs%selectNlevs
                do l=1, LVT_rc%strat_nlevels
                   call LVT_computeCI(stats%mean%obs_value_total(:,k,l),LVT_rc%ngrid,&
                        LVT_rc%pval_CI,stats%mean%obs_value_ci(k,l))
                enddo
             enddo
          endif

          if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
             call LVT_writeDataBasedStrat(model,obs,stats,metric,&
                  LVT_rc%ngrid, stats%mean%model_value_total, LVT_rc%ngrid,stats%mean%obs_value_total)
             
             if(metric%computeSC.eq.1) then 
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%mean%model_value_asc,stats%mean%count_model_value_asc,&
                     LVT_rc%ngrid,stats%mean%obs_value_asc,stats%mean%count_obs_value_asc)          
             endif
             if(metric%computeADC.eq.1) then 
                call LVT_writeAvgDiurnalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%mean%model_value_adc,stats%mean%count_model_value_adc,&
                     LVT_rc%ngrid,stats%mean%obs_value_adc,stats%mean%count_obs_value_adc)     
             endif
          else
             call LVT_writeDataBasedStrat(model,obs,stats,metric,&
                   LVT_rc%ngrid,stats%mean%model_value_total)

             if(metric%computeSC.eq.1) then 
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%mean%model_value_asc,stats%mean%count_model_value_asc)
             endif
             if(metric%computeADC.eq.1) then 
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,stats%mean%model_value_adc,stats%mean%count_model_value_adc)
             endif
          endif
       endif
    endif
  end subroutine computeSingleMEAN

!BOP
! 
! !ROUTINE: LVT_writeMetric_MEAN
! \label(LVT_writeMetric_MEAN)
!
! !INTERFACE:
  subroutine LVT_writeMetric_MEAN(pass,final,vlevels,stats,obs)
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

    if(pass.eq.LVT_metrics%mean%npass) then 
       if(final.ne.1) then
           if(stats%selectOpt.eq.1) then 
              do k=1,vlevels
                if(LVT_metrics%mean%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_gridded(LVT_metrics%mean%ftn_ts, &
                           stats%mean%tavg_model_value_ts(:,k,l),&
                           stats%vid_ts(LVT_MEANid,1),k)
                      
                      call LVT_writevar_gridded(LVT_metrics%mean%ftn_ts, &
                           real(stats%mean%tavg_count_model_value_ts(:,k,l)),&
                           stats%vid_count_ts(LVT_MEANid,1),k)
                   enddo

                   if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                      do l=1,LVT_rc%strat_nlevels
                         call LVT_writevar_gridded(LVT_metrics%mean%ftn_ts, &
                              stats%mean%tavg_obs_value_ts(:,k,l),&
                              stats%vid_ts(LVT_MEANid,2),k)
                         call LVT_writevar_gridded(LVT_metrics%mean%ftn_ts, &
                              real(stats%mean%tavg_count_obs_value_ts(:,k,l)),&
                              stats%vid_count_ts(LVT_MEANid,2),k)
                      enddo
                   endif
                endif
             enddo
          endif
       else
          if(pass.eq.1) then 
             if(stats%selectOpt.eq.1) then 
                do k=1,vlevels
                   if(LVT_metrics%mean%selectOpt.eq.1) then 
                      do l=1,LVT_rc%strat_nlevels                      
                         call LVT_writevar_gridded(LVT_metrics%mean%ftn_total, &
                              stats%mean%model_value_total(:,k,l),&
                              stats%vid_total(LVT_MEANid,1),k)
                         call LVT_writevar_gridded(LVT_metrics%mean%ftn_total, &
                              real(stats%mean%count_model_value_total(:,k,l)),&
                              stats%vid_count_total(LVT_MEANid,1),k)
                         
                      enddo
                      if(LVT_rc%obssource(2).ne."none".and.&
                           obs%selectNlevs.ge.1) then 
                         do l=1,LVT_rc%strat_nlevels                      
                            call LVT_writevar_gridded(LVT_metrics%mean%ftn_total, &
                                 stats%mean%obs_value_total(:,k,l),&
                                 stats%vid_total(LVT_MEANid,2),k)
                            call LVT_writevar_gridded(LVT_metrics%mean%ftn_total, &
                                 real(stats%mean%count_obs_value_total(:,k,l)),&
                                 stats%vid_count_total(LVT_MEANid,2),k)                         
                         enddo
                      endif
                   
                      if(LVT_metrics%mean%computeSC.eq.1) then 
                         do tind = 1,LVT_rc%nasc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%mean%ftn_total,&
                                 stats%mean%model_value_asc(:,k,tind),&
                                 stats%vid_sc_total(tind,LVT_MEANid,1),k)
                         enddo
                         if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                            do tind = 1,LVT_rc%nasc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%mean%ftn_total,&
                                    stats%mean%obs_value_asc(:,k,tind),&
                                    stats%vid_sc_total(tind,LVT_MEANid,2),k)
                            enddo
                         endif
                      endif
                      if(LVT_metrics%mean%computeADC.eq.1) then 
                         do tind = 1,LVT_rc%nadc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%mean%ftn_total,&
                                 stats%mean%model_value_adc(:,k,tind),&
                                 stats%vid_adc_total(tind,LVT_MEANid,1),k)
                         enddo
                         if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                            do tind = 1,LVT_rc%nadc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%mean%ftn_total,&
                                    stats%mean%obs_value_adc(:,k,tind),&
                                    stats%vid_adc_total(tind,LVT_MEANid,2),k)
                            enddo
                         endif
                      endif
                      call LVT_writeSummaryStats(&
                           LVT_metrics%mean%ftn_summ,&
                           LVT_metrics%mean%short_name,&
                           LVT_rc%ngrid,&
                           stats%mean%model_value_total(:,k,:), &
                           stats%mean%count_model_value_total(:,k,:),&
                           stats%standard_name,&
                           stats%mean%model_value_ci(k,:))
                      if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                         call LVT_writeSummaryStats(&
                              LVT_metrics%mean%ftn_summ,&
                              LVT_metrics%mean%short_name,&
                              LVT_rc%ngrid,&
                              stats%mean%obs_value_total(:,k,:), &
                              stats%mean%count_obs_value_total(:,k,:),&
                              "DS2_"//trim(stats%standard_name),&
                              stats%mean%obs_value_ci(k,:))
                      endif
                   endif
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_writeMetric_MEAN

!BOP
! 
! !ROUTINE: LVT_resetMetric_MEAN
! \label(LVT_resetMetric_MEAN)
!
! !INTERFACE:
  subroutine LVT_resetMetric_MEAN(alarm)
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
             if(LVT_metrics%mean%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   stats%mean%model_value_ts(:,k,l) = 0.0
                   stats%mean%count_model_value_ts(:,k,l)=0 
                   if(alarm) then 
                      stats%mean%tavg_model_value_ts(:,k,l) = 0.0
                      stats%mean%tavg_count_model_value_ts(:,k,l)=0 
                   endif
                enddo
                if(LVT_rc%obssource(2).ne."none".and.&
                     obs%selectNlevs.ge.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      stats%mean%obs_value_ts(:,k,l) = 0.0
                      stats%mean%count_obs_value_ts(:,k,l)=0 
                      if(alarm) then 
                         stats%mean%tavg_obs_value_ts(:,k,l) = 0.0
                         stats%mean%tavg_count_obs_value_ts(:,k,l)=0 
                      endif
                   enddo
                   if(obs%stdev_flag) then 
                      do l=1,LVT_rc%strat_nlevels
                         stats%mean%obs_value_stdev_ts(:,k,l) = 0.0
                         stats%mean%count_obs_value_stdev_ts(:,k,l)=0 
                         if(alarm) then 
                            stats%mean%tavg_obs_value_stdev_ts(:,k,l) = 0.0
                            stats%mean%tavg_count_obs_value_stdev_ts(:,k,l)=0 
                         endif
                      enddo
                      
                   endif
                endif
             endif
          enddo
       endif
       
       model => model%next
       obs => obs%next
       stats => stats%next
    enddo

  end subroutine LVT_resetMetric_MEAN

!BOP
! 
! !ROUTINE: LVT_writerestart_MEAN
! 
! !INTERFACE:
  subroutine LVT_writerestart_MEAN(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for MEAN metric computations
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
       if(LVT_metrics%mean%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   call LVT_writevar_restart(ftn,&
                        stats%mean%model_value_total(:,k,l))
                   call LVT_writevar_restart(ftn,&
                        stats%mean%count_model_value_total(:,k,l))
                enddo
             enddo
             if(LVT_metrics%mean%computeSC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc         
                      call LVT_writevar_restart(ftn,&
                           stats%mean%model_value_asc(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%mean%count_model_value_asc(:,k,l))
                   enddo
                enddo
             endif
             if(LVT_metrics%mean%computeADC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nadc         
                      call LVT_writevar_restart(ftn,&
                           stats%mean%model_value_adc(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%mean%count_model_value_adc(:,k,l))
                   enddo
                enddo
             endif
             
             if(LVT_rc%obssource(2).ne."none") then 
                if(obs%selectNlevs.ge.1) then 
                   do k=1,obs%selectNlevs
                      do l=1,LVT_rc%strat_nlevels
                         call LVT_writevar_restart(ftn,&
                              stats%mean%obs_value_total(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%mean%count_obs_value_total(:,k,l))
                      enddo
                   enddo
                   
                   if(LVT_metrics%mean%computeSC.eq.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%nasc
                            call LVT_writevar_restart(ftn,&
                                 stats%mean%obs_value_asc(:,k,l))
                            call LVT_writevar_restart(ftn,&
                                 stats%mean%count_obs_value_asc(:,k,l))
                         enddo
                      enddo
                   endif
                   if(LVT_metrics%mean%computeADC.eq.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%nasc
                            call LVT_writevar_restart(ftn,&
                                 stats%mean%obs_value_adc(:,k,l))
                            call LVT_writevar_restart(ftn,&
                                 stats%mean%count_obs_value_adc(:,k,l))
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
  end subroutine LVT_writerestart_MEAN


!BOP
! 
! !ROUTINE: LVT_readrestart_MEAN
! 
! !INTERFACE:
  subroutine LVT_readrestart_MEAN(ftn)
! !USES: 
! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for MEAN metric computations
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
       if(LVT_metrics%mean%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 

             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   call LVT_readvar_restart(ftn,&
                        stats%mean%model_value_total(:,k,l))
                   call LVT_readvar_restart(ftn,&
                        stats%mean%count_model_value_total(:,k,l))
                enddo
             enddo
             if(LVT_metrics%mean%computeSC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc         
                      call LVT_readvar_restart(ftn,&
                           stats%mean%model_value_asc(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%mean%count_model_value_asc(:,k,l))
                   enddo
                enddo
             endif
             if(LVT_metrics%mean%computeADC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nadc         
                      call LVT_readvar_restart(ftn,&
                           stats%mean%model_value_adc(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%mean%count_model_value_adc(:,k,l))
                   enddo
                enddo
             endif
             
             if(LVT_rc%obssource(2).ne."none") then 
                if(obs%selectNlevs.ge.1) then 
                   do k=1,obs%selectNlevs
                      do l=1,LVT_rc%strat_nlevels
                         call LVT_readvar_restart(ftn,&
                              stats%mean%obs_value_total(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%mean%count_obs_value_total(:,k,l))
                      enddo
                   enddo
                   
                   if(LVT_metrics%mean%computeSC.eq.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%nasc
                            call LVT_readvar_restart(ftn,&
                                 stats%mean%obs_value_asc(:,k,l))
                            call LVT_readvar_restart(ftn,&
                                 stats%mean%count_obs_value_asc(:,k,l))
                         enddo
                      enddo
                   endif
                   if(LVT_metrics%mean%computeADC.eq.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%nasc
                            call LVT_readvar_restart(ftn,&
                                 stats%mean%obs_value_adc(:,k,l))
                            call LVT_readvar_restart(ftn,&
                                 stats%mean%count_obs_value_adc(:,k,l))
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
  end subroutine LVT_readrestart_MEAN


end module LVT_MEANMod
