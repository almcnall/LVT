!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_ensMEANMod
! \label(LVT_ensMEANMod)
!
! !INTERFACE:
module LVT_ensMEANMod
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
!  !DESCRIPTION: 
!   This module handles the computations required to compute 
!   ensemble mean values
!   of desired variables from the LIS output. Note that the LIS
!   output is expected to be written in tile space format. 
! 
! !FILES USED:
!
!  !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initensMEAN
  public :: LVT_diagnoseensMEAN
  public :: LVT_computeensMEAN
  public :: LVT_writeMetric_ensMEAN
  public :: LVT_resetMetric_ensMEAN
  public :: LVT_writerestart_ensMEAN
  public :: LVT_readrestart_ensMEAN
!EOP
  
  private

contains
  subroutine LVT_initensMEAN(model, obs, stats,metric)
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
       allocate(stats%ensmean%model_value_total(LVT_LIS_rc(1)%ntiles, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensmean%model_value_total = 0.0
       allocate(stats%ensmean%count_model_value_total(LVT_LIS_rc(1)%ntiles, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensmean%count_model_value_total = 0
       allocate(stats%ensmean%model_value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%ensmean%model_value_ci = LVT_rc%udef
       
       if(metric%computeSC.eq.1) then 
          allocate(stats%ensmean%model_value_asc(LVT_LIS_rc(1)%ntiles, model%selectNlevs,&
               LVT_rc%nasc))
          stats%ensmean%model_value_asc = 0.0
          allocate(stats%ensmean%count_model_value_asc(LVT_LIS_rc(1)%ntiles, model%selectNlevs,&
               LVT_rc%nasc))
          stats%ensmean%count_model_value_asc = 0
       endif
       if(metric%computeADC.eq.1) then 
          allocate(stats%ensmean%model_value_adc(LVT_LIS_rc(1)%ntiles, obs%selectNlevs,&
               LVT_rc%nadc))
          stats%ensmean%model_value_adc = 0.0
          allocate(stats%ensmean%count_model_value_adc(LVT_LIS_rc(1)%ntiles, obs%selectNlevs,&
               LVT_rc%nadc))
          stats%ensmean%count_model_value_adc = 0
       endif
       
       if(trim(LVT_rc%obssource(2)).ne."none") then 
          if(obs%selectNlevs.ge.1) then 
             allocate(stats%ensmean%obs_value_total(LVT_rc%ngrid, obs%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%ensmean%obs_value_total = 0.0
             allocate(stats%ensmean%count_obs_value_total(LVT_rc%ngrid, obs%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%ensmean%count_obs_value_total = 0
             allocate(stats%ensmean%obs_value_ci(obs%selectNlevs,LVT_rc%strat_nlevels))
             stats%ensmean%obs_value_ci = LVT_rc%udef
             
             if(metric%computeSC.eq.1) then 
                allocate(stats%ensmean%obs_value_asc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nasc))
                stats%ensmean%obs_value_asc = 0.0
                allocate(stats%ensmean%count_obs_value_asc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nasc))
                stats%ensmean%count_obs_value_asc = 0
             endif
             if(metric%computeADC.eq.1) then 
                allocate(stats%ensmean%obs_value_adc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nadc))
                stats%ensmean%obs_value_adc = 0.0
                allocate(stats%ensmean%count_obs_value_adc(LVT_rc%ngrid, obs%selectNlevs,&
                     LVT_rc%nadc))
                stats%ensmean%count_obs_value_adc = 0
             endif
          endif
       endif
    endif
    if(metric%timeOpt.eq.1) then 
       allocate(stats%ensmean%model_value_ts(LVT_LIS_rc(1)%ntiles, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensmean%model_value_ts = 0.0
       allocate(stats%ensmean%count_model_value_ts(LVT_LIS_rc(1)%ntiles, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensmean%count_model_value_ts = 0

       allocate(stats%ensmean%tavg_model_value_ts(LVT_LIS_rc(1)%ntiles, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensmean%tavg_model_value_ts = 0.0
       allocate(stats%ensmean%tavg_count_model_value_ts(LVT_LIS_rc(1)%ntiles, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensmean%tavg_count_model_value_ts = 0

       
       if(trim(LVT_rc%obssource(2)).ne."none") then 
          if(obs%selectNlevs.ge.1) then 
             allocate(stats%ensmean%obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%ensmean%obs_value_ts = 0.0
             allocate(stats%ensmean%count_obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%ensmean%count_obs_value_ts = 0

             allocate(stats%ensmean%tavg_obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%ensmean%tavg_obs_value_ts = 0.0
             allocate(stats%ensmean%tavg_count_obs_value_ts(LVT_rc%ngrid, model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%ensmean%tavg_count_obs_value_ts = 0
          endif
       endif
    endif

!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1  
    if(obs%selectNlevs.ge.1) then   
       metric%obsData = .true. 
    else
       metric%obsData = .false. 
    endif
    metric%stdevFlag = .false. 

  end subroutine LVT_initensMEAN

!BOP
! 
! !ROUTINE: LVT_diagnoseensMEAN
! \label{LVT_diagnoseensMEAN}
!
! !INTERFACE: 
  subroutine LVT_diagnoseensMEAN(pass)
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
!    \item[diagnoseSingleModelensMEAN](\ref{diagnoseSingleModelensMEAN})
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
       if(LVT_metrics%ensmean%selectOpt.eq.1.or.&
            LVT_metrics%ensmean%timeOpt.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call diagnoseSingleensMEAN(&
                  model, obs, stats, &
                  LVT_metrics%ensmean)
             
             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseensMEAN

!BOP
! 
! !ROUTINE: diagnoseSingleensMEAN
! \label{diagnoseSingleensMEAN}
!
! !INTERFACE: 
  subroutine diagnoseSingleensMEAN(model, obs, stats,metric)
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
! 
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: model
    type(LVT_metadataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
!EOP
    integer    :: t,k,tind,gid

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then        
       do t=1,LVT_LIS_rc(1)%ntiles
          gid = LVT_domain%gindex(LVT_LIS_domain(1)%tile(t)%col,&
               LVT_LIS_domain(1)%tile(t)%row)
          do k=1,model%selectNlevs
             if(model%count(t,k).gt.0) then
                if(metric%selectOpt.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      stats%ensmean%model_value_total(t,k,1) = &
                           stats%ensmean%model_value_total(t,k,1) + &
                           model%value(t,k)
                      stats%ensmean%count_model_value_total(t,k,1) = &
                           stats%ensmean%count_model_value_total(t,k,1) + 1
                   endif
                endif
                
                if(metric%timeOpt.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      stats%ensmean%model_value_ts(t,k,1) = &
                           stats%ensmean%model_value_ts(t,k,1)+&
                           model%value(t,k)
                      stats%ensmean%count_model_value_ts(t,k,1) = & 
                           stats%ensmean%count_model_value_ts(t,k,1)+1
                   endif
                endif
                if(metric%computeSC.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                           tind)
                      stats%ensmean%model_value_asc(t,k,tind) = &
                           stats%ensmean%model_value_asc(t,k,tind)+&
                           model%value(t,k)
                      stats%ensmean%count_model_value_asc(t,k,tind) = &
                           stats%ensmean%count_model_value_asc(t,k,tind) + 1
                   endif
                endif
                if(metric%computeADC.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      call LVT_getADCTimeIndex(tind)
                      stats%ensmean%model_value_adc(t,k,tind) = &
                           stats%ensmean%model_value_adc(t,k,tind)+&
                           model%value(t,k)
                      stats%ensmean%count_model_value_adc(t,k,tind) = &
                           stats%ensmean%count_model_value_adc(t,k,tind) + 1
                   endif
                endif
                if(LVT_rc%strat_nlevels.gt.1) then 
                   if(LVT_stats%strat_var(gid,k).gt.&
                        LVT_rc%strat_var_threshold) then
                      if(metric%selectOpt.eq.1) then
                         if(model%value(t,k).ne.LVT_rc%udef) then  
                            stats%ensmean%model_value_total(t,k,2) = &
                                 stats%ensmean%model_value_total(t,k,2) + &
                                 model%value(t,k)
                            stats%ensmean%count_model_value_total(t,k,2) = &
                                 stats%ensmean%count_model_value_total(t,k,2) + 1
                         endif
                      endif
                      if(metric%timeOpt.eq.1) then 
                         if(model%value(t,k).ne.LVT_rc%udef) then 
                            stats%ensmean%model_value_ts(t,k,2) = &
                                 stats%ensmean%model_value_ts(t,k,2)+&
                                 model%value(t,k)
                            stats%ensmean%count_model_value_ts(t,k,2) = & 
                                 stats%ensmean%count_model_value_ts(t,k,2)+1
                         endif
                      endif
                   elseif(LVT_stats%strat_var(gid,k).le.&
                        LVT_rc%strat_var_threshold) then
                      if(metric%selectOpt.eq.1) then 
                         if(model%value(t,k).ne.LVT_rc%udef) then 
                            stats%ensmean%model_value_total(t,k,3) = &
                                 stats%ensmean%model_value_total(t,k,3) + &
                                 model%value(t,k)
                            stats%ensmean%count_model_value_total(t,k,3) = &
                                 stats%ensmean%count_model_value_total(t,k,3) + 1
                         endif
                      endif
                      if(metric%timeOpt.eq.1) then 
                         if(model%value(t,k).ne.LVT_rc%udef) then 
                            stats%ensmean%model_value_ts(t,k,3) = &
                                 stats%ensmean%model_value_ts(t,k,3)+&
                                 model%value(t,k)
                            stats%ensmean%count_model_value_ts(t,k,3) = & 
                                 stats%ensmean%count_model_value_ts(t,k,3)+1
                         endif
                      endif
                   endif
                endif
             endif
          enddo
       enddo
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             if(trim(LVT_rc%obssource(2)).ne."none") then
                if(obs%selectNlevs.ge.1) then
                   if(obs%count(t,k).gt.0) then 
                      if(metric%selectOpt.eq.1.and.&
                           obs%value(t,k).ne.LVT_rc%udef) then 
                         stats%ensmean%obs_value_total(t,k,1) = &
                              stats%ensmean%obs_value_total(t,k,1) + &
                              obs%value(t,k)
                         stats%ensmean%count_obs_value_total(t,k,1) = &
                              stats%ensmean%count_obs_value_total(t,k,1) + 1
                      endif
                      
                      if(metric%timeOpt.eq.1.and.obs%value(t,k).ne.LVT_rc%udef) then 
                         stats%ensmean%obs_value_ts(t,k,1) = stats%ensmean%obs_value_ts(t,k,1)+&
                              obs%value(t,k)
                         stats%ensmean%count_obs_value_ts(t,k,1) = & 
                              stats%ensmean%count_obs_value_ts(t,k,1)+1
                      endif
                      if(metric%computeSC.eq.1.and.obs%value(t,k).ne.LVT_rc%udef) then 
                         call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,tind)
                         stats%ensmean%obs_value_asc(t,k,tind) = &
                              stats%ensmean%obs_value_asc(t,k,tind)+&
                              obs%value(t,k)
                         stats%ensmean%count_obs_value_asc(t,k,tind) = &
                              stats%ensmean%count_obs_value_asc(t,k,tind) + 1
                      endif
                      if(metric%computeADC.eq.1.and.obs%value(t,k).ne.LVT_rc%udef) then 
                         call LVT_getADCTimeIndex(tind)
                         stats%ensmean%obs_value_adc(t,k,tind) = &
                              stats%ensmean%obs_value_adc(t,k,tind)+&
                              obs%value(t,k)
                         stats%ensmean%count_obs_value_adc(t,k,tind) = &
                              stats%ensmean%count_obs_value_adc(t,k,tind) + 1
                      endif
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(gid,k).gt.&
                              LVT_rc%strat_var_threshold) then
                            if(metric%selectOpt.eq.1 &
                                 .and.obs%value(t,k).ne.LVT_rc%udef) then 
                               stats%ensmean%obs_value_total(t,k,2) = &
                                    stats%ensmean%obs_value_total(t,k,2) + &
                                    obs%value(t,k)
                               stats%ensmean%count_obs_value_total(t,k,2) = &
                                    stats%ensmean%count_obs_value_total(t,k,2) + 1
                            endif
                            
                            if(metric%timeOpt.eq.1.and.&
                                 obs%value(t,k).ne.LVT_rc%udef) then 
                               stats%ensmean%obs_value_ts(t,k,2) = &
                                    stats%ensmean%obs_value_ts(t,k,2)+&
                                    obs%value(t,k)
                               stats%ensmean%count_obs_value_ts(t,k,2) = & 
                                    stats%ensmean%count_obs_value_ts(t,k,2)+1
                            endif
                         elseif(LVT_stats%strat_var(gid,k).le.&
                              LVT_rc%strat_var_threshold) then
                            if(metric%selectOpt.eq.1.and.&
                                 obs%value(t,k).ne.LVT_rc%udef) then 
                               stats%ensmean%obs_value_total(t,k,3) = &
                                    stats%ensmean%obs_value_total(t,k,3) + &
                                    obs%value(t,k)
                               stats%ensmean%count_obs_value_total(t,k,3) = &
                                    stats%ensmean%count_obs_value_total(t,k,3) + 1
                            endif
                            
                            if(metric%timeOpt.eq.1.and.&
                                 obs%value(t,k).ne.LVT_rc%udef) then 
                               stats%ensmean%obs_value_ts(t,k,3) = &
                                    stats%ensmean%obs_value_ts(t,k,3)+&
                                    obs%value(t,k)
                               stats%ensmean%count_obs_value_ts(t,k,3) = & 
                                    stats%ensmean%count_obs_value_ts(t,k,3)+1
                            endif
                            
                         endif
                      endif
                   endif
                endif
             endif
          enddo
       enddo
    endif
  end subroutine diagnoseSingleensMEAN

!BOP
! 
! !ROUTINE: LVT_computeensMEAN
! \label{LVT_computeensMEAN}
!
! !INTERFACE: 
  subroutine LVT_computeensMEAN(pass,alarm)
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
!    \item[computeSingleModelensMEAN](\ref{computeSingleModelensMEAN})
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
       if(LVT_metrics%ensmean%selectOpt.eq.1.or.&
            LVT_metrics%ensmean%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%ensmean%timeOpt.eq.1.and.&
                  LVT_metrics%ensmean%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ensmean%ftn_ts_loc(i),200,advance='no') &
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
             call computeSingleensMEAN(alarm,&
                  model, obs, stats, &
                  LVT_metrics%ensmean)

             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%ensmean%timeOpt.eq.1.and.&
                  LVT_metrics%ensmean%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ensmean%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_computeensMEAN
  

!BOP
! 
! !ROUTINE: computeSingleensMEAN
! \label{computeSingleensMEAN}
!
! !INTERFACE: 
  subroutine computeSingleensMEAN(alarm,model,obs,stats,metric)
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
          do t=1,LVT_LIS_rc(1)%ntiles
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%ensmean%count_model_value_ts(t,k,l).gt.0) then 
                      stats%ensmean%model_value_ts(t,k,l) = stats%ensmean%model_value_ts(t,k,l)/&
                           stats%ensmean%count_model_value_ts(t,k,l)                   
                      
                      stats%ensmean%tavg_model_value_ts(t,k,l) = & 
                           stats%ensmean%tavg_model_value_ts(t,k,l)  + & 
                           stats%ensmean%model_value_ts(t,k,l)
                      stats%ensmean%tavg_count_model_value_ts(t,k,l) = & 
                           stats%ensmean%tavg_count_model_value_ts(t,k,l)  + 1
                   else
                      stats%ensmean%model_value_ts(t,k,l) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(trim(LVT_rc%obssource(2)).ne."none") then 
                      if(obs%selectNlevs.ge.1) then 
                         if(stats%ensmean%count_obs_value_ts(t,k,l).gt.0) then 
                            stats%ensmean%obs_value_ts(t,k,l) = &
                                 stats%ensmean%obs_value_ts(t,k,l)/&
                                 stats%ensmean%count_obs_value_ts(t,k,l)                   
                         else
                            stats%ensmean%obs_value_ts(t,k,l) = LVT_rc%udef
                         endif
                      endif
                   endif
                enddo
             enddo
          enddo

          if(alarm) then 

             do t=1,LVT_LIS_rc(1)%ntiles
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%ensmean%tavg_count_model_value_ts(t,k,l).gt.0) then 
                         stats%ensmean%tavg_model_value_ts(t,k,l) = & 
                              stats%ensmean%tavg_model_value_ts(t,k,l)/ & 
                              stats%ensmean%tavg_count_model_value_ts(t,k,l)
                      endif
                   enddo
                enddo
             enddo

             if(metric%extractTS.eq.1) then 
                if(trim(LVT_rc%obssource(2)).ne."none".and.obs%selectNlevs.ge.1) then 
                   call LVT_writeTSinfo(metric%ftn_ts_loc,&
                        model,&
                        LVT_LIS_rc(1)%ntiles,&
                        stats%ensmean%tavg_model_value_ts,&
                        stats%ensmean%tavg_count_model_value_ts,&
                        LVT_rc%ngrid,&
                        stats%ensmean%tavg_obs_value_ts,&
                        stats%ensmean%tavg_count_obs_value_ts)
                   
                else
                   call LVT_writeTSinfo(metric%ftn_ts_loc,&
                        model,&
                        LVT_LIS_rc(1)%ntiles,&
                        stats%ensmean%tavg_model_value_ts,&
                        stats%ensmean%tavg_count_model_value_ts)
                endif
             endif
          endif
       endif
    endif

       
    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_LIS_rc(1)%ntiles
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%ensmean%count_model_value_total(t,k,l).gt.&
                        LVT_rc%obsCountThreshold) then 
                      stats%ensmean%model_value_total(t,k,l) = &
                           stats%ensmean%model_value_total(t,k,l)/&
                              stats%ensmean%count_model_value_total(t,k,l)           
                   else
                      stats%ensmean%model_value_total(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1, LVT_rc%nasc
                      if(stats%ensmean%count_model_value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then 
                         stats%ensmean%model_value_asc(t,k,l) = &
                              stats%ensmean%model_value_asc(t,k,l)/&
                              stats%ensmean%count_model_value_asc(t,k,l)           
                      else
                         stats%ensmean%model_value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif

                if(metric%computeADC.eq.1) then 
                   do l=1, LVT_rc%nadc
                      if(stats%ensmean%count_model_value_adc(t,k,l).gt.&
                           LVT_rc%ADCCountThreshold) then 
                         stats%ensmean%model_value_adc(t,k,l) = &
                              stats%ensmean%model_value_adc(t,k,l)/&
                              stats%ensmean%count_model_value_adc(t,k,l)           
                      else
                         stats%ensmean%model_value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
             enddo
          enddo
          if(trim(LVT_rc%obssource(2)).ne."none") then 
             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   if(obs%selectNlevs.ge.1) then
                      do l=1,LVT_rc%strat_nlevels                      
                         if(stats%ensmean%count_obs_value_total(t,k,l).gt.&
                              LVT_rc%obsCountThreshold) then 
                            stats%ensmean%obs_value_total(t,k,l) = &
                                 stats%ensmean%obs_value_total(t,k,l)/&
                                 stats%ensmean%count_obs_value_total(t,k,l)       
                         else
                            stats%ensmean%obs_value_total(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                      if(metric%computeSC.eq.1) then
                         do l=1, LVT_rc%nasc 
                            if(stats%ensmean%count_obs_value_asc(t,k,l).gt.&
                                 LVT_rc%SCCountThreshold) then 
                               stats%ensmean%obs_value_asc(t,k,l) = &
                                    stats%ensmean%obs_value_asc(t,k,l)/&
                                    stats%ensmean%count_obs_value_asc(t,k,l)           
                            else
                               stats%ensmean%obs_value_asc(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                      endif
                      
                      if(metric%computeADC.eq.1) then
                         do l=1, LVT_rc%nadc  
                            if(stats%ensmean%count_obs_value_adc(t,k,l).gt.&
                                 LVT_rc%ADCCountThreshold) then 
                               stats%ensmean%obs_value_adc(t,k,l) = &
                                    stats%ensmean%obs_value_adc(t,k,l)/&
                                    stats%ensmean%count_obs_value_adc(t,k,l)           
                            else
                               stats%ensmean%obs_value_adc(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                      endif
                   endif
                enddo
             enddo
          endif
             
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%ensmean%model_value_total(:,k,l),LVT_LIS_rc(1)%ntiles,&
                     LVT_rc%pval_CI,stats%ensmean%model_value_ci(k,l))
             enddo
          enddo
          do k=1,obs%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%ensmean%obs_value_total(:,k,l),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%ensmean%obs_value_ci(k,l))
             enddo
          enddo
          
          if(trim(LVT_rc%obssource(2)).ne."none".and.obs%selectNlevs.ge.1) then 
             call LVT_writeDataBasedStrat(model,obs,stats,metric,&
                  LVT_LIS_rc(1)%ntiles,stats%ensmean%model_value_total, &
                  LVT_rc%ngrid,stats%ensmean%obs_value_total)
             
             if(metric%computeSC.eq.1) then 
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_LIS_rc(1)%ntiles,stats%ensmean%model_value_asc,stats%ensmean%count_model_value_asc,&
                     LVT_rc%ngrid,stats%ensmean%obs_value_asc,stats%ensmean%count_obs_value_asc)          
             endif
             if(metric%computeADC.eq.1) then 
                call LVT_writeAvgDiurnalCycleInfo(model,obs,stats,metric,&
                     LVT_LIS_rc(1)%ntiles,stats%ensmean%model_value_adc,stats%ensmean%count_model_value_adc,&
                     LVT_rc%ngrid,stats%ensmean%obs_value_adc,stats%ensmean%count_obs_value_adc)     
             endif
          else
             call LVT_writeDataBasedStrat(model,obs,stats,metric,&
                  LVT_LIS_rc(1)%ntiles,stats%ensmean%model_value_total)

             if(metric%computeSC.eq.1) then 
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_LIS_rc(1)%ntiles,stats%ensmean%model_value_asc,stats%ensmean%count_model_value_asc)
             endif
             if(metric%computeADC.eq.1) then 
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_LIS_rc(1)%ntiles,stats%ensmean%model_value_adc,stats%ensmean%count_model_value_adc)
             endif
          endif
       endif
    endif
  end subroutine computeSingleensMEAN


  subroutine LVT_writeMetric_ensMEAN(pass,final,vlevels,stats,obs)
! !USES:
    use LVT_statsMod, only : LVT_writeSummaryStats
    use LVT_pluginIndices
! !ARGUMENTS:
    implicit none

    integer                 :: pass
    integer                 :: final
    integer                 :: vlevels
    type(LVT_statsEntry)    :: stats
    type(LVT_metaDataEntry) :: obs

    real                    :: dummy
    integer                 :: k,l,tind

    if(pass.eq.LVT_metrics%ensmean%npass) then 
       if(final.ne.1) then
           if(stats%selectOpt.eq.1) then 
             do k=1,vlevels
                if(LVT_metrics%ensmean%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_gridded(LVT_metrics%ensmean%ftn_ts, &
                           stats%ensmean%model_value_ts(:,k,l),&
                           stats%vid_ts(LVT_ensMEANid,1),dummy,k)
                      call LVT_writevar_gridded(LVT_metrics%ensmean%ftn_ts, &
                           real(stats%ensmean%count_model_value_ts(:,k,l)),&
                           stats%vid_count_ts(LVT_ensMEANid,1),dummy,k)
                   enddo
                   if(trim(LVT_rc%obssource(2)).ne."none") then 
                      do l=1,LVT_rc%strat_nlevels
                         call LVT_writevar_gridded(LVT_metrics%ensmean%ftn_ts, &
                              stats%ensmean%obs_value_ts(:,k,l),&
                              stats%vid_ts(LVT_ensMEANid,2),k)
                         call LVT_writevar_gridded(LVT_metrics%ensmean%ftn_ts, &
                              real(stats%ensmean%count_obs_value_ts(:,k,l)),&
                              stats%vid_count_ts(LVT_ensMEANid,2),k)
                      enddo
                   endif
                endif
             enddo
          endif
       else
          if(pass.eq.1) then 
             if(stats%selectOpt.eq.1) then 
                do k=1,vlevels
                   if(LVT_metrics%ensmean%selectOpt.eq.1) then 
                      do l=1,LVT_rc%strat_nlevels                      
                         call LVT_writevar_gridded(LVT_metrics%ensmean%ftn_total, &
                              stats%ensmean%model_value_total(:,k,l),&
                              stats%vid_total(LVT_ensMEANid,1),dummy,k)
                         call LVT_writevar_gridded(LVT_metrics%ensmean%ftn_total, &
                              real(stats%ensmean%count_model_value_total(:,k,l)),&
                              stats%vid_count_total(LVT_ensMEANid,1),dummy,k)
                         
                      enddo
                      if(trim(LVT_rc%obssource(2)).ne."none".and.&
                           obs%selectNlevs.ge.1) then 
                         do l=1,LVT_rc%strat_nlevels                      
                            call LVT_writevar_gridded(&
                                 LVT_metrics%ensmean%ftn_total, &
                                 stats%ensmean%obs_value_total(:,k,l),&
                                 stats%vid_total(LVT_ensMEANid,2),k)
                            call LVT_writevar_gridded(&
                                 LVT_metrics%ensmean%ftn_total, &
                                 real(stats%ensmean%count_obs_value_total(:,k,l)),&
                                 stats%vid_count_total(LVT_ensMEANid,2),k)                         
                         enddo
                      endif
                   
                      if(LVT_metrics%ensmean%computeSC.eq.1) then 
                         do tind = 1,LVT_rc%nasc
                            call LVT_writevar_gridded(LVT_metrics%ensmean%ftn_total,&
                                 stats%ensmean%model_value_asc(:,k,tind),&
                                 stats%vid_sc_total(tind,LVT_ensMEANid,1),dummy,1)
                         enddo
                         if(trim(LVT_rc%obssource(2)).ne."none".and.obs%selectNlevs.ge.1) then 
                            do tind = 1,LVT_rc%nasc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%ensmean%ftn_total,&
                                    stats%ensmean%obs_value_asc(:,k,tind),&
                                    stats%vid_sc_total(tind,LVT_ensMEANid,2),&
                                    1)
                            enddo
                         endif
                      endif
                      if(LVT_metrics%ensmean%computeADC.eq.1) then 
                         do tind = 1,LVT_rc%nadc
                            call LVT_writevar_gridded(LVT_metrics%ensmean%ftn_total,&
                                 stats%ensmean%model_value_adc(:,k,tind),&
                                 stats%vid_adc_total(tind,LVT_ensMEANid,1),&
                                 dummy,1)
                         enddo
                         if(trim(LVT_rc%obssource(2)).ne."none".and.obs%selectNlevs.ge.1) then 
                            do tind = 1,LVT_rc%nadc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%ensmean%ftn_total,&
                                    stats%ensmean%obs_value_adc(:,k,tind),&
                                    stats%vid_adc_total(tind,LVT_ensMEANid,2),&
                                    1)
                            enddo
                         endif
                      endif
                      call LVT_writeSummaryStats(&
                           LVT_metrics%ensmean%ftn_summ,&
                           LVT_metrics%ensmean%short_name,&
                           LVT_LIS_rc(1)%ntiles,&
                           stats%ensmean%model_value_total(:,k,:), &
                           stats%ensmean%count_model_value_total(:,k,:),&
                           stats%standard_name,&
                           stats%ensmean%model_value_ci(k,:))
                      if(trim(LVT_rc%obssource(2)).ne."none".and.obs%selectNlevs.ge.1) then 
                         call LVT_writeSummaryStats(&
                              LVT_metrics%ensmean%ftn_summ,&
                              LVT_metrics%ensmean%short_name,&
                              LVT_rc%ngrid,&
                              stats%ensmean%obs_value_total(:,k,:), &
                              stats%ensmean%count_obs_value_total(:,k,:),&
                              "OBS_"//trim(stats%standard_name),&
                              stats%ensmean%obs_value_ci(k,:))
                      endif
                   endif
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_writeMetric_ensMEAN

  subroutine LVT_resetMetric_ensMEAN(alarm)
    
    logical                :: alarm

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
             if(LVT_metrics%ensmean%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   stats%ensmean%model_value_ts(:,k,l) = 0.0
                   stats%ensmean%count_model_value_ts(:,k,l)=0 

                   if(alarm) then 
                      stats%ensmean%tavg_model_value_ts(:,k,l) = 0.0
                      stats%ensmean%tavg_count_model_value_ts(:,k,l)=0 
                   endif
                enddo
                if(trim(LVT_rc%obssource(2)).ne."none".and.&
                     obs%selectNlevs.ge.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      stats%ensmean%obs_value_ts(:,k,l) = 0.0
                      stats%ensmean%count_obs_value_ts(:,k,l)=0 

                      if(alarm) then 
                         stats%ensmean%tavg_obs_value_ts(:,k,l) = 0.0
                         stats%ensmean%tavg_count_obs_value_ts(:,k,l)=0 
                      endif
                   enddo
                endif
             endif
          enddo
       endif

       model => model%next
       obs   => obs%next
       stats => stats%next
    enddo

  end subroutine LVT_resetMetric_ensMEAN


!BOP
! 
! !ROUTINE: LVT_writerestart_ensMEAN
! 
! !INTERFACE:
  subroutine LVT_writerestart_ensMEAN(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for ensMEAN metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ensmean%selectOpt.eq.1) then 
       
       print*, 'The writerestart method is not implemented for ensMEAN'

    end if
    
  end subroutine LVT_writerestart_ensMEAN

!BOP
! 
! !ROUTINE: LVT_readrestart_ensMEAN
! 
! !INTERFACE:
  subroutine LVT_readrestart_ensMEAN(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for ensMEAN metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ensmean%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for ensMEAN'
       stop
    end if
    
  end subroutine LVT_readrestart_ensMEAN

end module LVT_ensMEANMod
