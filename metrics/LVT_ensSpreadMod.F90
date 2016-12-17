!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_ensSpreadMod
! \label(LVT_ensSpreadMod)
!
! !INTERFACE:
module LVT_ensSpreadMod
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
!   ensemble spread (defined as the difference between max and 
!   min values within the ensemble) of desired variables from the LIS output. 
!   
!   NOTES: 
!   * The LIS output should be written in a tile space format
! 
! !FILES USED:
!
!  !REVISION HISTORY: 
!  28 Feb 2012    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initensSpread
  public :: LVT_diagnoseensSpread
  public :: LVT_computeensSpread
  public :: LVT_writeMetric_ensSpread
  public :: LVT_resetMetric_ensSpread
  public :: LVT_writerestart_ensSpread
  public :: LVT_readrestart_ensSpread
  
!EOP
  
  private

contains
  subroutine LVT_initensSpread(model, obs, stats,metric)
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
       allocate(stats%ensspread%model_value_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensspread%model_value_total = 0.0
       allocate(stats%ensspread%count_model_value_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensspread%count_model_value_total = 0
       allocate(stats%ensspread%model_value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%ensspread%model_value_ci = LVT_rc%udef
       
       if(metric%computeSC.eq.1) then 
          allocate(stats%ensspread%model_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%ensspread%model_value_asc = 0.0
          allocate(stats%ensspread%count_model_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%ensspread%count_model_value_asc = 0
       endif
       if(metric%computeADC.eq.1) then 
          allocate(stats%ensspread%model_value_adc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nadc))
          stats%ensspread%model_value_adc = 0.0
          allocate(stats%ensspread%count_model_value_adc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nadc))
          stats%ensspread%count_model_value_adc = 0
       endif
    endif

    if(metric%timeOpt.eq.1) then 
       allocate(stats%ensspread%model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensspread%model_value_ts = 0.0
       allocate(stats%ensspread%count_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensspread%count_model_value_ts = 0

       allocate(stats%ensspread%tavg_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensspread%tavg_model_value_ts = 0.0
       allocate(stats%ensspread%tavg_count_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensspread%tavg_count_model_value_ts = 0
       
    endif

!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1    
    metric%obsData = .false. 
    metric%stdevFlag = .false. 

  end subroutine LVT_initensSpread

!BOP
! 
! !ROUTINE: LVT_diagnoseensSpread
! \label{LVT_diagnoseensSpread}
!
! !INTERFACE: 
  subroutine LVT_diagnoseensSpread(pass)
! 
! !USES:     

    implicit none

    integer       :: pass
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the computations for 
!   calculating the std of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelensSpread](\ref{diagnoseSingleModelensSpread})
!     updates the std computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%ensspread%selectOpt.eq.1.or.&
            LVT_metrics%ensspread%timeOpt.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleensSpread(model, obs, stats, &
                  LVT_metrics%ensspread)
             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseensSpread

!BOP
! 
! !ROUTINE: diagnoseSingleensSpread
! \label{diagnoseSingleensSpread}
!
! !INTERFACE: 
  subroutine diagnoseSingleensSpread(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the std computation of the 
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
!
!EOP
    integer    :: t,k,tind,g,m
    real       :: smax,smin,spread
    integer    :: nval

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then        
       do g=1,LVT_rc%ngrid
          do k=1, model%selectNlevs
             smax = -1E10 
             smin = 1E10
             do m=1,LVT_LIS_rc(1)%nensem
                t = (g-1)*LVT_LIS_rc(1)%nensem+m                
                if(model%count(t,k).gt.0) then 
                   if(metric%selectOpt.eq.1) then 
                      if(model%value(t,k).ne.LVT_rc%udef) then 
                         if(model%value(t,k).gt.smax) then 
                            smax = model%value(t,k)
                         endif
                         if(model%value(t,k).lt.smin) then 
                            smin = model%value(t,k)
                         endif
                      endif
                   endif
                else
                   smax = LVT_rc%udef
                   smin = 0 
                endif
             enddo
             spread = smax - smin
             stats%ensspread%model_value_total(g,k,1) = & 
                  stats%ensspread%model_value_total(g,k,1) + spread
             stats%ensspread%count_model_value_total(g,k,1) = & 
                  stats%ensspread%count_model_value_total(g,k,1) + 1
                   
             if(metric%timeOpt.eq.1) then 
                stats%ensspread%model_value_ts(g,k,1) = & 
                     stats%ensspread%model_value_ts(g,k,1) + spread
                stats%ensspread%count_model_value_ts(g,k,1) = &
                     stats%ensspread%count_model_value_ts(g,k,1) + 1
                
             endif
             if(metric%computeSC.eq.1) then 
                call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                     tind)      
                stats%ensspread%model_value_asc(g,k,tind) = & 
                     stats%ensspread%model_value_asc(g,k,tind) + spread
                stats%ensspread%count_model_value_asc(g,k,tind) = & 
                     stats%ensspread%count_model_value_asc(g,k,tind) + 1
             endif
             if(metric%computeADC.eq.1) then 
                call LVT_getADCTimeIndex(tind)      
                stats%ensspread%model_value_adc(g,k,tind) = & 
                     stats%ensspread%model_value_adc(g,k,tind) + spread
                stats%ensspread%count_model_value_adc(g,k,tind) = & 
                     stats%ensspread%count_model_value_adc(g,k,tind) + 1
             endif
             if(LVT_rc%strat_nlevels.gt.1) then 
                if(LVT_stats%strat_var(g,k).gt.&
                     LVT_rc%strat_var_threshold) then
                   stats%ensspread%model_value_total(g,k,2) = &
                        stats%ensspread%model_value_total(g,k,2) + spread
                   stats%ensspread%count_model_value_total(g,k,2) = &
                        stats%ensspread%count_model_value_total(g,k,2) + 1
                   
                   if(metric%timeOpt.eq.1) then 
                      stats%ensspread%model_value_ts(g,k,2) = & 
                           stats%ensspread%model_value_ts(g,k,2) + spread
                      stats%ensspread%count_model_value_ts(g,k,2) = &
                           stats%ensspread%count_model_value_ts(g,k,2) + 1
                      
                   endif
                elseif(LVT_stats%strat_var(g,k).le.&
                     LVT_rc%strat_var_threshold) then
                   stats%ensspread%model_value_total(g,k,3) = &
                        stats%ensspread%model_value_total(g,k,3) + spread
                   stats%ensspread%count_model_value_total(g,k,3) = &
                        stats%ensspread%count_model_value_total(g,k,3) + 1
                   
                   if(metric%timeOpt.eq.1) then 
                      stats%ensspread%model_value_ts(g,k,3) = & 
                           stats%ensspread%model_value_ts(g,k,3) + spread
                      stats%ensspread%count_model_value_ts(g,k,3) = &
                           stats%ensspread%count_model_value_ts(g,k,3) + 1
                      
                   endif
                endif
             endif
          enddo
       enddo
    endif
  end subroutine diagnoseSingleensSpread

!BOP
! 
! !ROUTINE: LVT_computeensSpread
! \label{LVT_computeensSpread}
!
! !INTERFACE: 
  subroutine LVT_computeensSpread(pass,alarm)
! 
! !USES: 

    implicit none

    integer               :: pass
    logical               :: alarm
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the std values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelensSpread](\ref{computeSingleModelensSpread})
!     computes the std values for a single variable
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats
    integer     :: i 
    
    if(pass.eq.1) then 
       if(LVT_metrics%ensspread%selectOpt.eq.1.or.&
            LVT_metrics%ensspread%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%ensspread%timeOpt.eq.1.and.&
                  LVT_metrics%ensspread%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ensspread%ftn_ts_loc(i),200,advance='no') &
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
             call computeSingleensSpread(alarm,&
                  model, obs, stats, &
                  LVT_metrics%ensspread)
             
             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%ensspread%timeOpt.eq.1.and.&
                  LVT_metrics%ensspread%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ensspread%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_computeensSpread
  

!BOP
! 
! !ROUTINE: computeSingleensSpread
! \label{computeSingleensSpread}
!
! !INTERFACE: 
  subroutine computeSingleensSpread(alarm,model,obs,stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the std values
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
    integer  :: dummy1,dummy2

    if(metric%timeOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%ensspread%count_model_value_ts(t,k,l).gt.0) then 
                      stats%ensspread%model_value_ts(t,k,l) = &
                           stats%ensspread%model_value_ts(t,k,l)/&
                           stats%ensspread%count_model_value_ts(t,k,l)                   


                      stats%ensspread%tavg_model_value_ts(t,k,l) = &
                           stats%ensspread%tavg_model_value_ts(t,k,l) + &
                           stats%ensspread%model_value_ts(t,k,l) 
                      stats%ensspread%tavg_count_model_value_ts(t,k,l) = &
                           stats%ensspread%tavg_count_model_value_ts(t,k,l) + 1

                   else
                      stats%ensspread%model_value_ts(t,k,l) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
          
          if(alarm) then 
             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%ensspread%tavg_count_model_value_ts(t,k,l).gt.0) then 
                         stats%ensspread%tavg_model_value_ts(t,k,l) = &
                              stats%ensspread%tavg_model_value_ts(t,k,l) / &
                              stats%ensspread%tavg_count_model_value_ts(t,k,l) 
                      else
                         stats%ensspread%tavg_model_value_ts(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                enddo
             enddo

             if(metric%extractTS.eq.1) then 
                call LVT_writeTSinfo(metric%ftn_ts_loc,&
                     model,&
                     LVT_rc%ngrid,&
                     stats%ensspread%tavg_model_value_ts,&
                     stats%ensspread%tavg_count_model_value_ts,dummy1,dummy2)
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
                   if(stats%ensspread%count_model_value_total(t,k,l).gt.&
                        LVT_rc%obsCountThreshold) then 
                      stats%ensspread%model_value_total(t,k,l) = &
                           stats%ensspread%model_value_total(t,k,l)/&
                              stats%ensspread%count_model_value_total(t,k,l)           
                   else
                      stats%ensspread%model_value_total(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1, LVT_rc%nasc
                      if(stats%ensspread%count_model_value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then 
                         stats%ensspread%model_value_asc(t,k,l) = &
                              stats%ensspread%model_value_asc(t,k,l)/&
                              stats%ensspread%count_model_value_asc(t,k,l)           
                      else
                         stats%ensspread%model_value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif

                if(metric%computeADC.eq.1) then 
                   do l=1, LVT_rc%nadc
                      if(stats%ensspread%count_model_value_adc(t,k,l).gt.&
                           LVT_rc%ADCCountThreshold) then 
                         stats%ensspread%model_value_adc(t,k,l) = &
                              stats%ensspread%model_value_adc(t,k,l)/&
                              stats%ensspread%count_model_value_adc(t,k,l)           
                      else
                         stats%ensspread%model_value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
             enddo
          enddo
             
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%ensspread%model_value_total(:,k,l),&
                     LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%ensspread%model_value_ci(k,l))
             enddo
          enddo
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid,stats%ensspread%model_value_total)
          
          if(metric%computeSC.eq.1) then 
             call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%ensspread%model_value_asc,stats%ensspread%count_model_value_asc)
          endif
          if(metric%computeADC.eq.1) then 
             call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%ensspread%model_value_adc,stats%ensspread%count_model_value_adc)
          endif
       endif
    endif
  end subroutine computeSingleensSpread


  subroutine LVT_writeMetric_ensSpread(pass,final,vlevels,stats,obs)
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

    if(pass.eq.LVT_metrics%ensspread%npass) then 
       if(final.ne.1) then
           if(stats%selectOpt.eq.1) then 
             do k=1,vlevels
                if(LVT_metrics%ensspread%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_gridded(LVT_metrics%ensspread%ftn_ts, &
                           stats%ensspread%model_value_ts(:,k,l),&
                           stats%vid_ts(LVT_ensSpreadid,1),1)
                      call LVT_writevar_gridded(LVT_metrics%ensspread%ftn_ts, &
                           real(stats%ensspread%count_model_value_ts(:,k,l)),&
                           stats%vid_count_ts(LVT_ensSpreadid,1),1)
                   enddo
                endif
             enddo
          endif
       else
          if(pass.eq.1) then 
             if(stats%selectOpt.eq.1) then 
                do k=1,vlevels
                   if(LVT_metrics%ensspread%selectOpt.eq.1) then 
                      do l=1,LVT_rc%strat_nlevels                      
                         call LVT_writevar_gridded(&
                              LVT_metrics%ensspread%ftn_total, &
                              stats%ensspread%model_value_total(:,k,l),&
                              stats%vid_total(LVT_ensSpreadid,1))
                         call LVT_writevar_gridded(&
                              LVT_metrics%ensspread%ftn_total, &
                              real(stats%ensspread%count_model_value_total(:,k,l)),&
                              stats%vid_count_total(LVT_ensSpreadid,1))
                         
                      enddo
                   
                      if(LVT_metrics%ensspread%computeSC.eq.1) then 
                         do tind = 1,LVT_rc%nasc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%ensspread%ftn_total,&
                                 stats%ensspread%model_value_asc(:,k,tind),&
                                 stats%vid_sc_total(tind,LVT_ensSpreadid,1))
                         enddo
                      endif
                      if(LVT_metrics%ensspread%computeADC.eq.1) then 
                         do tind = 1,LVT_rc%nadc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%ensspread%ftn_total,&
                                 stats%ensspread%model_value_adc(:,k,tind),&
                                 stats%vid_adc_total(tind,LVT_ensSpreadid,1))
                         enddo
                      endif
                      call LVT_writeSummaryStats(&
                           LVT_metrics%ensspread%ftn_summ,&
                           LVT_metrics%ensspread%short_name,&
                           LVT_rc%ngrid,&
                           stats%ensspread%model_value_total(:,k,:), &
                           stats%ensspread%count_model_value_total(:,k,:),&
                           stats%standard_name,&
                           stats%ensspread%model_value_ci(k,:))
                   endif
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_writeMetric_ensSpread

!BOP
! 
! !ROUTINE: LVT_resetMetric_ensSpread)
! \label(LVT_resetMetric_ensSpread)
!
! !INTERFACE:
  subroutine LVT_resetMetric_ensSpread(alarm)
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
    type(LVT_statsEntry)   , pointer :: stats

    call LVT_getDataStream1Ptr(model)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))      
       if(stats%selectOpt.eq.1) then 
          do k=1,model%selectNlevs
             if(LVT_metrics%ensspread%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   stats%ensspread%model_value_ts(:,k,l) = 0.0
                   stats%ensspread%count_model_value_ts(:,k,l)=0
                   if(alarm) then 
                      stats%ensspread%tavg_model_value_ts(:,k,l) = 0.0
                      stats%ensspread%tavg_count_model_value_ts(:,k,l)=0
                   endif
                enddo
             endif
          enddo
       endif
    enddo

  end subroutine LVT_resetMetric_ensSpread


!BOP
! 
! !ROUTINE: LVT_writerestart_ensSpread
! 
! !INTERFACE:
  subroutine LVT_writerestart_ensSpread(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for ensSpread metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ensSpread%selectOpt.eq.1) then 
       
       print*, 'The writerestart method is not implemented for ensSpread'

    end if
    
  end subroutine LVT_writerestart_ensSpread

!BOP
! 
! !ROUTINE: LVT_readrestart_ensSpread
! 
! !INTERFACE:
  subroutine LVT_readrestart_ensSpread(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for ensSpread metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ensSpread%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for ensSpread'
       stop
    end if
    
  end subroutine LVT_readrestart_ensSpread

end module LVT_ensSpreadMod
