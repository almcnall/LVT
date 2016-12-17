!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_ensStdevMod
! \label(LVT_ensStdevMod)
!
! !INTERFACE:
module LVT_ensStdevMod
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
!   ensemble standard deviation values
!   of desired variables from the LIS output. 
!   
!   NOTES: 
!   * The LIS output should be written in a tile space format
!   * Observation standard deviation is not computed through this 
!     module, as this routine is computing the ensemble standard
!     deviation. To compute observation standard deviation (which 
!     the standard deviation of observations in time), use the
!     LVT_StdevMod.F90 routines
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
  public :: LVT_initensStdev
  public :: LVT_diagnoseensStdev
  public :: LVT_computeensStdev
  public :: LVT_writeMetric_ensStdev
  public :: LVT_resetMetric_ensStdev
  public :: LVT_writerestart_ensStdev
  public :: LVT_readrestart_ensStdev

!EOP
  
  private

contains
  subroutine LVT_initensStdev(model, obs, stats,metric)
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
       allocate(stats%ensstdev%model_value_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensstdev%model_value_total = 0.0
       allocate(stats%ensstdev%count_model_value_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensstdev%count_model_value_total = 0
       allocate(stats%ensstdev%model_value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%ensstdev%model_value_ci = LVT_rc%udef
       
       if(metric%computeSC.eq.1) then 
          allocate(stats%ensstdev%model_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%ensstdev%model_value_asc = 0.0
          allocate(stats%ensstdev%count_model_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%ensstdev%count_model_value_asc = 0
       endif
       if(metric%computeADC.eq.1) then 
          allocate(stats%ensstdev%model_value_adc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nadc))
          stats%ensstdev%model_value_adc = 0.0
          allocate(stats%ensstdev%count_model_value_adc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nadc))
          stats%ensstdev%count_model_value_adc = 0
       endif
    endif

    if(metric%timeOpt.eq.1) then 
       allocate(stats%ensstdev%model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensstdev%model_value_ts = 0.0
       allocate(stats%ensstdev%count_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensstdev%count_model_value_ts = 0

       allocate(stats%ensstdev%tavg_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensstdev%tavg_model_value_ts = 0.0
       allocate(stats%ensstdev%tavg_count_model_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensstdev%tavg_count_model_value_ts = 0
       
    endif

!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1    
    metric%obsData = .false. 
    metric%stdevFlag = .false. 

  end subroutine LVT_initensStdev

!BOP
! 
! !ROUTINE: LVT_diagnoseensStdev
! \label{LVT_diagnoseensStdev}
!
! !INTERFACE: 
  subroutine LVT_diagnoseensStdev(pass)
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
!    \item[diagnoseSingleModelensStdev](\ref{diagnoseSingleModelensStdev})
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
       if(LVT_metrics%ensstdev%selectOpt.eq.1.or.&
            LVT_metrics%ensstdev%timeOpt.eq.1) then 
          
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleensStdev(model, obs, stats, &
                  LVT_metrics%ensstdev)

             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseensStdev

!BOP
! 
! !ROUTINE: diagnoseSingleensStdev
! \label{diagnoseSingleensStdev}
!
! !INTERFACE: 
  subroutine diagnoseSingleensStdev(model, obs, stats,metric)
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
    real       :: sx, sxx, std_v
    integer    :: nval

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then        
       do g=1,LVT_rc%ngrid
          do k=1, model%selectNlevs
             sx = 0 
             sxx = 0 
             nval = 0 
             do m=1,LVT_LIS_rc(1)%nensem
                t = (g-1)*LVT_LIS_rc(1)%nensem+m
                if(model%count(t,k).gt.0) then 
                   if(metric%selectOpt.eq.1) then 
                      if(model%value(t,k).ne.LVT_rc%udef) then 
                         sx = sx + model%value(t,k)
                         sxx = sxx + model%value(t,k)**2
                         nval = nval + 1
                      endif
                   endif
                endif
             enddo
             if(nval.gt.0) then 
                std_v = sxx/nval - (sx/nval)**2
                if(std_v.ge.0) then 
                   std_v = sqrt(std_v)
                   stats%ensstdev%model_value_total(g,k,1) = & 
                        stats%ensstdev%model_value_total(g,k,1) + std_v
                   stats%ensstdev%count_model_value_total(g,k,1) = & 
                        stats%ensstdev%count_model_value_total(g,k,1) + 1
                   
                   if(metric%timeOpt.eq.1) then 
                      stats%ensstdev%model_value_ts(g,k,1) = & 
                           stats%ensstdev%model_value_ts(g,k,1) + std_v
                      stats%ensstdev%count_model_value_ts(g,k,1) = &
                           stats%ensstdev%count_model_value_ts(g,k,1) + 1
                      
                   endif
                   if(metric%computeSC.eq.1) then 
                      call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                           tind)      
                      stats%ensstdev%model_value_asc(g,k,tind) = & 
                           stats%ensstdev%model_value_asc(g,k,tind) + std_v
                      stats%ensstdev%count_model_value_asc(g,k,tind) = & 
                           stats%ensstdev%count_model_value_asc(g,k,tind) + 1
                   endif
                   if(metric%computeADC.eq.1) then 
                      call LVT_getADCTimeIndex(tind)      
                      stats%ensstdev%model_value_adc(g,k,tind) = & 
                           stats%ensstdev%model_value_adc(g,k,tind) + std_v
                      stats%ensstdev%count_model_value_adc(g,k,tind) = & 
                           stats%ensstdev%count_model_value_adc(g,k,tind) + 1
                   endif
                   if(LVT_rc%strat_nlevels.gt.1) then 
                      if(LVT_stats%strat_var(g,k).gt.&
                           LVT_rc%strat_var_threshold) then
                         stats%ensstdev%model_value_total(g,k,2) = &
                              stats%ensstdev%model_value_total(g,k,2) + std_v
                         stats%ensstdev%count_model_value_total(g,k,2) = &
                              stats%ensstdev%count_model_value_total(g,k,2) + 1
                         
                         if(metric%timeOpt.eq.1) then 
                            stats%ensstdev%model_value_ts(g,k,2) = & 
                                 stats%ensstdev%model_value_ts(g,k,2) + std_v
                            stats%ensstdev%count_model_value_ts(g,k,2) = &
                                 stats%ensstdev%count_model_value_ts(g,k,2) + 1
                            
                         endif
                      elseif(LVT_stats%strat_var(g,k).le.&
                           LVT_rc%strat_var_threshold) then
                         stats%ensstdev%model_value_total(g,k,3) = &
                              stats%ensstdev%model_value_total(g,k,3) + std_v
                         stats%ensstdev%count_model_value_total(g,k,3) = &
                              stats%ensstdev%count_model_value_total(g,k,3) + 1
                         
                         if(metric%timeOpt.eq.1) then 
                            stats%ensstdev%model_value_ts(g,k,3) = & 
                                 stats%ensstdev%model_value_ts(g,k,3) + std_v
                            stats%ensstdev%count_model_value_ts(g,k,3) = &
                                 stats%ensstdev%count_model_value_ts(g,k,3) + 1
                            
                         endif
                      endif
                   endif
                endif
             endif
          enddo
       enddo
    endif
  end subroutine diagnoseSingleensStdev

!BOP
! 
! !ROUTINE: LVT_computeensStdev
! \label{LVT_computeensStdev}
!
! !INTERFACE: 
  subroutine LVT_computeensStdev(pass,alarm)
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
!    \item[computeSingleModelensStdev](\ref{computeSingleModelensStdev})
!     computes the std values for a single variable
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer     :: i 
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats
    
    if(pass.eq.1) then 
       if(LVT_metrics%ensstdev%selectOpt.eq.1.or.&
            LVT_metrics%ensstdev%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%ensstdev%timeOpt.eq.1.and.&
                  LVT_metrics%ensstdev%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ensstdev%ftn_ts_loc(i),200,advance='no') &
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
             call computeSingleensStdev(alarm,&
                  model, obs, stats, &
                  LVT_metrics%ensstdev)

             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%ensstdev%timeOpt.eq.1.and.&
                  LVT_metrics%ensstdev%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ensstdev%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_computeensStdev
  

!BOP
! 
! !ROUTINE: computeSingleensStdev
! \label{computeSingleensStdev}
!
! !INTERFACE: 
  subroutine computeSingleensStdev(alarm,model,obs,stats,metric)
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

    integer  :: t,l,k,dummy1,dummy2

    if(metric%timeOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%ensstdev%count_model_value_ts(t,k,l).gt.0) then 
                      stats%ensstdev%model_value_ts(t,k,l) = stats%ensstdev%model_value_ts(t,k,l)/&
                           stats%ensstdev%count_model_value_ts(t,k,l)                   

                      stats%ensstdev%tavg_model_value_ts(t,k,l) = &
                           stats%ensstdev%tavg_model_value_ts(t,k,l) + & 
                           stats%ensstdev%model_value_ts(t,k,l) 
                      stats%ensstdev%tavg_count_model_value_ts(t,k,l) = &
                           stats%ensstdev%tavg_count_model_value_ts(t,k,l) + 1
                   else
                      stats%ensstdev%model_value_ts(t,k,l) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
          
          if(alarm) then 
             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%ensstdev%tavg_count_model_value_ts(t,k,l).gt.0) then 
                         stats%ensstdev%tavg_model_value_ts(t,k,l) = &
                              stats%ensstdev%tavg_model_value_ts(t,k,l) / & 
                              stats%ensstdev%tavg_count_model_value_ts(t,k,l) 
                      else
                         stats%ensstdev%tavg_model_value_ts(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                enddo
             enddo

             if(metric%extractTS.eq.1) then 
                call LVT_writeTSinfo(metric%ftn_ts_loc,&
                     model,&
                     LVT_rc%ngrid,&
                     stats%ensstdev%tavg_model_value_ts,&
                     stats%ensstdev%tavg_count_model_value_ts,dummy1,dummy2)
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
                   if(stats%ensstdev%count_model_value_total(t,k,l).gt.&
                        LVT_rc%obsCountThreshold) then 
                      stats%ensstdev%model_value_total(t,k,l) = &
                           stats%ensstdev%model_value_total(t,k,l)/&
                              stats%ensstdev%count_model_value_total(t,k,l)           
                   else
                      stats%ensstdev%model_value_total(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1, LVT_rc%nasc
                      if(stats%ensstdev%count_model_value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then 
                         stats%ensstdev%model_value_asc(t,k,l) = &
                              stats%ensstdev%model_value_asc(t,k,l)/&
                              stats%ensstdev%count_model_value_asc(t,k,l)           
                      else
                         stats%ensstdev%model_value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif

                if(metric%computeADC.eq.1) then 
                   do l=1, LVT_rc%nadc
                      if(stats%ensstdev%count_model_value_adc(t,k,l).gt.&
                           LVT_rc%ADCCountThreshold) then 
                         stats%ensstdev%model_value_adc(t,k,l) = &
                              stats%ensstdev%model_value_adc(t,k,l)/&
                              stats%ensstdev%count_model_value_adc(t,k,l)           
                      else
                         stats%ensstdev%model_value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
             enddo
          enddo
             
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%ensstdev%model_value_total(:,k,l),&
                     LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%ensstdev%model_value_ci(k,l))
             enddo
          enddo
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid,stats%ensstdev%model_value_total)
          
          if(metric%computeSC.eq.1) then 
             call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%ensstdev%model_value_asc,stats%ensstdev%count_model_value_asc)
          endif
          if(metric%computeADC.eq.1) then 
             call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%ensstdev%model_value_adc,stats%ensstdev%count_model_value_adc)
          endif
       endif
    endif
  end subroutine computeSingleensStdev


  subroutine LVT_writeMetric_ensStdev(pass,final,vlevels,stats,obs)
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

    if(pass.eq.LVT_metrics%ensstdev%npass) then 
       if(final.ne.1) then
           if(stats%selectOpt.eq.1) then 
             do k=1,vlevels
                if(LVT_metrics%ensstdev%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_gridded(LVT_metrics%ensstdev%ftn_ts, &
                           stats%ensstdev%model_value_ts(:,k,l),&
                           stats%vid_ts(LVT_ensStdevid,1),k)
                      call LVT_writevar_gridded(LVT_metrics%ensstdev%ftn_ts, &
                           real(stats%ensstdev%count_model_value_ts(:,k,l)),&
                           stats%vid_count_ts(LVT_ensStdevid,1),k)
                   enddo
                endif
             enddo
          endif
       else
          if(pass.eq.1) then 
             if(stats%selectOpt.eq.1) then 
                do k=1,vlevels
                   if(LVT_metrics%ensstdev%selectOpt.eq.1) then 
                      do l=1,LVT_rc%strat_nlevels                      
                         call LVT_writevar_gridded(&
                              LVT_metrics%ensstdev%ftn_total, &
                              stats%ensstdev%model_value_total(:,k,l),&
                              stats%vid_total(LVT_ensStdevid,1),&
                              k)
                         call LVT_writevar_gridded(&
                              LVT_metrics%ensstdev%ftn_total, &
                              real(stats%ensstdev%count_model_value_total(:,k,l)),&
                              stats%vid_count_total(LVT_ensStdevid,1),&
                              k)
                         
                      enddo
                   
                      if(LVT_metrics%ensstdev%computeSC.eq.1) then 
                         do tind = 1,LVT_rc%nasc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%ensstdev%ftn_total,&
                                 stats%ensstdev%model_value_asc(:,k,tind),&
                                 stats%vid_sc_total(tind,LVT_ensStdevid,1),&
                                 1)
                         enddo
                      endif
                      if(LVT_metrics%ensstdev%computeADC.eq.1) then 
                         do tind = 1,LVT_rc%nadc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%ensstdev%ftn_total,&
                                 stats%ensstdev%model_value_adc(:,k,tind),&
                                 stats%vid_adc_total(tind,LVT_ensStdevid,1),&
                                 1)
                         enddo
                      endif
                      call LVT_writeSummaryStats(&
                           LVT_metrics%ensstdev%ftn_summ,&
                           LVT_metrics%ensstdev%short_name,&
                           LVT_rc%ngrid,&
                           stats%ensstdev%model_value_total(:,k,:), &
                           stats%ensstdev%count_model_value_total(:,k,:),&
                           stats%standard_name,&
                           stats%ensstdev%model_value_ci(k,:))
                   endif
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_writeMetric_ensStdev

!BOP
! 
! !ROUTINE: LVT_resetMetric_ensStdev)
! \label(LVT_resetMetric_ensStdev)
!
! !INTERFACE:
  subroutine LVT_resetMetric_ensStdev(alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
    logical           :: alarm 
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
             if(LVT_metrics%ensstdev%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   stats%ensstdev%model_value_ts(:,k,l) = 0.0
                   stats%ensstdev%count_model_value_ts(:,k,l)=0 
                   if(alarm) then 
                      stats%ensstdev%tavg_model_value_ts(:,k,l) = 0.0
                      stats%ensstdev%tavg_count_model_value_ts(:,k,l)=0 
                   endif
                enddo
             endif
          enddo
       endif

       model => model%next
       stats => stats%next
       
    enddo

  end subroutine LVT_resetMetric_ensStdev


!BOP
! 
! !ROUTINE: LVT_writerestart_ensStdev
! 
! !INTERFACE:
  subroutine LVT_writerestart_ensStdev(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for ensStdev metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ensstdev%selectOpt.eq.1) then 
       
       print*, 'The writerestart method is not implemented for ensStdev'

    end if
    
  end subroutine LVT_writerestart_ensStdev

!BOP
! 
! !ROUTINE: LVT_readrestart_ensStdev
! 
! !INTERFACE:
  subroutine LVT_readrestart_ensStdev(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for ensStdev metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ensstdev%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for ensStdev'
       stop
    end if
    
  end subroutine LVT_readrestart_ensStdev

end module LVT_ensStdevMod
