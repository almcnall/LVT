!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_ensMeanBiasMod
! \label(LVT_ensMeanBiasMod)
!
! !INTERFACE:
module LVT_ensMeanBiasMod
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
!   ensemble Mean Error (K) defined as follows:
!     
!    K = [ sqrt(1/N sum ( Xi - o)^2 ) ]
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
  public :: LVT_initensMeanBias
  public :: LVT_diagnoseensMeanBias
  public :: LVT_computeensMeanBias
  public :: LVT_writeMetric_ensMeanBias
  public :: LVT_resetMetric_ensMeanBias
  public :: LVT_writerestart_ensMeanBias
  public :: LVT_readrestart_ensMeanBias
  
!EOP
  
  private

contains
  subroutine LVT_initensMeanBias(model, obs, stats,metric)
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
       allocate(stats%ensmbias%value_final(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensmbias%value_final = 0.0
       allocate(stats%ensmbias%count_value_final(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensmbias%count_value_final = 0
       allocate(stats%ensmbias%value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%ensmbias%value_ci = LVT_rc%udef
       
       if(metric%computeSC.eq.1) then 
          allocate(stats%ensmbias%value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%ensmbias%value_asc = 0.0
          allocate(stats%ensmbias%count_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%ensmbias%count_value_asc = 0
       endif
       if(metric%computeADC.eq.1) then 
          allocate(stats%ensmbias%value_adc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nadc))
          stats%ensmbias%value_adc = 0.0
          allocate(stats%ensmbias%count_value_adc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nadc))
          stats%ensmbias%count_value_adc = 0
       endif
    endif

    if(metric%timeOpt.eq.1) then 
       allocate(stats%ensmbias%value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensmbias%value_ts = 0.0
       allocate(stats%ensmbias%count_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensmbias%count_value_ts = 0

       allocate(stats%ensmbias%tavg_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensmbias%tavg_value_ts = 0.0
       allocate(stats%ensmbias%tavg_count_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensmbias%tavg_count_value_ts = 0
       
    endif

!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1    
    metric%obsData = .false. 
    metric%stdevFlag = .false. 

  end subroutine LVT_initensMeanBias

!BOP
! 
! !ROUTINE: LVT_diagnoseensMeanBias
! \label{LVT_diagnoseensMeanBias}
!
! !INTERFACE: 
  subroutine LVT_diagnoseensMeanBias(pass)
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
!    \item[diagnoseSingleModelensMeanBias](\ref{diagnoseSingleModelensMeanBias})
!     updates the std computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    integer       :: index
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%ensmeanbias%selectOpt.eq.1.or.&
            LVT_metrics%ensmeanbias%timeOpt.eq.1) then 
          
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call diagnoseSingleensMeanBias(&
                  model, obs, stats, &
                  LVT_metrics%ensmeanbias)
             model => model%next
             obs   => obs%next
             stats => stats%next

          enddo
       endif
    endif
  end subroutine LVT_diagnoseensMeanBias

!BOP
! 
! !ROUTINE: diagnoseSingleensMeanBias
! \label{diagnoseSingleensMeanBias}
!
! !INTERFACE: 
  subroutine diagnoseSingleensMeanBias(model, obs, stats,metric)
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
    real       :: sumval,me,ens_mean
    integer    :: nsumval

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1.and.&
         obs%selectNlevs.ge.1) then        
       if(trim(obs%units).eq.trim(model%units)) then              
          do g=1,LVT_rc%ngrid
             do k=1, model%selectNlevs
                sumval = 0 
                nsumval = 0 
                do m=1,LVT_LIS_rc(1)%nensem
                   t = (g-1)*LVT_LIS_rc(1)%nensem+m   
                   if(obs%count(g,k).ne.0.and. &
                        model%count(g,k).ne.0) then      
                      if(metric%selectOpt.eq.1) then 
                         if(model%value(g,k).ne.LVT_rc%udef) then 
                            sumval = sumval + model%value(t,k)
                            nsumval = nsumval + 1
                         endif
                      endif
                   endif
                enddo
                if(nsumval.gt.0) then 
                   ens_mean = sumval/nsumval
                   me = (ens_mean - obs%value(g,k))
                   stats%ensmbias%value_final(g,k,1) = & 
                        stats%ensmbias%value_final(g,k,1) +&
                        me
                   stats%ensmbias%count_value_final(g,k,1) = & 
                        stats%ensmbias%count_value_final(g,k,1) + 1
                   
                   if(metric%timeOpt.eq.1) then 
                      stats%ensmbias%value_ts(g,k,1) = & 
                           stats%ensmbias%value_ts(g,k,1) + me
                      stats%ensmbias%count_value_ts(g,k,1) = &
                           stats%ensmbias%count_value_ts(g,k,1) + 1
                      
                   endif
                   if(metric%computeSC.eq.1) then 
                      call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                           tind)      
                      stats%ensmbias%value_asc(g,k,tind) = & 
                           stats%ensmbias%value_asc(g,k,tind) + me
                      stats%ensmbias%count_value_asc(g,k,tind) = & 
                           stats%ensmbias%count_value_asc(g,k,tind) + 1
                   endif
                   if(metric%computeADC.eq.1) then 
                      call LVT_getADCTimeIndex(tind)      
                      stats%ensmbias%value_adc(g,k,tind) = & 
                           stats%ensmbias%value_adc(g,k,tind) + me
                      stats%ensmbias%count_value_adc(g,k,tind) = & 
                           stats%ensmbias%count_value_adc(g,k,tind) + 1
                   endif
                   if(LVT_rc%strat_nlevels.gt.1) then 
                      if(LVT_stats%strat_var(g,k).gt.&
                           LVT_rc%strat_var_threshold) then
                         stats%ensmbias%value_final(g,k,2) = &
                              stats%ensmbias%value_final(g,k,2) + me
                         stats%ensmbias%count_value_final(g,k,2) = &
                              stats%ensmbias%count_value_final(g,k,2) + 1
                         
                         if(metric%timeOpt.eq.1) then 
                            stats%ensmbias%value_ts(g,k,2) = & 
                                 stats%ensmbias%value_ts(g,k,2) + me
                            stats%ensmbias%count_value_ts(g,k,2) = &
                                 stats%ensmbias%count_value_ts(g,k,2) + 1
                            
                         endif
                      elseif(LVT_stats%strat_var(g,k).le.&
                           LVT_rc%strat_var_threshold) then
                         stats%ensmbias%value_final(g,k,3) = &
                              stats%ensmbias%value_final(g,k,3) + me
                         stats%ensmbias%count_value_final(g,k,3) = &
                              stats%ensmbias%count_value_final(g,k,3) + 1
                         
                         if(metric%timeOpt.eq.1) then 
                            stats%ensmbias%value_ts(g,k,3) = & 
                                 stats%ensmbias%value_ts(g,k,3) + me
                            stats%ensmbias%count_value_ts(g,k,3) = &
                                 stats%ensmbias%count_value_ts(g,k,3) + 1
                            
                         endif
                      endif
                   endif
                end if
             enddo
          enddo
       endif
    endif
  end subroutine diagnoseSingleensMeanBias

!BOP
! 
! !ROUTINE: LVT_computeensMeanBias
! \label{LVT_computeensMeanBias}
!
! !INTERFACE: 
  subroutine LVT_computeensMeanBias(pass,alarm)
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
!    \item[computeSingleModelensMeanBias](\ref{computeSingleModelensMeanBias})
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
       if(LVT_metrics%ensmeanbias%selectOpt.eq.1.or.&
            LVT_metrics%ensmeanbias%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%ensmeanbias%timeOpt.eq.1.and.&
                  LVT_metrics%ensmeanbias%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ensmeanbias%ftn_ts_loc(i),200,advance='no') &
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

             call computeSingleensMeanBias(alarm,&
                  model, obs, stats, &
                  LVT_metrics%ensmeanbias)
          enddo
          
          if(alarm) then 
             if(LVT_metrics%ensmeanbias%timeOpt.eq.1.and.&
                  LVT_metrics%ensmeanbias%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ensmeanbias%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_computeensMeanBias
  

!BOP
! 
! !ROUTINE: computeSingleensMeanBias
! \label{computeSingleensMeanBias}
!
! !INTERFACE: 
  subroutine computeSingleensMeanBias(alarm,model,obs,stats,metric)
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

    integer  :: t,l,k,tind
    integer  :: dummy1,dummy2

    if(metric%timeOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%ensmbias%count_value_ts(t,k,l).gt.0) then 
                      stats%ensmbias%value_ts(t,k,l) = &
                           stats%ensmbias%value_ts(t,k,l)/&
                           stats%ensmbias%count_value_ts(t,k,l)

                      stats%ensmbias%tavg_value_ts(t,k,l) = &
                           stats%ensmbias%tavg_value_ts(t,k,l) + & 
                           stats%ensmbias%count_value_ts(t,k,l)
                      stats%ensmbias%tavg_count_value_ts(t,k,l) = &
                           stats%ensmbias%tavg_count_value_ts(t,k,l) + 1
                   else
                      stats%ensmbias%value_ts(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1, LVT_rc%nasc
                      if(stats%ensmbias%count_value_ts(t,k,l).ne.0) then
                         call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                              tind) 
                         stats%ensmbias%value_asc(t,k,tind) = &
                              stats%ensmbias%value_asc(t,k,tind)+ &
                              stats%ensmbias%value_ts(t,k,tind)
                         stats%ensmbias%count_value_asc(t,k,tind) = & 
                              stats%ensmbias%count_value_asc(t,k,tind) + 1
                      endif
                   enddo
                endif

                if(metric%computeADC.eq.1) then 
                   do l=1, LVT_rc%nadc
                      if(stats%ensmbias%count_value_ts(t,k,l).ne.0) then 
                         call LVT_getADCTimeIndex(tind)
                         stats%ensmbias%value_adc(t,k,tind) = &
                              stats%ensmbias%value_adc(t,k,tind) + &
                              stats%ensmbias%value_ts(t,k,tind) 
                         stats%ensmbias%count_value_adc(t,k,tind) = & 
                              stats%ensmbias%count_value_adc(t,k,tind) + 1
                      endif
                   enddo
                endif
             enddo
          enddo
          
          if(alarm) then 
             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if( stats%ensmbias%tavg_count_value_ts(t,k,l).gt.0) then 
                         stats%ensmbias%tavg_value_ts(t,k,l) = &
                              stats%ensmbias%tavg_value_ts(t,k,l)/ & 
                              stats%ensmbias%tavg_count_value_ts(t,k,l)
                      endif
                   enddo
                enddo
             enddo

             if(metric%extractTS.eq.1) then 
                call LVT_writeTSinfo(metric%ftn_ts_loc,&
                     model,&
                     LVT_rc%ngrid,&
                     stats%ensmbias%tavg_value_ts,&
                     stats%ensmbias%tavg_count_value_ts,dummy1,dummy2)
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
                   if(stats%ensmbias%count_value_final(t,k,l).gt.&
                        LVT_rc%obsCountThreshold) then 
                      stats%ensmbias%value_final(t,k,l) = &
                           stats%ensmbias%value_final(t,k,l)/&
                           stats%ensmbias%count_value_final(t,k,l)           
                   else
                      stats%ensmbias%value_final(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1, LVT_rc%nasc
                      if(stats%ensmbias%count_value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then 
                         stats%ensmbias%value_asc(t,k,l) = &
                              stats%ensmbias%value_asc(t,k,l)/ &
                              stats%ensmbias%count_value_asc(t,k,l)
                      else
                         stats%ensmbias%value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
                   
                if(metric%computeADC.eq.1) then 
                   do l=1, LVT_rc%nadc
                      if(stats%ensmbias%count_value_adc(t,k,l).gt.&
                           LVT_rc%ADCCountThreshold) then 
                         stats%ensmbias%value_adc(t,k,l) = &
                              stats%ensmbias%value_adc(t,k,l)/ &
                              stats%ensmbias%count_value_adc(t,k,l) 
                      else
                         stats%ensmbias%value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif                   
             enddo
          enddo
             
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%ensmbias%value_final(:,k,l),&
                     LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%ensmbias%value_ci(k,l))
             enddo
          enddo
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid,stats%ensmbias%value_final)
          
          if(metric%computeSC.eq.1) then 
             call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%ensmbias%value_asc,stats%ensmbias%count_value_asc)
          endif
          if(metric%computeADC.eq.1) then 
             call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%ensmbias%value_adc,stats%ensmbias%count_value_adc)
          endif
       endif
    endif
  end subroutine computeSingleensMeanBias


  subroutine LVT_writeMetric_ensMeanBias(pass,final,vlevels,stats,obs)
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

    if(pass.eq.LVT_metrics%ensmeanbias%npass) then 
       if(final.ne.1) then
           if(stats%selectOpt.eq.1) then 
             do k=1,vlevels
                if(LVT_metrics%ensmeanbias%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_gridded(LVT_metrics%ensmeanbias%ftn_ts, &
                           stats%ensmbias%value_ts(:,k,l),&
                           stats%vid_ts(LVT_ensMeanBiasid,1),k)
                      call LVT_writevar_gridded(LVT_metrics%ensmeanbias%ftn_ts, &
                           real(stats%ensmbias%count_value_ts(:,k,l)),&
                           stats%vid_count_ts(LVT_ensMeanBiasid,1),k)
                   enddo
                endif
             enddo
          endif
       else
          if(pass.eq.1) then 
             if(stats%selectOpt.eq.1) then 
                do k=1,vlevels
                   if(LVT_metrics%ensmeanbias%selectOpt.eq.1) then 
                      do l=1,LVT_rc%strat_nlevels                      
                         call LVT_writevar_gridded(&
                              LVT_metrics%ensmeanbias%ftn_total, &
                              stats%ensmbias%value_final(:,k,l),&
                              stats%vid_total(LVT_ensMeanBiasid,1),&
                              k)
                         call LVT_writevar_gridded(&
                              LVT_metrics%ensmeanbias%ftn_total, &
                              real(stats%ensmbias%count_value_final(:,k,l)),&
                              stats%vid_count_total(LVT_ensMeanBiasid,1),&
                              k)
                         
                      enddo
                   
                      if(LVT_metrics%ensmeanbias%computeSC.eq.1) then 
                         do tind = 1,LVT_rc%nasc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%ensmeanbias%ftn_total,&
                                 stats%ensmbias%value_asc(:,k,tind),&
                                 stats%vid_sc_total(tind,LVT_ensMeanBiasid,1),&
                                 1)
                         enddo
                      endif
                      if(LVT_metrics%ensmeanbias%computeADC.eq.1) then 
                         do tind = 1,LVT_rc%nadc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%ensmeanbias%ftn_total,&
                                 stats%ensmbias%value_adc(:,k,tind),&
                                 stats%vid_adc_total(tind,LVT_ensMeanBiasid,1),&
                                 1)
                         enddo
                      endif
                      call LVT_writeSummaryStats(&
                           LVT_metrics%ensmeanbias%ftn_summ,&
                           LVT_metrics%ensmeanbias%short_name,&
                           LVT_rc%ngrid,&
                           stats%ensmbias%value_final(:,k,:), &
                           stats%ensmbias%count_value_final(:,k,:),&
                           stats%standard_name,&
                           stats%ensmbias%value_ci(k,:))
                   endif
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_writeMetric_ensMeanBias

!BOP
! 
! !ROUTINE: LVT_resetMetric_ensMeanBias)
! \label(LVT_resetMetric_ensMeanBias)
!
! !INTERFACE:
  subroutine LVT_resetMetric_ensMeanBias(alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
    logical              :: alarm 
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
             if(LVT_metrics%ensmeanbias%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   stats%ensmbias%value_ts(:,k,l) = 0.0
                   stats%ensmbias%count_value_ts(:,k,l)=0 
                   if(alarm) then 
                      stats%ensmbias%tavg_value_ts(:,k,l) = 0.0
                      stats%ensmbias%tavg_count_value_ts(:,k,l)=0 
                   endif
                enddo
             endif
          enddo
       endif
       
       model => model%next
       stats => stats%next
    enddo

  end subroutine LVT_resetMetric_ensMeanBias


!BOP
! 
! !ROUTINE: LVT_writerestart_ensMeanBias
! 
! !INTERFACE:
  subroutine LVT_writerestart_ensMeanBias(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for ensMeanBias metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ensMeanBias%selectOpt.eq.1) then 
       
       print*, 'The writerestart method is not implemented for ensMeanBias'

    end if
    
  end subroutine LVT_writerestart_ensMeanBias

!BOP
! 
! !ROUTINE: LVT_readrestart_ensMeanBias
! 
! !INTERFACE:
  subroutine LVT_readrestart_ensMeanBias(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for ensMeanBias metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ensMeanBias%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for ensMeanBias'
       stop
    end if
    
  end subroutine LVT_readrestart_ensMeanBias

end module LVT_ensMeanBiasMod
