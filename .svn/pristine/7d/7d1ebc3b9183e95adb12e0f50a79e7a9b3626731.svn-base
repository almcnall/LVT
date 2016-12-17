!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_MAEMod
! \label(LVT_MAEMod)
!
! !INTERFACE:
module LVT_MAEMod
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
!  This module handles the MAE computations by comparing 
!  the LIS output to the specified observations. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP

  private

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initMAE
  public :: LVT_diagnoseMAE
  public :: LVT_computeMAE 
  public :: LVT_writeMetric_MAE
  public :: LVT_resetMetric_MAE
  public :: LVT_writerestart_MAE
  public :: LVT_readrestart_MAE
contains
  
!BOP
! 
! !ROUTINE: LVT_initMAE
! \label{LVT_initMAE}
!
! !INTERFACE: 
  subroutine LVT_initMAE(model,obs,stats,metric)
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
       allocate(stats%mae%value_total(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))
       stats%mae%value_total = 0.0
       allocate(stats%mae%count_value_total(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       
       stats%mae%count_value_total = 0
       allocate(stats%mae%value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%mae%value_ci = LVT_rc%udef
    endif
    
    if(metric%timeopt.eq.1) then 
       allocate(stats%mae%value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%mae%value_ts = 0.0
       allocate(stats%mae%count_value_ts(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       
       stats%mae%count_value_ts = 0 

       allocate(stats%mae%tavg_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%mae%tavg_value_ts = 0.0
       allocate(stats%mae%tavg_count_value_ts(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       
       stats%mae%tavg_count_value_ts = 0 

       if(metric%computeSC.eq.1) then 
          allocate(stats%mae%value_asc(LVT_rc%ngrid, model%selectNlevs,LVT_rc%nasc))
          stats%mae%value_asc = 0.0
          allocate(stats%mae%count_value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%mae%count_value_asc = 0
       endif
       if(metric%computeADC.eq.1) then 
          allocate(stats%mae%value_adc(LVT_rc%ngrid, model%selectNlevs,LVT_rc%nadc))
          stats%mae%value_adc = 0.0
          allocate(stats%mae%count_value_adc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nadc))
          stats%mae%count_value_adc = 0
       endif
    endif
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1
    metric%obsData = .false. 
    metric%stdevFlag = .false. 

  end subroutine LVT_initMAE
  
!BOP
! 
! !ROUTINE: LVT_diagnoseMAE
! \label{LVT_diagnoseMAE}
!
! !INTERFACE: 
  subroutine LVT_diagnoseMAE(pass)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the MAE calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleMAE](\ref{diagnoseSingleMAE})
!     updates the MAE computation for a single variable 
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
       if(LVT_metrics%mae%selectOpt.eq.1.or.&
            LVT_metrics%mae%timeOpt.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             
             call diagnoseSingleMAE(obs,model,stats,&
                  LVT_metrics%mae)

             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
       endif
    endif
  end subroutine LVT_diagnoseMAE

!BOP
! 
! !ROUTINE: diagnoseSingleMAE
! \label{diagnoseSingleMAE}
!
! !INTERFACE: 
  subroutine diagnoseSingleMAE(obs, model, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the MAE computation (updates the running 
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
             if(trim(obs%units).eq.trim(model%units)) then 
                if(obs%count(t,k).ne.0.and. &
                     model%count(t,k).ne.0) then      
                   if(metric%selectOpt.eq.1) then
                      stats%mae%value_total(t,k,1) = stats%mae%value_total(t,k,1) + &
                           abs((-obs%value(t,k) +model%value(t,k)))
                      stats%mae%count_value_total(t,k,1) = &
                           stats%mae%count_value_total(t,k,1) + 1
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then 
                            stats%mae%value_total(t,k,2) = &
                                 stats%mae%value_total(t,k,2) + &
                                 abs((-obs%value(t,k) +model%value(t,k)))
                            stats%mae%count_value_total(t,k,2) = &
                                 stats%mae%count_value_total(t,k,2) + 1
                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then 
                            stats%mae%value_total(t,k,3) = &
                                 stats%mae%value_total(t,k,3) + &
                                 abs((-obs%value(t,k) +model%value(t,k)))
                            stats%mae%count_value_total(t,k,3) = & 
                                 stats%mae%count_value_total(t,k,3) + 1
                         endif
                      endif
                   endif
                   if(metric%timeOpt.eq.1) then 
                      stats%mae%value_ts(t,k,1) = stats%mae%value_ts(t,k,1) + &
                           abs((-obs%value(t,k) +model%value(t,k)))
                      stats%mae%count_value_ts(t,k,1) = stats%mae%count_value_ts(t,k,1)+1
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then
                            stats%mae%value_ts(t,k,2) = stats%mae%value_ts(t,k,2) + &
                                 abs((-obs%value(t,k) +model%value(t,k)))
                            stats%mae%count_value_ts(t,k,2) = stats%mae%count_value_ts(t,k,2)+1
                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then 
                            stats%mae%value_ts(t,k,3) = stats%mae%value_ts(t,k,3) + &
                                 abs((-obs%value(t,k) + model%value(t,k)))
                            stats%mae%count_value_ts(t,k,3) = stats%mae%count_value_ts(t,k,3)+1
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
    
  end subroutine diagnoseSingleMAE


!BOP
! 
! !ROUTINE: LVT_computeMAE
! \label{LVT_computeMAE}
!
! !INTERFACE: 
  subroutine LVT_computeMAE(pass,alarm)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute MAE values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleMAE](\ref{computeSingleMAE})
!     updates the MAE computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     MAE computation has been reached
!   \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: pass

    logical     :: alarm

    integer     :: i
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%mae%selectOpt.eq.1.or.LVT_metrics%mae%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%mae%timeOpt.eq.1.and.&
                  LVT_metrics%mae%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%mae%ftn_ts_loc(i),200,advance='no') &
                        LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                        LVT_rc%hr,'',LVT_rc%mn, '' 
                enddo
             endif
          endif
200    format(I4, a1, I2.2, a1, I2.2, a1, I2.2, a1, I2.2,a1)


          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call computeSingleMAE(alarm,obs,model,stats,&
                  LVT_metrics%mae)

             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
          
          if(alarm) then 
             if(LVT_metrics%mae%timeOpt.eq.1.and.&
                  LVT_metrics%mae%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%mae%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif

  end subroutine LVT_ComputeMAE

!BOP
! 
! !ROUTINE: computeSingleMAE
! \label{computeSingleMAE}
!
! !INTERFACE: 
  subroutine computeSingleMAE(alarm,obs, model,stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the MAE values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     MAE computation has been reached
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

    if(metric%timeOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%mae%count_value_ts(t,k,l).ne.0) then 
                      stats%mae%value_ts(t,k,l) = stats%mae%value_ts(t,k,l)/&
                           stats%mae%count_value_ts(t,k,l)

                      stats%mae%tavg_value_ts(t,k,l) = stats%mae%tavg_value_ts(t,k,l)/&
                           stats%mae%value_ts(t,k,l)
                      stats%mae%tavg_count_value_ts(t,k,l) = & 
                           stats%mae%tavg_count_value_ts(t,k,l) + 1

                   else
                      stats%mae%value_ts(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   if(stats%mae%count_value_ts(t,k,1).ne.0) then 
                      call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                           tind)
                      stats%mae%value_asc(t,k,tind) = &
                           stats%mae%value_asc(t,k,tind)+ &
                           stats%mae%value_ts(t,k,1)
                      stats%mae%count_value_asc(t,k,tind) = &
                           stats%mae%count_value_asc(t,k,tind)+ 1
                   endif
                endif
                if(metric%computeADC.eq.1) then 
                   if(stats%mae%count_value_ts(t,k,1).ne.0) then 
                      call LVT_getADCTimeIndex(tind)
                      stats%mae%value_adc(t,k,tind) = &
                           stats%mae%value_adc(t,k,tind)+ &
                           stats%mae%value_ts(t,k,1)
                      stats%mae%count_value_adc(t,k,tind) = &
                           stats%mae%count_value_adc(t,k,tind)+ 1
                   endif                
                endif                
             enddo
          enddo
          if(alarm) then 
             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%mae%tavg_count_value_ts(t,k,l).ne.0) then 
                         stats%mae%tavg_value_ts(t,k,l) = stats%mae%tavg_value_ts(t,k,l)/&
                              stats%mae%tavg_count_value_ts(t,k,l) 
                      else
                         stats%mae%tavg_value_ts(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                enddo
             enddo

             if(metric%extractTS.eq.1) then 
                call LVT_writeTSinfo(metric%ftn_ts_loc,&
                     model,&
                     LVT_rc%ngrid,&
                     stats%mae%tavg_value_ts,&
                     stats%mae%tavg_count_value_ts)
             endif
          endif
       endif
    endif


    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%mae%count_value_total(t,k,l).ne.0.and.&
                        stats%mae%count_value_total(t,k,l)&
                        .gt.LVT_rc%obsCountThreshold) then 
                      
                      stats%mae%value_total(t,k,l) = stats%mae%value_total(t,k,l)/&
                           stats%mae%count_value_total(t,k,l)
                   else
                      stats%mae%value_total(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1,LVT_rc%nasc
                      if(stats%mae%count_value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then
                         stats%mae%value_asc(t,k,l) = stats%mae%value_asc(t,k,l)/&
                              stats%mae%count_value_asc(t,k,l)
                      else
                         stats%mae%value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
                if(metric%computeADC.eq.1) then 
                   do l=1,LVT_rc%nadc
                      if(stats%mae%count_value_adc(t,k,l).gt.&
                           LVT_rc%ADCCountThreshold) then
                         stats%mae%value_adc(t,k,l) = stats%mae%value_adc(t,k,l)/&
                              stats%mae%count_value_adc(t,k,l)
                      else
                         stats%mae%value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
             enddo
          enddo
          
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%mae%value_total(:,k,l),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%mae%value_ci(k,l))
             enddo
          enddo
       
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid,stats%mae%value_total)      
          if(metric%computeSC.eq.1) then 
             call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%mae%value_asc,stats%mae%count_value_asc)
          endif
          if(metric%computeADC.eq.1) then 
             call LVT_writeAvgDiurnalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,stats%mae%value_adc,stats%mae%count_value_adc)
          endif
       endif
    endif

  end subroutine computeSingleMAE

!BOP
! 
! !ROUTINE: LVT_writeMetric_MAE
! \label{LVT_writeMetric_MAE}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_MAE(pass,final,vlevels,stats,obs)
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

    if(pass.eq.LVT_metrics%mae%npass) then 
       if(final.ne.1) then 
          if(stats%selectOpt.eq.1) then 
             do k=1,vlevels
                if(LVT_metrics%mae%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_gridded(LVT_metrics%mae%ftn_ts, &
                           stats%mae%value_ts(:,k,l),stats%vid_ts(LVT_MAEid,1),1)
                      call LVT_writevar_gridded(LVT_metrics%mae%ftn_ts, &
                           real(stats%mae%count_value_ts(:,k,l)),&
                           stats%vid_count_ts(LVT_MAEid,1),1)
                   enddo
                endif
             enddo
          endif
       else
          if(pass.eq.1) then 
             if(LVT_metrics%mae%selectOpt.eq.1) then
                if(stats%selectOpt.eq.1) then 
                   do k=1,vlevels
                      do l=1,LVT_rc%strat_nlevels
                         
                         call LVT_writevar_gridded(LVT_metrics%mae%ftn_total, &
                              stats%mae%value_total(:,k,l),stats%vid_total(LVT_MAEid,1))
                         call LVT_writevar_gridded(LVT_metrics%mae%ftn_total, &
                              real(stats%mae%count_value_total(:,k,l)),&
                              stats%vid_count_total(LVT_MAEid,1))
                         
                      enddo
                      if(LVT_metrics%mae%computeSC.eq.1) then 
                         do tind = 1,LVT_rc%nasc
                            call LVT_writevar_gridded(LVT_metrics%mae%ftn_total,&
                                 stats%mae%value_asc(:,k,tind),&
                                 stats%vid_sc_total(tind,LVT_MAEid,1))
                         enddo
                      endif
                      
                      if(LVT_metrics%mae%computeADC.eq.1) then 
                         do tind = 1,LVT_rc%nadc
                            call LVT_writevar_gridded(LVT_metrics%mae%ftn_total,&
                                 stats%mae%value_adc(:,k,tind),&
                                 stats%vid_adc_total(tind,LVT_MAEid,1))
                         enddo
                      endif
                      
                      call LVT_writeSummaryStats(&
                           LVT_metrics%mae%ftn_summ,&
                           LVT_metrics%mae%short_name,&
                           LVT_rc%ngrid,&
                           stats%mae%value_total(:,k,:), &
                           stats%mae%count_value_total(:,k,:),stats%standard_name,&
                           stats%mae%value_ci(k,:))
                   enddo
                endif
             endif
          endif
       end if
    endif
  end subroutine LVT_writeMetric_MAE


!BOP
! 
! !ROUTINE: LVT_resetMetric_MAE
! \label(LVT_resetMetric_MAE)
!
! !INTERFACE:
  subroutine LVT_resetMetric_MAE(alarm)
! 
! !INPUT PARAMETERS: 
    logical         :: alarm
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
    integer                :: i,k,l
    type(LVT_metadataEntry), pointer :: model
    type(LVT_statsEntry)   , pointer :: stats


    call LVT_getDataStream1Ptr(model)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(stats%selectOpt.eq.1) then 
          do k=1,model%selectNlevs
             if(LVT_metrics%mae%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   stats%mae%value_ts(:,k,l) = 0.0
                   stats%mae%count_value_ts(:,k,l)=0 
                enddo
             endif
          enddo
       endif

       model => model%next
       stats => stats%next
       
    enddo

  end subroutine LVT_resetMetric_MAE


!BOP
! 
! !ROUTINE: LVT_writerestart_MAE
! 
! !INTERFACE:
  subroutine LVT_writerestart_MAE(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for MAE metric computations
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
    
    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))      
       if(LVT_metrics%mae%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   
                   call LVT_writevar_restart(ftn,&
                        stats%mae%value_total(:,k,l))                   
                   call LVT_writevar_restart(ftn,&
                        stats%mae%count_value_total(:,k,l))
                enddo
             enddo
             if(LVT_metrics%mae%computeSC.eq.1) then
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc  
                      call LVT_writevar_restart(ftn,stats%mae%value_asc(:,k,l))
                      call LVT_writevar_restart(ftn,stats%mae%count_value_asc(:,k,l))
                   end do
                end do
             endif
             if(LVT_metrics%mae%computeADC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nadc  
                      call LVT_writevar_restart(ftn,stats%mae%value_adc(:,k,l))
                      call LVT_writevar_restart(ftn,stats%mae%count_value_adc(:,k,l))
                   enddo
                enddo
             endif
          endif
       endif

       model => model%next
       obs => obs%next
       stats => stats%next
    enddo
  end subroutine LVT_writerestart_MAE


!BOP
! 
! !ROUTINE: LVT_readrestart_MAE
! 
! !INTERFACE:
  subroutine LVT_readrestart_MAE(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for MAE metric computations
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
    
    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))      
       if(LVT_metrics%mae%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   
                   call LVT_readvar_restart(ftn,&
                        stats%mae%value_total(:,k,l))
                   call LVT_readvar_restart(ftn,&
                        stats%mae%count_value_total(:,k,l))
                enddo
             enddo
             if(LVT_metrics%mae%computeSC.eq.1) then
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc  
                      call LVT_readvar_restart(ftn,stats%mae%value_asc(:,k,l))
                      call LVT_readvar_restart(ftn,stats%mae%count_value_asc(:,k,l))
                   end do
                end do
             endif
             if(LVT_metrics%mae%computeADC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nadc  
                      call LVT_readvar_restart(ftn,stats%mae%value_adc(:,k,l))
                      call LVT_readvar_restart(ftn,stats%mae%count_value_adc(:,k,l))
                   enddo
                enddo
             endif
          endif
       endif
       model => model%next
       obs => obs%next
       stats => stats%next
    enddo
  end subroutine LVT_readrestart_MAE
end module LVT_MAEMod
