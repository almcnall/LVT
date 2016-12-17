!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_ensSkillMod
! \label(LVT_ensSkillMod)
!
! !INTERFACE:
module LVT_ensSkillMod
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
!   ensemble Skill (K) defined as follows:
!     
!    K = 1/[ sqrt(1/N sum ( Xi - o)^2 ) ]
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
  public :: LVT_initensSkill
  public :: LVT_diagnoseensSkill
  public :: LVT_computeensSkill
  public :: LVT_writeMetric_ensSkill
  public :: LVT_resetMetric_ensSkill
  public :: LVT_writerestart_ensSkill
  public :: LVT_readrestart_ensSkill
  
!EOP
  
  private

contains
  subroutine LVT_initensSkill(model, obs, stats,metric)
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
       allocate(stats%ensskill%value_final(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensskill%value_final = 0.0
       allocate(stats%ensskill%value_final(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensskill%value_final = 0
       allocate(stats%ensskill%value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%ensskill%value_ci = LVT_rc%udef
       
       if(metric%computeSC.eq.1) then 
          allocate(stats%ensskill%value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%ensskill%value_asc = 0.0
          allocate(stats%ensskill%value_asc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nasc))
          stats%ensskill%value_asc = 0
       endif
       if(metric%computeADC.eq.1) then 
          allocate(stats%ensskill%value_adc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nadc))
          stats%ensskill%value_adc = 0.0
          allocate(stats%ensskill%value_adc(LVT_rc%ngrid, model%selectNlevs,&
               LVT_rc%nadc))
          stats%ensskill%value_adc = 0
       endif
    endif

    if(metric%timeOpt.eq.1) then 
       allocate(stats%ensskill%value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensskill%value_ts = 0.0
       allocate(stats%ensskill%value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensskill%value_ts = 0

       allocate(stats%ensskill%tavg_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensskill%tavg_value_ts = 0.0
       allocate(stats%ensskill%tavg_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%ensskill%tavg_value_ts = 0
       
    endif

!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1    
    metric%obsData = .false. 
    metric%stdevFlag = .false. 

  end subroutine LVT_initensSkill

!BOP
! 
! !ROUTINE: LVT_diagnoseensSkill
! \label{LVT_diagnoseensSkill}
!
! !INTERFACE: 
  subroutine LVT_diagnoseensSkill(pass)
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
!    \item[diagnoseSingleModelensSkill](\ref{diagnoseSingleModelensSkill})
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
       if(LVT_metrics%ensskill%selectOpt.eq.1.or.&
            LVT_metrics%ensskill%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call diagnoseSingleensSkill(model, obs, stats, &
                  LVT_metrics%ensskill)

             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseensSkill

!BOP
! 
! !ROUTINE: diagnoseSingleensSkill
! \label{diagnoseSingleensSkill}
!
! !INTERFACE: 
  subroutine diagnoseSingleensSkill(model, obs, stats,metric)
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
    real       :: sumval,skill,ens_mean
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
                   skill = (ens_mean - obs%value(g,k))**2
                   stats%ensskill%value_final(g,k,1) = & 
                        stats%ensskill%value_final(g,k,1) +&
                        skill
                   stats%ensskill%value_final(g,k,1) = & 
                        stats%ensskill%value_final(g,k,1) + 1
                   
                   if(metric%timeOpt.eq.1) then 
                      stats%ensskill%value_ts(g,k,1) = & 
                           stats%ensskill%value_ts(g,k,1) + skill
                      stats%ensskill%value_ts(g,k,1) = &
                           stats%ensskill%value_ts(g,k,1) + 1
                      
                   endif
                   if(metric%computeSC.eq.1) then 
                      call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                           tind)      
                      stats%ensskill%value_asc(g,k,tind) = & 
                           stats%ensskill%value_asc(g,k,tind) + skill
                      stats%ensskill%value_asc(g,k,tind) = & 
                           stats%ensskill%value_asc(g,k,tind) + 1
                   endif
                   if(metric%computeADC.eq.1) then 
                      call LVT_getADCTimeIndex(tind)      
                      stats%ensskill%value_adc(g,k,tind) = & 
                           stats%ensskill%value_adc(g,k,tind) + skill
                      stats%ensskill%value_adc(g,k,tind) = & 
                           stats%ensskill%value_adc(g,k,tind) + 1
                   endif
                   if(LVT_rc%strat_nlevels.gt.1) then 
                      if(LVT_stats%strat_var(g,k).gt.&
                           LVT_rc%strat_var_threshold) then
                         stats%ensskill%value_final(g,k,2) = &
                              stats%ensskill%value_final(g,k,2) + skill
                         stats%ensskill%value_final(g,k,2) = &
                              stats%ensskill%value_final(g,k,2) + 1
                         
                         if(metric%timeOpt.eq.1) then 
                            stats%ensskill%value_ts(g,k,2) = & 
                                 stats%ensskill%value_ts(g,k,2) + skill
                            stats%ensskill%value_ts(g,k,2) = &
                                 stats%ensskill%value_ts(g,k,2) + 1
                            
                         endif
                      elseif(LVT_stats%strat_var(g,k).le.&
                           LVT_rc%strat_var_threshold) then
                         stats%ensskill%value_final(g,k,3) = &
                              stats%ensskill%value_final(g,k,3) + skill
                         stats%ensskill%value_final(g,k,3) = &
                              stats%ensskill%value_final(g,k,3) + 1
                         
                         if(metric%timeOpt.eq.1) then 
                            stats%ensskill%value_ts(g,k,3) = & 
                                 stats%ensskill%value_ts(g,k,3) + skill
                            stats%ensskill%value_ts(g,k,3) = &
                                 stats%ensskill%value_ts(g,k,3) + 1
                            
                         endif
                      endif
                   endif
                end if
             enddo
          enddo
       endif
    endif
  end subroutine diagnoseSingleensSkill

!BOP
! 
! !ROUTINE: LVT_computeensSkill
! \label{LVT_computeensSkill}
!
! !INTERFACE: 
  subroutine LVT_computeensSkill(pass,alarm)
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
!    \item[computeSingleModelensSkill](\ref{computeSingleModelensSkill})
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
       if(LVT_metrics%ensskill%selectOpt.eq.1.or.&
            LVT_metrics%ensskill%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%ensskill%timeOpt.eq.1.and.&
                  LVT_metrics%ensskill%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ensskill%ftn_ts_loc(i),200,advance='no') &
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

             call computeSingleensSkill(alarm,&
                  model, obs, stats, &
                  LVT_metrics%ensskill)
          enddo
          
          if(alarm) then 
             if(LVT_metrics%ensskill%timeOpt.eq.1.and.&
                  LVT_metrics%ensskill%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ensskill%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_computeensSkill
  

!BOP
! 
! !ROUTINE: computeSingleensSkill
! \label{computeSingleensSkill}
!
! !INTERFACE: 
  subroutine computeSingleensSkill(alarm,model,obs,stats,metric)
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
    real     :: merr
    integer  :: dummy1,dummy2

    if(metric%timeOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%ensskill%value_ts(t,k,l).gt.0.and.&
                        stats%ensskill%value_ts(t,k,l).ne.0) then 
                      stats%ensskill%value_ts(t,k,l) = &
                           1/(sqrt(stats%ensskill%value_ts(t,k,l)/&
                           stats%ensskill%value_ts(t,k,l) ))                  

                      stats%ensskill%tavg_value_ts(t,k,l) = & 
                           stats%ensskill%tavg_value_ts(t,k,l) + & 
                           stats%ensskill%value_ts(t,k,l)
                      stats%ensskill%tavg_count_value_ts(t,k,l) = & 
                           stats%ensskill%tavg_count_value_ts(t,k,l) + 1
                   else
                      stats%ensskill%value_ts(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1, LVT_rc%nasc
                      if(stats%ensskill%value_ts(t,k,l).ne.0) then
                         call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                              tind) 
                         stats%ensskill%value_asc(t,k,tind) = &
                              stats%ensskill%value_asc(t,k,tind)+ &
                              stats%ensskill%value_ts(t,k,tind)
                         stats%ensskill%value_asc(t,k,tind) = & 
                              stats%ensskill%value_asc(t,k,tind) + 1
                      endif
                   enddo
                endif

                if(metric%computeADC.eq.1) then 
                   do l=1, LVT_rc%nadc
                      if(stats%ensskill%value_ts(t,k,l).ne.0) then 
                         call LVT_getADCTimeIndex(tind)
                         stats%ensskill%value_adc(t,k,tind) = &
                              stats%ensskill%value_adc(t,k,tind) + &
                              stats%ensskill%value_ts(t,k,tind) 
                         stats%ensskill%value_adc(t,k,tind) = & 
                              stats%ensskill%value_adc(t,k,tind) + 1
                      endif
                   enddo
                endif
             enddo
          enddo
          
          if(alarm) then 

             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%ensskill%tavg_count_value_ts(t,k,l).gt.0) then 
                         stats%ensskill%tavg_value_ts(t,k,l) = & 
                              stats%ensskill%tavg_value_ts(t,k,l) / & 
                              stats%ensskill%tavg_count_value_ts(t,k,l)
                      endif
                   enddo
                enddo
             enddo

             if(metric%extractTS.eq.1) then 
                call LVT_writeTSinfo(metric%ftn_ts_loc,&
                     model,&
                     LVT_rc%ngrid,&
                     stats%ensskill%tavg_value_ts,&
                     stats%ensskill%tavg_count_value_ts,dummy1,dummy2)
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
                   if(stats%ensskill%value_final(t,k,l).gt.&
                        LVT_rc%obsCountThreshold) then 
                      merr = sqrt(stats%ensskill%value_final(t,k,l)/&
                              stats%ensskill%value_final(t,k,l))
                      if(merr.ne.0) then 
                         stats%ensskill%value_final(t,k,l) = &
                              1/merr
                      else
                         stats%ensskill%value_final(t,k,l) = LVT_rc%udef
                      endif
                   else
                      stats%ensskill%value_final(t,k,l) = LVT_rc%udef
                   endif
                enddo
                if(metric%computeSC.eq.1) then 
                   do l=1, LVT_rc%nasc
                      if(stats%ensskill%value_asc(t,k,l).gt.&
                           LVT_rc%SCCountThreshold) then 
                         stats%ensskill%value_asc(t,k,l) = &
                              stats%ensskill%value_asc(t,k,l)/ &
                              stats%ensskill%value_asc(t,k,l)
                      else
                         stats%ensskill%value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif
                   
                if(metric%computeADC.eq.1) then 
                   do l=1, LVT_rc%nadc
                      if(stats%ensskill%value_adc(t,k,l).gt.&
                           LVT_rc%ADCCountThreshold) then 
                         stats%ensskill%value_adc(t,k,l) = &
                              stats%ensskill%value_adc(t,k,l)/ &
                              stats%ensskill%value_adc(t,k,l) 
                      else
                         stats%ensskill%value_adc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                endif                   
             enddo
          enddo
             
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%ensskill%value_final(:,k,l),&
                     LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%ensskill%value_ci(k,l))
             enddo
          enddo
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid,stats%ensskill%value_final)
          
          if(metric%computeSC.eq.1) then 
             print*, 'average seasonal cycle computations in ens skill metric need to be implemented'
             stop
!             call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
!                  LVT_rc%ngrid,stats%ensskill%value_asc,stats%ensskill%value_asc)
          endif
          if(metric%computeADC.eq.1) then 
             print*, 'average diurnal cycle computations in ens skill metric need to be implemented'
             stop
!             call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
!                  LVT_rc%ngrid,stats%ensskill%value_adc,stats%ensskill%value_adc)
          endif
       endif
    endif
  end subroutine computeSingleensSkill


  subroutine LVT_writeMetric_ensSkill(pass,final,vlevels,stats,obs)
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

    if(pass.eq.LVT_metrics%ensskill%npass) then 
       if(final.ne.1) then
           if(stats%selectOpt.eq.1) then 
             do k=1,vlevels
                if(LVT_metrics%ensskill%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_gridded(LVT_metrics%ensskill%ftn_ts, &
                           stats%ensskill%value_ts(:,k,l),&
                           stats%vid_ts(LVT_ensSkillid,1),k)
                      call LVT_writevar_gridded(LVT_metrics%ensskill%ftn_ts, &
                           real(stats%ensskill%value_ts(:,k,l)),&
                           stats%vid_count_ts(LVT_ensSkillid,1),&
                           k)
                   enddo
                endif
             enddo
          endif
       else
          if(pass.eq.1) then 
             if(stats%selectOpt.eq.1) then 
                do k=1,vlevels
                   if(LVT_metrics%ensskill%selectOpt.eq.1) then 
                      do l=1,LVT_rc%strat_nlevels                      
                         call LVT_writevar_gridded(&
                              LVT_metrics%ensskill%ftn_total, &
                              stats%ensskill%value_final(:,k,l),&
                              stats%vid_total(LVT_ensSkillid,1),&
                              k)
                         call LVT_writevar_gridded(&
                              LVT_metrics%ensskill%ftn_total, &
                              real(stats%ensskill%value_final(:,k,l)),&
                              stats%vid_count_total(LVT_ensSkillid,1),&
                              k)                         
                      enddo
                   
                      if(LVT_metrics%ensskill%computeSC.eq.1) then 
                         do tind = 1,LVT_rc%nasc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%ensskill%ftn_total,&
                                 stats%ensskill%value_asc(:,k,tind),&
                                 stats%vid_sc_total(tind,LVT_ensSkillid,1),&
                                 1)
                         enddo
                      endif
                      if(LVT_metrics%ensskill%computeADC.eq.1) then 
                         do tind = 1,LVT_rc%nadc
                            call LVT_writevar_gridded(&
                                 LVT_metrics%ensskill%ftn_total,&
                                 stats%ensskill%value_adc(:,k,tind),&
                                 stats%vid_adc_total(tind,LVT_ensSkillid,1),&
                                 1)
                         enddo
                      endif
                      call LVT_writeSummaryStats(&
                           LVT_metrics%ensskill%ftn_summ,&
                           LVT_metrics%ensskill%short_name,&
                           LVT_rc%ngrid,&
                           stats%ensskill%value_final(:,k,:), &
                           stats%ensskill%count_value_final(:,k,:),&
                           stats%standard_name,&
                           stats%ensskill%value_ci(k,:))
                   endif
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_writeMetric_ensSkill

!BOP
! 
! !ROUTINE: LVT_resetMetric_ensSkill)
! \label(LVT_resetMetric_ensSkill)
!
! !INTERFACE:
  subroutine LVT_resetMetric_ensSkill(alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
    logical            :: alarm 
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
             if(LVT_metrics%ensskill%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   stats%ensskill%value_ts(:,k,l) = 0.0
                   stats%ensskill%value_ts(:,k,l)=0 
                   if(alarm) then 
                      stats%ensskill%tavg_value_ts(:,k,l) = 0.0
                      stats%ensskill%tavg_value_ts(:,k,l)=0 
                   endif
                enddo
             endif
          enddo
       endif
    enddo

  end subroutine LVT_resetMetric_ensSkill


!BOP
! 
! !ROUTINE: LVT_writerestart_ensSkill
! 
! !INTERFACE:
  subroutine LVT_writerestart_ensSkill(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for ensSkill metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ensSkill%selectOpt.eq.1) then 
       
       print*, 'The writerestart method is not implemented for ensSkill'

    end if
    
  end subroutine LVT_writerestart_ensSkill

!BOP
! 
! !ROUTINE: LVT_readrestart_ensSkill
! 
! !INTERFACE:
  subroutine LVT_readrestart_ensSkill(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for ensSkill metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ensSkill%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for ensSkill'
       stop
    end if
    
  end subroutine LVT_readrestart_ensSkill

end module LVT_ensSkillMod
