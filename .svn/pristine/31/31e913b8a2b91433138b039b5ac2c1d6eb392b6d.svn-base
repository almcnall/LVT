!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_NSEMod
! \label(LVT_NSEMod)
!
! !INTERFACE:
module LVT_NSEMod
! 
! !USES:   
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
  use LVT_TSMod
  use LVT_logMod
  use LVT_CIMod

  private
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This module handles the Nash-SutCliffe-Efficiency (NSE)
!   computations by comparing the LIS output to the specified 
!   observations. 
!   
!   NSE is defined as : 
!    
!    NSE = (1-(sum((obs-model)^2)/sum((obs-mean(obs))^2)))
!
!  NOTES: Stratification by seasons is not supported for NSE computation
!  TODO: add codes for stratification, confidence intervals
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 2 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initNSE
  public :: LVT_diagnoseNSE
  public :: LVT_computeNSE
  public :: LVT_writeMetric_NSE
  public :: LVT_resetMetric_NSE
  public :: LVT_writerestart_NSE
  public :: LVT_readrestart_NSE


contains
  
!BOP
! 
! !ROUTINE: LVT_initNSE
! \label{LVT_initNSE}
!
! !INTERFACE: 
  subroutine LVT_initNSE(model,obs,stats,metric)
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
       allocate(stats%nse%value(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))
       stats%nse%value = 0.0
       allocate(stats%nse%value_numer(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       
       stats%nse%value_numer = 0
       allocate(stats%nse%value_denom(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       
       stats%nse%value_denom = 0
       allocate(stats%nse%value_obs_mean(LVT_rc%ngrid,model%selectNlevs,&
            LVT_rc%strat_nlevels))
       allocate(stats%nse%count_value_obs_mean(LVT_rc%ngrid,model%selectNlevs,&
            LVT_rc%strat_nlevels))
       stats%nse%value_obs_mean = 0 
       stats%nse%count_value_obs_mean = 0 

       allocate(stats%nse%value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%nse%value_ci = LVT_rc%udef
    endif

    if(metric%timeopt.eq.1) then 
       allocate(stats%nse%value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%nse%value_ts = 0.0

       allocate(stats%nse%tavg_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%nse%tavg_value_ts = 0.0
       allocate(stats%nse%tavg_count_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%strat_nlevels))
       stats%nse%tavg_count_value_ts = 0

       allocate(stats%nse%value_numer_ts(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       
       stats%nse%value_numer_ts = 0
       allocate(stats%nse%value_denom_ts(LVT_rc%ngrid, model%selectNlevs,&
            LVT_rc%strat_nlevels))       
       stats%nse%value_denom_ts = 0
       allocate(stats%nse%value_obs_mean_ts(LVT_rc%ngrid,model%selectNlevs,&
            LVT_rc%strat_nlevels))
       allocate(stats%nse%count_value_obs_mean_ts(LVT_rc%ngrid,model%selectNlevs,&
            LVT_rc%strat_nlevels))
       stats%nse%value_obs_mean_ts = 0 
       stats%nse%count_value_obs_mean_ts = 0 

    end if

!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 2
    metric%obsData = .false. 
    metric%stdevFlag = .false. 

  end subroutine LVT_initNSE

!BOP
! 
! !ROUTINE: LVT_diagnoseNSE
! \label{LVT_diagnoseNSE}
!
! !INTERFACE: 
  subroutine LVT_diagnoseNSE(pass)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the NSE calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleNSE](\ref{diagnoseSingleNSE})
!     updates the anomaly correlation computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer :: pass
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(LVT_metrics%nse%selectOpt.eq.1.or.&
         LVT_metrics%nse%timeOpt.eq.1) then 
       if(pass.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleNSEobsMean(&
                  obs,model,stats)
             
             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
       elseif(pass.eq.2) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)
          
          do while(associated(model))
             call diagnoseSingleNSE(obs,model,stats,&
                  LVT_metrics%nse)
             
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif

  end subroutine LVT_diagnoseNSE

!BOP
! 
! !ROUTINE: diagnoseSingleNSEobsMean
! \label{diagnoseSingleNSEobsMean}
!
! !INTERFACE: 
  subroutine diagnoseSingleNSEobsMean(obs, model, stats)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the computation of the 
!   observation mean for a specified variable. 
!
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
!BOP
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry) :: stats
!EOP
    integer    :: t,k
    integer    :: c,r

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             if(trim(obs%units).eq.trim(model%units)) then                 
                if(obs%count(t,k).ne.0) then  
                   stats%nse%value_obs_mean(t,k,1) = &
                        stats%nse%value_obs_mean(t,k,1) + &
                        obs%value(t,k)
                   stats%nse%count_value_obs_mean(t,k,1) = &
                        stats%nse%count_value_obs_mean(t,k,1) + 1
                   if(LVT_rc%strat_nlevels.gt.1) then 
                      if(LVT_stats%strat_var(t,k).gt.&
                           LVT_rc%strat_var_threshold) then 
                         stats%nse%value_obs_mean(t,k,2) = &
                              stats%nse%value_obs_mean(t,k,2) + &
                              obs%value(t,k)
                         stats%nse%count_value_obs_mean(t,k,2) = &
                              stats%nse%count_value_obs_mean(t,k,2) + 1
                      elseif(LVT_stats%strat_var(t,k).le.&
                           LVT_rc%strat_var_threshold) then 
                         stats%nse%value_obs_mean(t,k,3) = &
                              stats%nse%value_obs_mean(t,k,3) + &
                              obs%value(t,k)
                         stats%nse%count_value_obs_mean(t,k,3) = &
                              stats%nse%count_value_obs_mean(t,k,3) + 1
                      endif
                   endif
                   if(LVT_metrics%nse%timeOpt.eq.1) then 
                      stats%nse%value_obs_mean_ts(t,k,1) = &
                           stats%nse%value_obs_mean_ts(t,k,1) + &
                           obs%value(t,k)
                      stats%nse%count_value_obs_mean_ts(t,k,1) = &
                           stats%nse%count_value_obs_mean_ts(t,k,1) + 1
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then 
                            stats%nse%value_obs_mean_ts(t,k,2) = &
                                 stats%nse%value_obs_mean_ts(t,k,2) + &
                                 obs%value(t,k)
                            stats%nse%count_value_obs_mean_ts(t,k,2) = &
                                 stats%nse%count_value_obs_mean_ts(t,k,2) + 1
                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then 
                            stats%nse%value_obs_mean_ts(t,k,3) = &
                                 stats%nse%value_obs_mean_ts(t,k,3) + &
                                 obs%value(t,k)
                            stats%nse%count_value_obs_mean_ts(t,k,3) = &
                                 stats%nse%count_value_obs_mean_ts(t,k,3) + 1
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
    
  end subroutine diagnoseSingleNSEobsMean
!BOP
! 
! !ROUTINE: diagnoseSingleNSE
! \label{diagnoseSingleNSE}
!
! !INTERFACE: 
  subroutine diagnoseSingleNSE(obs, model, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the NSE computation for a single variable. 
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
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer    :: t,k,l
    integer    :: c,r

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
       do t=1, LVT_rc%ngrid
          do k=1,model%selectNlevs
             if(trim(obs%units).eq.trim(model%units)) then 
                if(model%count(t,k).ne.0.and.&
                     obs%count(t,k).ne.0) then 
                   if(metric%selectOpt.eq.1) then 
                      if(stats%nse%count_value_obs_mean(t,k,1).ne.0) then 
                         stats%nse%value_numer(t,k,1) = stats%nse%value_numer(t,k,1)+&
                              (obs%value(t,k)-model%value(t,k))**2
                         stats%nse%value_denom(t,k,1) = stats%nse%value_denom(t,k,1)+&
                              (obs%value(t,k)-stats%nse%value_obs_mean(t,k,1))**2
                      endif
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then 
                            stats%nse%value_numer(t,k,2) = stats%nse%value_numer(t,k,2)+&
                                 (obs%value(t,k)-model%value(t,k))**2
                            stats%nse%value_denom(t,k,2) = stats%nse%value_denom(t,k,2)+&
                                 (obs%value(t,k)-stats%nse%value_obs_mean(t,k,2))**2
                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then 
                            stats%nse%value_numer(t,k,3) = stats%nse%value_numer(t,k,3)+&
                                 (obs%value(t,k)-model%value(t,k))**2
                            stats%nse%value_denom(t,k,3) = stats%nse%value_denom(t,k,3)+&
                                 (obs%value(t,k)-stats%nse%value_obs_mean(t,k,3))**2
                         endif
                      endif
                   endif
                   if(metric%timeOpt.eq.1) then 
                      if(stats%nse%count_value_obs_mean_ts(t,k,1).ne.0) then 
                         stats%nse%value_numer_ts(t,k,1) = stats%nse%value_numer_ts(t,k,1)+&
                              (obs%value(t,k)-model%value(t,k))**2
                         stats%nse%value_denom_ts(t,k,1) = stats%nse%value_denom_ts(t,k,1)+&
                              (obs%value(t,k)-stats%nse%value_obs_mean(t,k,1))**2
                      endif
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then 
                            stats%nse%value_numer_ts(t,k,2) = stats%nse%value_numer_ts(t,k,2)+&
                                 (obs%value(t,k)-model%value(t,k))**2
                            stats%nse%value_denom_ts(t,k,2) = stats%nse%value_denom_ts(t,k,2)+&
                                 (obs%value(t,k)-stats%nse%value_obs_mean(t,k,2))**2
                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then 
                            stats%nse%value_numer_ts(t,k,3) = stats%nse%value_numer_ts(t,k,3)+&
                                 (obs%value(t,k)-model%value(t,k))**2
                            stats%nse%value_denom_ts(t,k,3) = stats%nse%value_denom_ts(t,k,3)+&
                                 (obs%value(t,k)-stats%nse%value_obs_mean(t,k,3))**2
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
  end subroutine diagnoseSingleNSE


!BOP
! 
! !ROUTINE: LVT_computeNSE
! \label{LVT_computeNSE}
!
! !INTERFACE: 
  subroutine LVT_computeNSE(pass,alarm)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine invokes the method to compute NSE values, 
!  for each specified variable. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    
    integer :: pass
    logical :: alarm
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats
    integer                          :: i 

    if(LVT_metrics%nse%selectOpt.eq.1) then
       if(pass.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call computeSingleNSEobsMean(alarm, obs,model,stats)

             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
       elseif(pass.eq.2) then 
          if(LVT_metrics%nse%selectOpt.eq.1.or.&
               LVT_metrics%nse%timeOpt.eq.1) then 
             if(alarm) then 
                if(LVT_metrics%nse%timeOpt.eq.1.and.&
                     LVT_metrics%nse%extractTS.eq.1) then 
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%nse%ftn_ts_loc(i),200,advance='no') &
                           LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                           LVT_rc%hr,'',LVT_rc%mn, '' 
                   enddo
                endif
             endif
200          format(I4, a1, I2.2, a1, I2.2, a1, I2.2, a1, I2.2,a1)
             
             call LVT_getDataStream1Ptr(model)
             call LVT_getDataStream2Ptr(obs)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(model))
                
                call computeSingleNSE(alarm,obs,model,stats,&
                     LVT_metrics%nse)
                
                model => model%next
                obs => obs%next
                stats => stats%next
             enddo
             
             if(alarm) then 
                if(LVT_metrics%nse%timeOpt.eq.1.and.&
                     LVT_metrics%nse%extractTS.eq.1) then 
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%nse%ftn_ts_loc(i),fmt='(a1)') ''
                   enddo
                endif
             end if
          endif
       endif
    endif
  end subroutine LVT_computeNSE


!BOP
! 
! !ROUTINE: computeSingleNSEobsMean
! \label{computeSingleNSEobsMean}
!
! !INTERFACE: 
  subroutine computeSingleNSEobsMean(alarm, obs, model, stats)
! 
! !USES: 
    use LVT_coreMod,    only : LVT_rc, LVT_domain
    use LVT_histDataMod 
    use LVT_statsDataMod
    use LVT_logMod,     only : LVT_logunit

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the mean observation value for a specified variable
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
!BOP
! !ARGUMENTS: 
    logical                 :: alarm
    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry) :: stats
!EOP
    integer    :: t,k,l,m

    integer    :: c,r

    if(alarm) then 
       if(LVT_metrics%nse%timeOpt.eq.1.and.&
            stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1, LVT_rc%strat_nlevels
                   if(stats%nse%count_value_obs_mean_ts(t,k,l).ne.0) then 
                      stats%nse%value_obs_mean_ts(t,k,l) = &
                           stats%nse%value_obs_mean_ts(t,k,l) /&
                           stats%nse%count_value_obs_mean_ts(t,k,l)
                   endif
                enddo
             enddo
          enddo
       endif
    endif

    if(LVT_rc%endtime.eq.1) then 
       write(LVT_logunit,*) 'Computing NSE obs Mean '
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1, LVT_rc%strat_nlevels
                   if(stats%nse%count_value_obs_mean(t,k,l).ne.0) then 
                      stats%nse%value_obs_mean(t,k,l) = &
                           stats%nse%value_obs_mean(t,k,l) /&
                           stats%nse%count_value_obs_mean(t,k,l)
                   endif
                enddo
             enddo
          enddo
       endif
    endif
    
  end subroutine computeSingleNSEobsMean

!BOP
! 
! !ROUTINE: computeSingleNSE
! \label{computeSingleNSE}
!
! !INTERFACE: 
  subroutine computeSingleNSE(alarm,obs, model, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the NSE values for a single variable
!  The arguments are: 
!
!  \begin{description}
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
    real     :: numer,denom

    if(metric%timeOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%nse%count_value_obs_mean_ts(t,k,l).ne.0) then 
                      numer = stats%nse%value_numer_ts(t,k,l)
                      denom = stats%nse%value_denom_ts(t,k,l)
                      if(denom.ne.0) then 
                         stats%nse%value_ts(t,k,l) = 1 - numer/denom
                         
                         stats%nse%tavg_value_ts(t,k,l) = & 
                              stats%nse%tavg_value_ts(t,k,l) + & 
                              stats%nse%value_ts(t,k,l)
                         stats%nse%tavg_count_value_ts(t,k,l) = & 
                              stats%nse%tavg_count_value_ts(t,k,l) + 1

                      else
                         stats%nse%value_ts(t,k,l) = LVT_rc%udef
                      endif
                   endif
                enddo
             enddo
          enddo
          if(alarm) then 
             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%nse%tavg_count_value_ts(t,k,l).gt.0) then 
                         stats%nse%tavg_value_ts(t,k,l) = & 
                              stats%nse%tavg_value_ts(t,k,l)/&
                              stats%nse%tavg_count_value_ts(t,k,l)
                      endif
                   enddo
                enddo
             enddo
             call LVT_writeTSinfo(metric%ftn_ts_loc,&
                  model,&
                  LVT_rc%ngrid,&
                  stats%nse%tavg_value_ts,&
                  stats%nse%tavg_count_value_ts)
          endif
       endif
    endif

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%nse%count_value_obs_mean(t,k,l).ne.0) then 
                      numer = stats%nse%value_numer(t,k,l)
                      denom = stats%nse%value_denom(t,k,l)
                      if(denom.ne.0) then 
                         stats%nse%value(t,k,l) = 1 - numer/denom
                      else
                         stats%nse%value(t,k,l) = LVT_rc%udef
                      endif
                   endif
                enddo
             enddo
          enddo
          
          do k=1,model%selectNlevs
             do l=1,LVT_rc%strat_nlevels
                call LVT_computeCI(stats%nse%value(:,k,l),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%nse%value_ci(k,l))
             enddo
          enddo
          
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid,stats%nse%value)      
       endif
    endif

  end subroutine computeSingleNSE

!BOP
! 
! !ROUTINE: LVT_writeMetric_NSE
! \label{LVT_writeMetric_NSE}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_NSE(pass,final,vlevels,stats,obs)
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
    integer                 :: l
    integer                 :: k

    if(LVT_metrics%nse%selectOpt.eq.1) then
       if(pass.eq.LVT_metrics%nse%npass) then 
          if(stats%selectOpt.eq.1) then 
             do k=1,vlevels
                do l=1,LVT_rc%strat_nlevels
                   call LVT_writevar_gridded(LVT_metrics%nse%ftn_total, &
                        stats%nse%value(:,k,l),stats%vid_total(LVT_NSEid,1))
                   call LVT_writevar_gridded(LVT_metrics%nse%ftn_total, &
                        real(stats%nse%count_value_obs_mean(:,k,l)),&
                        stats%vid_count_total(LVT_NSEid,1))                
                enddo
                call LVT_writeSummaryStats(&
                     LVT_metrics%nse%ftn_summ,&
                     LVT_metrics%nse%short_name,&
                     LVT_rc%ngrid,&
                     stats%nse%value(:,k,:), &
                     stats%nse%count_value_obs_mean(:,k,:),stats%standard_name,&
                     stats%nse%value_ci(k,:))
             enddo
          endif
       endif
    endif

  end subroutine LVT_writeMetric_NSE


  subroutine LVT_resetMetric_NSE(alarm)

    logical                 :: alarm

    integer                 :: vlevels
    
    integer                :: i,k,l
    type(LVT_metadataEntry), pointer :: model
    type(LVT_statsEntry)   , pointer :: stats


    call LVT_getDataStream1Ptr(model)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       
       if(stats%selectOpt.eq.1) then 
          do k=1,model%selectNlevs
             if(LVT_metrics%nse%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   stats%nse%value_ts(:,k,l) = 0
                   stats%nse%value_numer_ts(:,k,l) = 0
                   stats%nse%value_denom_ts(:,k,l) = 0
                   stats%nse%value_obs_mean_ts(:,k,l) = 0 
                   stats%nse%count_value_obs_mean_ts(:,k,l) = 0 
                   if(alarm) then 
                      stats%nse%tavg_value_ts(:,k,l) = 0 
                      stats%nse%tavg_count_value_ts(:,k,l) = 0 
                   endif
                enddo
             endif
             
          enddo
          
       endif
       model => model%next
       stats => stats%next

    enddo


  end subroutine LVT_resetMetric_NSE


!BOP
! 
! !ROUTINE: LVT_writerestart_NSE
! 
! !INTERFACE:
  subroutine LVT_writerestart_NSE(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for NSE metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%nse%selectOpt.eq.1) then 
       
       print*, 'WARNING: The writerestart method is not implemented for NSE'

    end if
    
  end subroutine LVT_writerestart_NSE

!BOP
! 
! !ROUTINE: LVT_readrestart_NSE
! 
! !INTERFACE:
  subroutine LVT_readrestart_NSE(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for NSE metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%nse%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for NSE'
       stop
    end if
    
  end subroutine LVT_readrestart_NSE

end module LVT_NSEMod
