!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_KStestMod
! \label(LVT_KStestMod)
!
! !INTERFACE:
module LVT_KStestMod
! 
! !USES:
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
  use LVT_TSMod
  use LVT_logMod
  use LVT_CIMod
  use LVT_NumericalRecipesMod

  implicit none

  private
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
!  !DESCRIPTION: 
!   This module performs the Kolmogorov Smirnov (K-S) test to 
!   compare the probability distributions of model and observational 
!   data. 
! 
! !FILES USED:
!
!  !REVISION HISTORY: 
!  20 Nov 2014    Sujay Kumar  Initial Specification
! 
!EOP


!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initKStest
  public :: LVT_diagnoseKStest
  public :: LVT_computeKStest 
  public :: LVT_writeMetric_KStest
  public :: LVT_resetMetric_KStest
  public :: LVT_writerestart_KStest
  public :: LVT_readrestart_KStest

  integer :: CDF_nbins = 100
contains
  
!BOP
! 
! !ROUTINE: LVT_initKStest
! \label{LVT_initKStest}
!
! !INTERFACE: 
  subroutine LVT_initKStest(model,obs,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes the data structures required for 
!  K-S test calculations 
!  
!  The arguments are:  
!  \begin{description}
!   \item[model]  model variable object
!   \item[obs]    observation object
!   \item[stats]  object to hold the updated statistics
!   \item[metric] object representing the K-S test outputs
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
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!EOP
    if(metric%selectOpt.eq.1) then 
       allocate(stats%kstest%value_model_max(LVT_rc%ngrid, model%selectNlevs))
       allocate(stats%kstest%value_model_min(LVT_rc%ngrid, model%selectNlevs))
       allocate(stats%kstest%value_model_bincounts(LVT_rc%ngrid, model%selectNlevs, &
            CDF_nbins))
       allocate(stats%kstest%value_model_xrange(LVT_rc%ngrid, model%selectNlevs, &
            CDF_nbins))
       allocate(stats%kstest%value_model_cdf(LVT_rc%ngrid, model%selectNlevs, &
            CDF_nbins))

       allocate(stats%kstest%value_obs_max(LVT_rc%ngrid, obs%vlevels))
       allocate(stats%kstest%value_obs_min(LVT_rc%ngrid, obs%vlevels))
       allocate(stats%kstest%value_obs_bincounts(LVT_rc%ngrid, obs%vlevels, &
            CDF_nbins))
       allocate(stats%kstest%value_obs_xrange(LVT_rc%ngrid, obs%vlevels, &
            CDF_nbins))
       allocate(stats%kstest%value_obs_cdf(LVT_rc%ngrid, obs%vlevels, &
            CDF_nbins))

       allocate(stats%kstest%value_d(LVT_rc%ngrid, obs%vlevels))
       allocate(stats%kstest%value_prob(LVT_rc%ngrid, obs%vlevels))
       
       stats%kstest%value_model_max       = -1000000 
       stats%kstest%value_model_min       = 1000000 
       stats%kstest%value_model_bincounts = 0 
       stats%kstest%value_model_xrange    = 0 
       stats%kstest%value_model_cdf       = 0 

       stats%kstest%value_obs_max       = -1000000 
       stats%kstest%value_obs_min       = 1000000 
       stats%kstest%value_obs_bincounts = 0 
       stats%kstest%value_obs_xrange    = 0 
       stats%kstest%value_obs_cdf       = 0 

       stats%kstest%value_d    = LVT_rc%udef
       stats%kstest%value_prob = LVT_rc%udef

    endif
    
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 2
    metric%obsData = .false. 
    metric%stdevFlag = .false. 
    metric%customNames = .true. 
    metric%nfields = 2
    allocate(metric%mName(metric%nfields))

    metric%mName(1) = 'KS_D'
    metric%mName(2) = 'KS_prob'

  end subroutine LVT_initKStest
  
!BOP
! 
! !ROUTINE: LVT_diagnoseKStest
! \label{LVT_diagnoseKStest}
!
! !INTERFACE: 
  subroutine LVT_diagnoseKStest(pass)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the K-S test 
!   calculation for desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleDrange](\ref{diagnoseSingleDrange})
!     updates the dynamic range computation for a single variable 
!    \item[diagnoseSingleCDF](\ref{diagnoseSingleCDF})
!     updates the Cumulative Distribution Function (CDF) computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none
    integer                 :: pass
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%kstest%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleDrange(&
                  obs, model, stats)

             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo
       endif
    elseif(pass.eq.2) then 
       if(LVT_metrics%kstest%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call diagnoseSingleCDF(&
                  obs, model, stats, &
                  LVT_metrics%kstest)

             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseKStest

!BOP
! 
! !ROUTINE: diagnoseSingleDrange
! \label{diagnoseSingleDrange}
!
! !INTERFACE: 
  subroutine diagnoseSingleDrange(obs, model, stats)
! 
! !USES: 
    use LVT_coreMod,    only : LVT_rc, LVT_domain
    use LVT_histDataMod 
    use LVT_statsDataMod
    use LVT_logMod,     only : LVT_logunit, LVT_endrun

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the dynamic range (max/min) computation of the 
!   specified variable. 
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
!
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry) :: stats
! 
!EOP
    integer    :: t,k

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             if(trim(obs%units).eq.trim(model%units)) then 
                if(model%count(t,k).ne.0) then 
                   if(model%value(t,k).gt.stats%kstest%value_model_max(t,k)) then 
                      stats%kstest%value_model_max(t,k) = &
                           model%value(t,k)
                   endif
                   if(model%value(t,k).lt.stats%kstest%value_model_min(t,k)) then 
                      stats%kstest%value_model_min(t,k) = &
                           model%value(t,k)
                   endif
                endif
                if(obs%count(t,k).ne.0) then  
                   if(obs%value(t,k).gt.stats%kstest%value_obs_max(t,k)) then 
                      stats%kstest%value_obs_max(t,k) = &
                           obs%value(t,k)
                   endif
                   if(obs%value(t,k).lt.stats%kstest%value_obs_min(t,k)) then 
                      stats%kstest%value_obs_min(t,k) = &
                           obs%value(t,k)
                   endif

                endif
                if(LVT_rc%strat_nlevels.gt.1) then 
                   if(LVT_stats%strat_var(t,k).gt.&
                        LVT_rc%strat_var_threshold) then 
                      !not implemented yet
                   elseif(LVT_stats%strat_var(t,k).le.&
                        LVT_rc%strat_var_threshold) then 
                      !not implemented yet
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

  end subroutine diagnoseSingleDrange

!BOP
! 
! !ROUTINE: diagnoseSingleCDF
! \label{diagnoseSingleCDF}
!
! !INTERFACE: 
  subroutine diagnoseSingleCDF(obs, model, stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the Cumulative Distribution Function (CDF)
!  calculations
!
!  The arguments are: 
!  \begin{description}
!   \item[obs] observation object
!   \item[model] model variable object
!   \item[stats] object to hold the updated statistics
!   \item[metric] object representing K-S test metric
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none

    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry) :: stats
    type(LVT_metricEntry)   :: metric

    integer    :: t,k, binval

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             if(trim(obs%units).eq.trim(model%units)) then 
                if(model%count(t,k).ne.0) then 
                   binval = nint((model%value(t,k) - &
                        stats%kstest%value_model_min(t,k))*(CDF_nbins-1)/&
                        (stats%kstest%value_model_max(t,k) - & 
                        stats%kstest%value_model_min(t,k))) + 1
                   stats%kstest%value_model_bincounts(t,k,binval) = & 
                        stats%kstest%value_model_bincounts(t,k,binval) + 1
                   
                endif
                if(obs%count(t,k).ne.0) then 
                   binval = nint((obs%value(t,k) - &
                        stats%kstest%value_obs_min(t,k))*(CDF_nbins-1)/&
                        (stats%kstest%value_obs_max(t,k) - & 
                        stats%kstest%value_obs_min(t,k))) + 1
                   stats%kstest%value_obs_bincounts(t,k,binval) = & 
                        stats%kstest%value_obs_bincounts(t,k,binval) + 1
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
    
  end subroutine diagnoseSingleCDF


!BOP
! 
! !ROUTINE: LVT_computeKStest
! \label{LVT_computeKStest}
!
! !INTERFACE: 
  subroutine LVT_computeKStest(pass,alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to perform the K-S test 
!   for the desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleCDF](\ref{computeSingleCDF})
!     updates the CDF computation for a single variable and 
!     perform the K-S test
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[pass]
!     current pass number over the data points
!    \item[alarm]
!     flag indicating if the temporal output interval has been reached
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS:
    implicit none

    integer     :: pass
    logical     :: alarm

! 
!EOP
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(LVT_metrics%kstest%selectOpt.eq.1) then 
       if(pass.eq.2) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call computeSingleCDF(&
                  obs, model, stats)
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif

  end subroutine LVT_ComputeKStest

!BOP
! 
! !ROUTINE: computeSingleCDF
! \label{computeSingleCDF}
!
! !INTERFACE: 
  subroutine computeSingleCDF(obs, model, stats)
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
!  This routine computes the CDF for a specified variable
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
    integer    :: t,k,m,l,i
    real       :: delta

    if(LVT_rc%endtime.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                if(sum(stats%kstest%value_model_bincounts(t,k,:)).gt.0) then 
                   stats%kstest%value_model_cdf(t,k,1) = & 
                        stats%kstest%value_model_bincounts(t,k,1)/&
                        sum(stats%kstest%value_model_bincounts(t,k,:))
                  
                   stats%kstest%value_model_xrange(t,k,1) = &
                        stats%kstest%value_model_min(t,k)
                   delta  = (stats%kstest%value_model_max(t,k) - & 
                        stats%kstest%value_model_min(t,k))/float(CDF_nbins)

                   do i=2,CDF_nbins
                      stats%kstest%value_model_xrange(t,k,i) = &
                           stats%kstest%value_model_xrange(t,k,i-1) +delta

                      stats%kstest%value_model_cdf(t,k,i) = & 
                           stats%kstest%value_model_cdf(t,k,i-1) + & 
                           stats%kstest%value_model_bincounts(t,k,i)/&
                           sum(stats%kstest%value_model_bincounts(t,k,:))
                   enddo
                endif
                
                if(sum(stats%kstest%value_obs_bincounts(t,k,:)).gt.0) then 
                   stats%kstest%value_obs_cdf(t,k,1) = & 
                        stats%kstest%value_obs_bincounts(t,k,1)/&
                        sum(stats%kstest%value_obs_bincounts(t,k,:))
                   
                   stats%kstest%value_obs_xrange(t,k,1) = &
                        stats%kstest%value_obs_min(t,k)
                   delta  = (stats%kstest%value_obs_max(t,k) - & 
                        stats%kstest%value_obs_min(t,k))/float(CDF_nbins)
                   
                   do i=2,CDF_nbins
                      stats%kstest%value_obs_xrange(t,k,i) = &
                           stats%kstest%value_obs_xrange(t,k,i-1) +delta

                      stats%kstest%value_obs_cdf(t,k,i) = & 
                           stats%kstest%value_obs_cdf(t,k,i-1) + & 
                           stats%kstest%value_obs_bincounts(t,k,i)/&
                           sum(stats%kstest%value_obs_bincounts(t,k,:))
                   enddo
                   
                endif
             enddo
          enddo
          
          do t=1,LVT_rc%ngrid
             do k=1, obs%vlevels
                !             do i=1,CDF_nbins
                !                print*, stats%kstest%value_model_xrange(t,k,i), &
                !                     stats%kstest%value_model_cdf(t,k,i),&
                !                     stats%kstest%value_obs_xrange(t,k,i),&
                !                     stats%kstest%value_obs_cdf(t,k,i)
                !             enddo
                call ks2d2s(stats%kstest%value_model_xrange(t,k,:), &
                     stats%kstest%value_model_cdf(t,k,:),&
                     CDF_nbins, &
                     stats%kstest%value_obs_xrange(t,k,:),&
                     stats%kstest%value_obs_cdf(t,k,:),&
                     CDF_nbins, &
                     stats%kstest%value_d(t,k),&
                     stats%kstest%value_prob(t,k))
             enddo
          enddo
       endif
    endif
    
  end subroutine computeSingleCDF


!BOP
! 
! !ROUTINE: LVT_writeMetric_KStest
! \label{LVT_writeMetric_KStest}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_KStest(pass,final,vlevels,stats,obs)
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
!   This subroutine writes the computed anomaly correlation values to an 
!   external file. 
! 
!  \begin{description}
!   \item[pass]  current pass number over the data points
!   \item[final] boolean value indicating if computation is for the final
!                timestep (=1) or not (=0)
!   \item[vlevels] number of vertical levels of the variable being written 
!   \item[stats] object to hold the updated statistics
!  \end{description}
!   The methods invoked are: 
!   \begin{description}
!    \item[LVT\_writevar\_gridded](\ref{LVT_writevar_gridded})
!    writes the variable to a gridded output file. 
!    \item[LVT\_writeSummaryStats](\ref{LVT_writeSummaryStats})
!    writes the domain averaged summary statistics to a text file
!   \end{description}
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

! The D statistic is written to the variable field and
! the probability is written to the count field
    if(pass.eq.LVT_metrics%kstest%npass) then 
       if(LVT_metrics%kstest%selectOpt.eq.1) then
          if(stats%selectOpt.eq.1) then 
             do k=1,vlevels
                call LVT_writevar_gridded(LVT_metrics%kstest%ftn_total, &
                     stats%kstest%value_d(:,k),stats%vid_total(LVT_Kstestid,1),k)
                call LVT_writevar_gridded(LVT_metrics%kstest%ftn_total, &
                     stats%kstest%value_prob(:,k),&
                     stats%vid_count_total(LVT_Kstestid,1),k)                
             enddo
          endif
       endif
    endif

  end subroutine LVT_writeMetric_KStest

!BOP
! 
! !ROUTINE: LVT_resetMetric_KStest
! \label{LVT_resetMetric_KStest}
!
! !INTERFACE:
  subroutine LVT_resetMetric_KStest()
! 
! !USES:
    use LVT_coreMod,   only : LVT_rc
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This subroutine resets the relevant variables during the anomaly 
!  correlation computation between each temporal output times
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

!nothing to do
  end subroutine LVT_resetMetric_KStest

!BOP
! 
! !ROUTINE: LVT_writerestart_KStest
! 
! !INTERFACE:
  subroutine LVT_writerestart_KStest(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for KStest metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    integer              :: k,m,l
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

  end subroutine LVT_writerestart_KStest

!BOP
! 
! !ROUTINE: LVT_readrestart_KStest
! 
! !INTERFACE:
  subroutine LVT_readrestart_KStest(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for KStest metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    integer              :: k,m,l
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

  end subroutine LVT_readrestart_KStest


end module LVT_KStestMod
