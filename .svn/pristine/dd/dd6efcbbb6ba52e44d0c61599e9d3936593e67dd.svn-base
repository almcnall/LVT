!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
module LVT_MetricEntropyMod
!
!BOP
! !MODULE: LVT_MetricEntropyMod
! 
!  !DESCRIPTION: 
!   This module handles the computation of the information-theory
!   metric called 'metric entropy'. The program works by converting
!   time series of variables to symbolic string sequences (binary
!   sequences). 
!  
!   References: 
!   Shannon, A mathematical theory of communication, 1984
!   Pachepsky et al., Geoderma, 134, 253-266, 2006. 
!
!  !NOTES: 
!  The computation of the metric in time, averge seasonal cycles, 
!  average diurnal cycles are not available for this metric. 
!
!  !REVISION HISTORY: 
!  1 Aug 2011    Sujay Kumar  Initial Specification
!
!EOP
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
  use LVT_TSMod
  use LVT_logMod
  use LVT_CIMod

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initMetricEntropy
  public :: LVT_diagnoseMetricEntropy
  public :: LVT_computeMetricEntropy 
  public :: LVT_writeMetric_MetricEntropy
  public :: LVT_resetMetric_MetricEntropy
  public :: LVT_writerestart_MetricEntropy
  public :: LVT_readrestart_MetricEntropy

contains
  
!BOP
! !ROUTINE: LVT_initMetricEntropy 
! \label{LVT_initMetricEntropy}
!
! !INTERFACE: 
  subroutine LVT_initMetricEntropy(model,obs,stats,metric)
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
!  This is the initialization routine for the metric entropy calculations. 
!  This routine allocates memory for the required data structures
!  and sets the number of required passes through the data (=1).
!  
!  The subroutine arguments are: 
!  \begin{description}
!    \item[model]
!      object to hold model variable information    
!    \item[obs]
!      object to hold observation variable information
!    \item[stats]
!     object to hold statistics computation related information
!    \item[metric]
!     object to hold metric-specific information
!   \end{description}
!EOP    
    if(metric%selectOpt.eq.1) then 
       allocate(stats%mentropy%value_model_ts(LVT_rc%ngrid,model%selectNlevs,&
            LVT_rc%nts))
       allocate(stats%mentropy%value_final(LVT_rc%ngrid,model%selectNlevs))
       allocate(stats%mentropy%value_ci(model%selectNlevs))

       stats%mentropy%value_model_ts = LVT_rc%udef
       stats%mentropy%value_final = LVT_rc%udef
       stats%mentropy%value_ci = LVT_rc%udef
    endif

!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------
    metric%npass = 1    
    metric%obsData = .true. 
    metric%stdevFlag = .false. 

!-------------------------------------------------------------------------
! These options are not supported
!-------------------------------------------------------------------------
    metric%timeOpt = 0 
    metric%extractTS = 0 
    metric%computeSC = 0 
    metric%computeADC = 0 


  end subroutine LVT_initMetricEntropy

!BOP
! !ROUTINE: LVT_diagnoseMetricEntropy
! \label{LVT_diagnoseMetricEntropy}
! 
! !INTERFACE: 
  subroutine LVT_diagnoseMetricEntropy(pass)

    implicit none
! !ARGUMENTS:
    integer                 :: pass
! !DESCRIPTION:
!   This routine invokes the call to update the metric entropy calculation
!   for each variable. 
! 
!   The arguments are: 
!   \begin{description}
!    \item[pass]
!     current pass index number over the data points
!   \end{description}
!EOP

    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%mentropy%selectOpt.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleMetricEntropy(&
                  obs,model,stats,&
                  LVT_metrics%mentropy)

             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
       endif
    endif

  end subroutine LVT_diagnoseMetricEntropy

!BOP
! !ROUTINE: diagnoseSingleMetricEntropy
! \label{diagnoseSingleMetricEntropy}
! 
! !INTERFACE: 
  subroutine diagnoseSingleMetricEntropy(obs,model,stats,metric)
! !USES: 
    use LVT_coreMod,    only : LVT_rc, LVT_domain
    use LVT_timeMgrMod, only : LVT_clock, LVT_calendar
    use LVT_histDataMod 
    use LVT_statsDataMod
    use LVT_logMod,     only : LVT_logunit, LVT_endrun

    implicit none
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
! !DESCRIPTION:  
!  This subroutine gathers the relevant information for the computation
!  of metric entropy during the pass through the data. The routine
!  simply stores the variable values into an array for computing
!  the metric at the end of the analysis time. 
!
!  The subroutine arguments are: 
!  \begin{description}
!    \item[model]
!      object to hold model variable information    
!    \item[obs]
!      object to hold observation variable information
!    \item[stats]
!     object to hold statistics computation related information
!    \item[metric]
!     object to hold metric-specific information
!   \end{description}
!EOP 
    type(ESMF_Time)         :: currTime,startTime
    type(ESMF_TimeInterval) :: timestep
    integer                 :: tindex
    integer                 :: status
    integer                 :: t,k
    
    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             if(model%count(t,k).gt.0) then
                if(metric%selectOpt.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then

                      call ESMF_TimeSet(currTime, yy = LVT_rc%yr, &
                           mm = LVT_rc%mo, &
                           dd = LVT_rc%da, &
                           h  = LVT_rc%hr, &
                           m  = LVT_rc%mn, & 
                           s  = LVT_rc%ss, &
                           calendar = LVT_calendar, & 
                           rc = status)
                      call LVT_verify(status,&
                           'Error in ESMF_TimeSet diagnoseEntropy')
                      call ESMF_ClockGet(LVT_clock,&
                           startTime = starttime,&
                           rc=status)
                      call LVT_verify(status,&
                           'Error in ESMF_TimeGet diagnoseEntropy')
                      call ESMF_TimeIntervalSet(timestep,&
                           s=LVT_rc%ts,rc=status)
                      call LVT_verify(status,&
                           'Error in ESMF_TimeIntervalSet diagnoseEntropy')
                      tindex = (nint((currTime - starttime)/timestep )+1)-1
                      
                      stats%mentropy%value_model_ts(t,k,tindex) = model%value(t,k)
                   endif
                endif
             endif
          enddo
       enddo
    endif

  end subroutine diagnoseSingleMetricEntropy

!BOP
! !ROUTINE: LVT_computeMetricEntropy
! \label{LVT_computeMetricEntropy}
! 
! !INTERFACE: 
  subroutine LVT_computeMetricEntropy(pass,alarm)

    implicit none
! !ARGUMENTS:
    integer     :: pass
    logical     :: alarm
!
! !DESCRIPTION: 
! This subroutine invokes the call to compute the metric entropy 
!  values for each variable. 
!
!   The arguments are: 
!   \begin{description}
!    \item[pass]
!     current pass index number over the data points
!    \item[alarm]
!     flag indicating if the temporal output interval has been reached
!   \end{description}
! 
!EOP
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%mentropy%selectOpt.eq.1.or.&
            LVT_metrics%mentropy%timeOpt.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call computeSingleMetricEntropy(&
                  alarm,model,obs,stats,&
                  LVT_metrics%mentropy)

             model => model%next
             obs => obs%next
             stats => stats%next             

          enddo
       endif
    endif
  end subroutine LVT_computeMetricEntropy

!BOP
! !ROUTINE: computeSingleMetricEntropy
! \label{computeSingleMetricEntropy}
! 
! !INTERFACE: 
  subroutine computeSingleMetricEntropy(alarm,model,obs,stats,metric)
! !USES: 
    use LVT_informationContentMod

    implicit none
! !ARGUMENTS: 
    logical                 :: alarm
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
!  This subroutine performs the computation of the metric entropy values 
!  for each variable. The time series values are converted to a binary 
!  string based on the computed median (All values above the median are 
!  given a value of 1 and all values below the median are given a value
!  of 0). From this binary string, words are defined basec on a given 
!  word length. The metric entropy is calculcated based on 
!  Pachepsky et al. (2006)
! 
!  The subroutine arguments are: 
!  \begin{description}
!    \item[alarm]
!     flag indicating if the temporal output interval has been reached
!    \item[model]
!      object to hold model variable information    
!    \item[obs]
!      object to hold observation variable information
!    \item[stats]
!     object to hold statistics computation related information
!    \item[metric]
!     object to hold metric-specific information
!   \end{description}
!EOP
    integer :: t,k,kk,i,l,iprev
    integer :: binval
    real    :: med_val(LVT_rc%ngrid,model%selectNlevs)
    character*1 :: bitstr(LVT_rc%ngrid,model%selectNlevs,LVT_rc%nts) 
!number of words in the string : N-L+1
    character*1 :: words(LVT_rc%nts-LVT_ICwordlength+1,LVT_ICwordlength)
    integer     :: nbins(2**LVT_ICwordlength)
    real        :: pL(2**LVT_ICwordlength)
    real        :: mentropy

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                med_val(t,k) = median(stats%mentropy%value_model_ts(t,k,:),&
                     LVT_rc%nts)
             enddo
          enddo
          
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs              
                do l=1,LVT_rc%nts
                   if(stats%mentropy%value_model_ts(t,k,l).ne.LVT_rc%udef) then 
                      if(stats%mentropy%value_model_ts(t,k,l).gt.med_val(t,k)) then 
                         bitstr(t,k,l) = '1'
                      else
                         bitstr(t,k,l) = '0'
                      endif
                   endif
                enddo
             enddo
          enddo

          do k=1,model%selectNlevs   
             do t=1,LVT_rc%ngrid
                nbins = 0 
                do kk=1,LVT_rc%nts-LVT_ICwordlength+1
                   do i=1,LVT_ICwordLength    
                      if(kk.gt.1) then 
                         l = (i+iprev)
                         if(i.eq.3) iprev = iprev+1
                      else
                         l= i
                         iprev=1
                      endif
                      words(kk,i) = bitstr(t,k,l)
                   enddo
                   call LVT_findWordIndex(words(kk,:),binval)
                   if(binval.gt.0) &
                        nbins(binval) = nbins(binval)+1
                enddo
                do l=1,2**LVT_ICwordlength
                   pL(l) = (float(nbins(l)))/float(sum(nbins))
                enddo
!metric entropy 
                mentropy = 0 
                do l=1,2**LVT_ICwordlength
                   if(pL(l).ne.0) then 
                      mentropy = mentropy -pL(l)*(log(pL(l))/log(2.0))
                   end if
                enddo

                stats%mentropy%value_final(t,k) = mentropy/float(LVT_ICwordLength)
             enddo
          enddo
          
          do k=1,model%selectNlevs
             call LVT_computeCI(stats%mentropy%value_final(:,k),&
                  LVT_rc%ngrid, LVT_rc%pval_CI,&
                  stats%mentropy%value_ci(k))
          enddo
       endif
    endif
                   
  end subroutine computeSingleMetricEntropy

  real function median(x,n)

    use LVT_SortingMod, only : LVT_sort

    IMPLICIT  NONE

    REAL,    DIMENSION(1:), INTENT(IN) :: X
    INTEGER, INTENT(IN)                :: N
    REAL,    DIMENSION(1:N)            :: Temp
    INTEGER                            :: i
    
    DO i = 1, N                       ! make a copy
       Temp(i) = X(i)
    END DO
    CALL  LVT_Sort(Temp, N)               ! sort the copy
    IF (MOD(N,2) == 0) THEN           ! compute the median
       Median = (Temp(N/2) + Temp(N/2+1)) / 2.0
    ELSE
       Median = Temp(N/2+1)
    END IF
  end function median

!BOP
! !ROUTINE: LVT_writeMetric_MetricEntropy
! \label{LVT_writeMetric_MetricEntropy}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_MetricEntropy(pass,final,vlevels,stats,obs)
! !USES:
    use LVT_statsMod, only : LVT_writeSummaryStats2
    use LVT_pluginIndices
! !ARGUMENTS: 
    integer                 :: pass
    integer                 :: final
    integer                 :: vlevels
    type(LVT_statsEntry)    :: stats
    type(LVT_metaDataEntry) :: obs
!
! !DESCRIPTION:
!   This subroutine writes the computed metric entropy values to an 
!   external file
!
!  The subroutine arguments are: 
!  \begin{description}
!    \item[pass]
!     current pass index number over the data points
!    \item[final]
!     integer flag indicating if the end of the analysis period is reached
!    \item[vlevels]
!     number of vertical levels in the current variable
!    \item[obs]
!      object to hold observation variable information
!    \item[stats]
!     object to hold statistics computation related information
!   \end{description}
!EOP
    integer                 :: count_mentropy_final(LVT_rc%ngrid,1)
    integer                 :: k,dummy

    count_mentropy_final = 1
    if(pass.eq.LVT_metrics%mentropy%npass) then 
       if(stats%selectOpt.eq.1) then 
          if(LVT_metrics%mentropy%selectOpt.eq.1) then 
             do k=1,vlevels
                
                call LVT_writevar_gridded(LVT_metrics%mentropy%ftn_total, &
                     stats%mentropy%value_final(:,k),&
                     stats%vid_total(LVT_MENTROPYid,1))
                call LVT_writevar_gridded(LVT_metrics%mentropy%ftn_total, &
                     real(count_mentropy_final(:,1)),&
                     stats%vid_count_total(LVT_MENTROPYid,1))
                
                call LVT_writeSummaryStats2(&
                     LVT_metrics%mentropy%ftn_summ,&
                     LVT_metrics%mentropy%short_name,&
                     LVT_rc%ngrid,&
                     stats%mentropy%value_final(:,k), &
                     count_mentropy_final,&
                     stats%standard_name,&
                     stats%mentropy%value_ci(k))
                
             enddo
          endif
       endif
    endif

  end subroutine LVT_writeMetric_MetricEntropy

!BOP
! 
! !ROUTINE: LVT_resetMetric_MetricEntropy
! \label{LVT_resetMetric_MetricEntropy}
!
! !INTERFACE: 
  subroutine LVT_resetMetric_MetricEntropy
!
! !DESCRIPTION: 
!  This routine resets the relevant variables between each temporal averaging
!  interval. 
! 
!EOP
  end subroutine LVT_resetMetric_MetricEntropy

!BOP
! 
! !ROUTINE: LVT_writerestart_MetricEntropy
! 
! !INTERFACE:
  subroutine LVT_writerestart_MetricEntropy(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for MetricEntropy metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%mentropy%selectOpt.eq.1) then 
       
       print*, 'The writerestart method is not implemented for MetricEntropy'
       stop
    end if
    
  end subroutine LVT_writerestart_MetricEntropy

!BOP
! 
! !ROUTINE: LVT_readrestart_MetricEntropy
! 
! !INTERFACE:
  subroutine LVT_readrestart_MetricEntropy(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for MetricEntropy metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%mentropy%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for MetricEntropy'
       stop
    end if
    
  end subroutine LVT_readrestart_MetricEntropy

end module LVT_MetricEntropyMod
