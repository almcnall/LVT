!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
module LVT_FluctuationComplexityMod
!
!BOP
! !MODULE: LVT_FluctuationComplexityMod
! 
!  !DESCRIPTION: 
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
  public :: LVT_initFluctuationComplexity
  public :: LVT_diagnoseFluctuationComplexity
  public :: LVT_computeFluctuationComplexity 
  public :: LVT_writeMetric_FluctuationComplexity
  public :: LVT_resetMetric_FluctuationComplexity
  public :: LVT_writerestart_FluctuationComplexity
  public :: LVT_readrestart_FluctuationComplexity

contains
  
  subroutine LVT_initFluctuationComplexity(model,obs,stats,metric)
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
    
    if(metric%selectOpt.eq.1) then 
       allocate(stats%fcomplexity%value_model_ts(LVT_rc%ngrid,model%selectNlevs,&
            LVT_rc%nts))
       allocate(stats%fcomplexity%value_final(LVT_rc%ngrid,model%selectNlevs))
       allocate(stats%fcomplexity%value_ci(model%selectNlevs))

       stats%fcomplexity%value_model_ts = LVT_rc%udef
       stats%fcomplexity%value_final = LVT_rc%udef
       stats%fcomplexity%value_ci = LVT_rc%udef
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


  end subroutine LVT_initFluctuationComplexity

  subroutine LVT_diagnoseFluctuationComplexity(pass)

    implicit none
    integer                 :: pass
    integer       :: i
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%fcomplexity%selectOpt.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleFluctuationComplexity(&
                  obs, model, stats, &
                  LVT_metrics%fcomplexity)

             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
       endif
    endif

  end subroutine LVT_diagnoseFluctuationComplexity

  subroutine diagnoseSingleFluctuationComplexity(obs,model,stats,metric)
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
! 
!  while diagnosing, simply gather the soil moisture values into an array
! 
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
                      
                      stats%fcomplexity%value_model_ts(t,k,tindex) = model%value(t,k)
                   endif
                endif
             endif
          enddo
       enddo
    endif

  end subroutine diagnoseSingleFluctuationComplexity

  subroutine LVT_computeFluctuationComplexity(pass,alarm)
! !ARGUMENTS:
    implicit none

    integer     :: pass
    logical     :: alarm

    integer     :: i
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%fcomplexity%selectOpt.eq.1.or.&
            LVT_metrics%fcomplexity%timeOpt.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)
          
          do while(associated(model))

             call computeSingleFluctuationComplexity(&
                  alarm,model,obs,stats,&
                  LVT_metrics%fcomplexity)

             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
       endif
    endif
  end subroutine LVT_computeFluctuationComplexity

!BOP
! !ROUTINE: computeSingleFluctuationComplexity
! \label{computeSingleFluctuationComplexity}
! 
! !INTERFACE: 
  subroutine computeSingleFluctuationComplexity(alarm,model,obs,stats,metric)
! !USES: 
    use LVT_informationContentMod

    implicit none
! !ARGUMENTS: 
    logical                 :: alarm
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer :: t,k,kk,i,j,l,iprev
    integer :: col,row
    real    :: med_val(LVT_rc%ngrid,model%selectNlevs)
    character*1 :: bitstr(LVT_rc%ngrid,model%selectNlevs,LVT_rc%nts) 
!number of words in the string : N-L+1
    character*1 :: words(LVT_rc%nts-LVT_ICwordlength+1,LVT_ICwordlength)
    integer     :: nbinsij(2**LVT_ICwordlength,2**LVT_ICwordlength)
    real        :: pLij(2**LVT_ICwordlength,2**LVT_ICwordlength)
    integer     :: nbins(2**LVT_ICwordlength)
    real        :: pL(2**LVT_ICwordlength)
    integer     :: binval
    real        :: fcomplexity

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                med_val(t,k) = median(stats%fcomplexity%value_model_ts(t,k,:),&
                     LVT_rc%nts)
             enddo
          enddo
          
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs              
                do l=1,LVT_rc%nts
                   if(stats%fcomplexity%value_model_ts(t,k,l).ne.LVT_rc%udef) then 
                      if(stats%fcomplexity%value_model_ts(t,k,l).gt.med_val(t,k)) then 
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
                nbinsij = 0 
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
                enddo
                do kk=1,LVT_rc%nts-LVT_ICwordlength+1
                   call LVT_findWordIndex(words(kk,:),binval)
                   if(binval.gt.0) &
                        nbins(binval) = nbins(binval)+1
                enddo
                do l=1,2**LVT_ICwordlength
                   pL(l) = (float(nbins(l)))/float(sum(nbins))
                enddo
                do kk=1,LVT_rc%nts-LVT_ICwordlength+1
                   if(kk.gt.1) then
                      call LVT_findWordTransitionIndex(&
                           words(kk,:),words(kk-1,:), col,row)
                      nbinsij(col,row) = nbinsij(col,row) + 1
                   endif
                enddo
!probability of transition from the ith to the jth word
                do i=1,(2**LVT_ICwordlength)
                   do j=1,(2**LVT_ICwordlength)
                      if(sum(nbinsij).ne.0) then 
                         pLij(i,j) = (float(nbinsij(i,j)))/float(sum(nbinsij))
                      else
                         pLij(i,j) = 0.0
                      endif
                   enddo
                enddo

                fcomplexity = 0 
                do i=1,2**LVT_ICwordlength
                   do j=1,2**LVT_ICwordlength
                      if(pL(j).ne.0.and.pL(i).ne.0) then
                         if(log(pL(j)).ne.0) then 
                            fcomplexity = fcomplexity +pLij(i,j)*(log(pL(i))/log(pL(j)))**2
                         endif
                      end if
                   enddo
                enddo

                stats%fcomplexity%value_final(t,k) = fcomplexity
             enddo
          enddo
          
          do k=1,model%selectNlevs
             call LVT_computeCI(stats%fcomplexity%value_final(:,k),&
                  LVT_rc%ngrid, LVT_rc%pval_CI,&
                  stats%fcomplexity%value_ci(k))
          enddo
       endif
    endif
                   
  end subroutine computeSingleFluctuationComplexity

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


  subroutine LVT_writeMetric_FluctuationComplexity(pass,final,vlevels,stats,obs)
! !USES:
    use LVT_statsMod, only : LVT_writeSummaryStats2
    use LVT_pluginIndices
! !ARGUMENTS: 
    integer                 :: pass
    integer                 :: final
    integer                 :: vlevels
    type(LVT_statsEntry)    :: stats
    type(LVT_metaDataEntry) :: obs
    integer                 :: count_fcomplexity_final(LVT_rc%ngrid,1)
    integer                 :: k,dummy

    count_fcomplexity_final = 1
    if(pass.eq.LVT_metrics%fcomplexity%npass) then 
       if(stats%selectOpt.eq.1) then 
          if(LVT_metrics%fcomplexity%selectOpt.eq.1) then 
             do k=1,vlevels
                call LVT_writevar_gridded(LVT_metrics%fcomplexity%ftn_total, &
                     stats%fcomplexity%value_final(:,k),&
                     stats%vid_total(LVT_FCOMPLEXITYid,1))
                call LVT_writevar_gridded(LVT_metrics%fcomplexity%ftn_total, &
                     real(count_fcomplexity_final(:,1)),&
                     stats%vid_count_total(LVT_FCOMPLEXITYid,1))
                
                call LVT_writeSummaryStats2(&
                     LVT_metrics%fcomplexity%ftn_summ,&
                     LVT_metrics%fcomplexity%short_name,&
                     LVT_rc%ngrid,&
                     stats%fcomplexity%value_final(:,k), &
                     count_fcomplexity_final,&
                     stats%standard_name,&
                     stats%fcomplexity%value_ci(k))
                
             enddo
          endif
       endif
    endif

  end subroutine LVT_writeMetric_FluctuationComplexity

  subroutine LVT_resetMetric_FluctuationComplexity

  end subroutine LVT_resetMetric_FluctuationComplexity


!BOP
! 
! !ROUTINE: LVT_writerestart_FluctuationComplexity
! 
! !INTERFACE:
  subroutine LVT_writerestart_FluctuationComplexity(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for FluctuationComplexity metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%fcomplexity%selectOpt.eq.1) then 
       
       print*, 'The writerestart method is not implemented for FluctuationComplexity'
       stop
    end if
    
  end subroutine LVT_writerestart_FluctuationComplexity

!BOP
! 
! !ROUTINE: LVT_readrestart_FluctuationComplexity
! 
! !INTERFACE:
  subroutine LVT_readrestart_FluctuationComplexity(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for FluctuationComplexity metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%fcomplexity%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for FluctuationComplexity'
       stop
    end if
    
  end subroutine LVT_readrestart_FluctuationComplexity

end module LVT_FluctuationComplexityMod
