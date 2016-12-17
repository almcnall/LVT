!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
module LVT_EffectiveComplexityMod
!
!BOP
! !MODULE: LVT_EffectiveComplexityMod
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
  public :: LVT_initEffectiveComplexity
  public :: LVT_diagnoseEffectiveComplexity
  public :: LVT_computeEffectiveComplexity 
  public :: LVT_writeMetric_EffectiveComplexity
  public :: LVT_resetMetric_EffectiveComplexity
  public :: LVT_writerestart_EffectiveComplexity
  public :: LVT_readrestart_EffectiveComplexity

contains
  
  subroutine LVT_initEffectiveComplexity(model,obs,stats,metric)
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
    
    if(metric%selectOpt.eq.1) then 
       allocate(stats%ecomplexity%value_model_ts(LVT_rc%ngrid,model%selectNlevs,&
            LVT_rc%nts))
       allocate(stats%ecomplexity%value_final(LVT_rc%ngrid,model%selectNlevs))
       allocate(stats%ecomplexity%value_ci(model%selectNlevs))

       stats%ecomplexity%value_model_ts = LVT_rc%udef
       stats%ecomplexity%value_final = LVT_rc%udef
       stats%ecomplexity%value_ci = LVT_rc%udef
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


  end subroutine LVT_initEffectiveComplexity

  subroutine LVT_diagnoseEffectiveComplexity(pass)

    implicit none
    integer                 :: pass
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%ecomplexity%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call diagnoseSingleEffectiveComplexity(&
                  obs, model, stats, &
                  LVT_metrics%ecomplexity)
             
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif

  end subroutine LVT_diagnoseEffectiveComplexity

  subroutine diagnoseSingleEffectiveComplexity(obs,model,stats,metric)
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
                      
                      stats%ecomplexity%value_model_ts(t,k,tindex) = model%value(t,k)
                   endif
                endif
             endif
          enddo
       enddo
    endif

  end subroutine diagnoseSingleEffectiveComplexity

  subroutine LVT_computeEffectiveComplexity(pass,alarm)
! !ARGUMENTS:
    implicit none

    integer     :: pass
    logical     :: alarm
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats


    if(pass.eq.1) then 
       if(LVT_metrics%ecomplexity%selectOpt.eq.1.or.&
            LVT_metrics%ecomplexity%timeOpt.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call computeSingleEffectiveComplexity(&
                  alarm,model, obs, stats, &
                  LVT_metrics%ecomplexity)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_computeEffectiveComplexity

!BOP
! !ROUTINE: computeSingleEffectiveComplexity
! \label{computeSingleEffectiveComplexity}
! 
! !INTERFACE: 
  subroutine computeSingleEffectiveComplexity(alarm,model,obs,stats,metric)
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
    real        :: pLcij(2**LVT_ICwordlength,2**LVT_ICwordlength)    
    integer     :: nbins(2**LVT_ICwordlength)
    real        :: pL(2**LVT_ICwordlength)
    integer     :: binval
    real        :: ecomplexity

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                med_val(t,k) = median(stats%ecomplexity%value_model_ts(t,k,:),&
                     LVT_rc%nts)
             enddo
          enddo
          
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs              
                do l=1,LVT_rc%nts
                   if(stats%ecomplexity%value_model_ts(t,k,l).ne.LVT_rc%udef) then 
                      if(stats%ecomplexity%value_model_ts(t,k,l).gt.med_val(t,k)) then 
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
!conditional probability from the ith to the jth word, given ith word
                do i=1,(2**LVT_ICwordlength)
                   do j=1,(2**LVT_ICwordlength)
                      if(sum(nbinsij(i,:)).ne.0) then 
                         pLcij(i,j) = (float(nbinsij(i,j)))/float(sum(nbinsij(i,:)))
                      else
                         pLcij(i,j) = 0.0
                      endif
                   enddo
                enddo

                ecomplexity = 0 
                do i=1,2**LVT_ICwordlength
                   do j=1,2**LVT_ICwordlength
                      if(pLcij(i,j).ne.0.and.pL(i).ne.0) then 
                         if(log(pL(i)).ne.0) then 
                            ecomplexity = ecomplexity +pLij(i,j)*(log(pLcij(i,j))/log(pL(i)))
                         endif
                      end if
                   enddo
                enddo

                stats%ecomplexity%value_final(t,k) = ecomplexity
             enddo
          enddo
          
          do k=1,model%selectNlevs
             call LVT_computeCI(stats%ecomplexity%value_final(:,k),&
                  LVT_rc%ngrid, LVT_rc%pval_CI,&
                  stats%ecomplexity%value_ci(k))
          enddo
       endif
    endif
                   
  end subroutine computeSingleEffectiveComplexity

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


  subroutine LVT_writeMetric_EffectiveComplexity(pass,final,vlevels,stats,obs)
! !USES:
    use LVT_statsMod, only : LVT_writeSummaryStats2
    use LVT_pluginIndices
! !ARGUMENTS: 
    integer                 :: pass
    integer                 :: final
    integer                 :: vlevels
    type(LVT_statsEntry)    :: stats
    type(LVT_metaDataEntry) :: obs
    integer                 :: count_ecomplexity_final(LVT_rc%ngrid,1)
    integer                 :: k,dummy

    count_ecomplexity_final = 1
    if(pass.eq.LVT_metrics%ecomplexity%npass) then 
       if(stats%selectOpt.eq.1) then 
          if(LVT_metrics%ecomplexity%selectOpt.eq.1) then 
             do k=1,vlevels
                call LVT_writevar_gridded(LVT_metrics%ecomplexity%ftn_total, &
                     stats%ecomplexity%value_final(:,k),&
                     stats%vid_total(LVT_ECOMPLEXITYid,1))
                call LVT_writevar_gridded(LVT_metrics%ecomplexity%ftn_total, &
                     real(count_ecomplexity_final(:,1)),&
                     stats%vid_count_total(LVT_ECOMPLEXITYid,1))
                
                call LVT_writeSummaryStats2(&
                     LVT_metrics%ecomplexity%ftn_summ,&
                     LVT_metrics%ecomplexity%short_name,&
                     LVT_rc%ngrid,&
                     stats%ecomplexity%value_final(:,k), &
                     count_ecomplexity_final,&
                     stats%standard_name,&
                     stats%ecomplexity%value_ci(k))
                
             enddo
          endif
       endif
    endif

  end subroutine LVT_writeMetric_EffectiveComplexity

  subroutine LVT_resetMetric_EffectiveComplexity

  end subroutine LVT_resetMetric_EffectiveComplexity


!BOP
! 
! !ROUTINE: LVT_writerestart_EffectiveComplexity
! 
! !INTERFACE:
  subroutine LVT_writerestart_EffectiveComplexity(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for EffectiveComplexity metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ecomplexity%selectOpt.eq.1) then 
       
       print*, 'The writerestart method is not implemented for EffectiveComplexity'
       stop
    end if
    
  end subroutine LVT_writerestart_EffectiveComplexity

!BOP
! 
! !ROUTINE: LVT_readrestart_EffectiveComplexity
! 
! !INTERFACE:
  subroutine LVT_readrestart_EffectiveComplexity(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for EffectiveComplexity metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ecomplexity%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for EffectiveComplexity'
       stop
    end if
    
  end subroutine LVT_readrestart_EffectiveComplexity

end module LVT_EffectiveComplexityMod
