!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_ensXcorrMod
! \label(LVT_ensXcorrMod)
!
! !INTERFACE:
module LVT_ensXcorrMod
! 
! !USES:   
  use ESMF
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
!  !DESCRIPTION: 
!   This module handles the computations required to compute 
!   ensemble mean values
!   of desired variables from the LIS output. Note that the LIS
!   output is expected to be written in tile space format. 
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
  public :: LVT_initensXcorr
  public :: LVT_diagnoseensXcorr
  public :: LVT_computeensXcorr
  public :: LVT_writeMetric_ensXcorr
  public :: LVT_resetMetric_ensXcorr
  public :: LVT_writerestart_ensXcorr
  public :: LVT_readrestart_ensXcorr

!EOP
  
  private

contains
  subroutine LVT_initensXcorr(model, obs, stats,metric)    
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
! 
!EOP
    character*100   :: filename
    integer         :: t,k,iter_f
    integer         :: ftn
    integer         :: rc
    character*100   :: pname
    real            :: dummy(LVT_LIS_rc(1)%ntiles)

    if(metric%selectOpt.eq.1) then 
!read the parameter ensemble
       call ESMF_ConfigGetAttribute(LVT_config,filename,&
            label="LIS OptUE restart file:",rc=rc)
       call LVT_verify(rc,'LIS OptUE restart file: option not specified in the config file')

       call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%nparam,&
            label="LIS OptUE number of model parameters:",rc=rc)
       call LVT_verify(rc,"LIS OptUE number of model parameters: option not defined in the config file")
       allocate(stats%ensxcorr%value_param(LVT_LIS_rc(1)%ntiles,LVT_rc%nparam))

       ftn = LVT_getNextUnitNumber()
       write(LVT_logunit,*) 'Reading ',trim(filename)
       open(ftn,file=trim(filename),form='unformatted')
       read(ftn) iter_f
       do t=1,LVT_rc%nparam+1
          if(t==LVT_rc%nparam+1) then 
             read(ftn) dummy
          else
             read(ftn) pname
             read(ftn) dummy
!             do k=1,LVT_LIS_rc(1)%ntiles
!                print*, t, LVT_rc%nparam, k,trim(pname), dummy(k)
!             enddo
             stats%ensxcorr%value_param(:,t) = dummy(:)
          endif
       enddo

       call LVT_releaseUnitNumber(ftn)
       write(LVT_logunit,*) 'Finished reading ',trim(filename)

       allocate(stats%ensxcorr%value_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%nparam))
       stats%ensxcorr%value_total = 0.0
       allocate(stats%ensxcorr%count_value_total(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%nparam))
       stats%ensxcorr%count_value_total = 0
       allocate(stats%ensxcorr%value_ci(model%selectNlevs,LVT_rc%nparam))
       stats%ensxcorr%value_ci = LVT_rc%udef       
    endif
    if(metric%timeOpt.eq.1) then 
       allocate(stats%ensxcorr%value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%nparam))
       stats%ensxcorr%value_ts = 0.0
       allocate(stats%ensxcorr%count_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%nparam))
       stats%ensxcorr%count_value_ts = 0  

       allocate(stats%ensxcorr%tavg_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%nparam))
       stats%ensxcorr%tavg_value_ts = 0.0
       allocate(stats%ensxcorr%tavg_count_value_ts(LVT_rc%ngrid, model%selectNlevs, &
            LVT_rc%nparam))
       stats%ensxcorr%tavg_count_value_ts = 0       
    endif

!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1    
    metric%obsData = .true. 

  end subroutine LVT_initensXcorr

!BOP
! 
! !ROUTINE: LVT_diagnoseensXcorr
! \label{LVT_diagnoseensXcorr}
!
! !INTERFACE: 
  subroutine LVT_diagnoseensXcorr(pass)
! 
! !USES:     

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the computations for 
!   calculating the mean of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelensXcorr](\ref{diagnoseSingleModelensXcorr})
!     updates the mean computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
    integer       :: pass
!EOP

    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats


    if(pass.eq.1) then 
       if(LVT_metrics%ensxcorr%selectOpt.eq.1.or.&
            LVT_metrics%ensxcorr%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleensXcorr(model, obs, stats, &
                  LVT_metrics%ensxcorr)
             
             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
       endif
    endif
  end subroutine LVT_diagnoseensXcorr

!BOP
! 
! !ROUTINE: diagnoseSingleensXcorr
! \label{diagnoseSingleensXcorr}
!
! !INTERFACE: 
  subroutine diagnoseSingleensXcorr(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the mean computation of the 
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
!EOP
    integer    :: t,k,kk,tind,gid
    real       :: sxy(LVT_rc%ngrid,model%selectNlevs,LVT_rc%nparam)
    real       :: sx(LVT_rc%ngrid,model%selectNlevs,LVT_rc%nparam)
    real       :: sy(LVT_rc%ngrid,model%selectNlevs,LVT_rc%nparam)
    real       :: sxx(LVT_rc%ngrid,model%selectNlevs,LVT_rc%nparam)
    real       :: syy(LVT_rc%ngrid,model%selectNlevs,LVT_rc%nparam)
    integer    :: ncount(LVT_rc%ngrid,model%selectNlevs,LVT_rc%nparam)
    real       :: xcorr(LVT_rc%ngrid,model%selectNlevs,LVT_rc%nparam)
    real       :: numer, denom, denom_sq

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then        
       sxy = 0
       sx = 0 
       sy = 0 
       sxx = 0 
       syy = 0 
       ncount = 0 

       do t=1,LVT_LIS_rc(1)%ntiles
          gid = LVT_domain%gindex(LVT_LIS_domain(1)%tile(t)%col,&
               LVT_LIS_domain(1)%tile(t)%row)
          do k=1,model%selectNlevs
             if(model%count(t,k).gt.0) then
                if(metric%selectOpt.eq.1) then 
                   if(model%value(t,k).ne.LVT_rc%udef) then 
                      do kk=1,LVT_rc%nparam
                         sxy(gid,k,kk) = sxy(gid,k,kk) +model%value(t,k)*&
                              stats%ensxcorr%value_param(t,kk)
                         sx(gid,k,kk) = sx(gid,k,kk) + model%value(t,k)
                         sy(gid,k,kk) = sy(gid,k,kk) + &
                              stats%ensxcorr%value_param(t,kk)
                         sxx(gid,k,kk) = sxx(gid,k,kk) + model%value(t,k)*&
                              model%value(t,k)
                         syy(gid,k,kk) = syy(gid,k,kk) + &
                              stats%ensxcorr%value_param(t,kk)*&
                              stats%ensxcorr%value_param(t,kk)
                         ncount(gid,k,kk) = ncount(gid,k,kk) + 1
                      enddo
                   endif
                endif
             endif
          enddo
       enddo
       do gid=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             do kk=1,LVT_rc%nparam
                numer = (float(ncount(gid,k,kk))* &
                     sxy(gid,k,kk) - sx(gid,k,kk)*sy(gid,k,kk))
                denom_sq =  (float(ncount(gid,k,kk))*sxx(gid,k,kk)-&
                     sx(gid,k,kk)**2)* &
                     (float(ncount(gid,k,kk))*syy(gid,k,kk)- sy(gid,k,kk)**2)
                if(denom_sq.gt.0) then 
                   denom = sqrt(denom_sq)
                else
                   denom = 0.0
                endif
                if(denom.ne.0) then
                   xcorr(gid,k,kk) = numer/denom
                else
                   xcorr(gid,k,kk) = LVT_rc%udef
                endif
             enddo
          enddo
       enddo
       do gid=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             do kk=1,LVT_rc%nparam
                if(xcorr(gid,k,kk).ne.LVT_rc%udef) then 
                   stats%ensxcorr%value_total(gid,k,kk) = &
                        stats%ensxcorr%value_total(gid,k,kk)+&
                        xcorr(gid,k,kk)
                   stats%ensxcorr%count_value_total(gid,k,kk) = & 
                        stats%ensxcorr%count_value_total(gid,k,kk) +1
                endif
             enddo
          enddo
       enddo
       if(metric%timeOpt.eq.1) then 
          do gid=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do kk=1,LVT_rc%nparam
                   if(xcorr(gid,k,kk).ne.LVT_rc%udef) then 
                      stats%ensxcorr%value_ts(gid,k,kk) = &
                           stats%ensxcorr%value_ts(gid,k,kk)+&
                           xcorr(gid,k,kk)
                      stats%ensxcorr%count_value_ts(gid,k,kk) = & 
                           stats%ensxcorr%count_value_ts(gid,k,kk) +1
                   endif
                enddo
             enddo
          enddo
       endif
    endif
  end subroutine diagnoseSingleensXcorr

!BOP
! 
! !ROUTINE: LVT_computeensXcorr
! \label{LVT_computeensXcorr}
!
! !INTERFACE: 
  subroutine LVT_computeensXcorr(pass,alarm)
! 
! !USES: 

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the mean values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelensXcorr](\ref{computeSingleModelensXcorr})
!     computes the mean values for a single variable
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
    integer               :: pass
    logical               :: alarm
!EOP
    integer     :: i 
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats
    
    if(pass.eq.1) then 
       if(LVT_metrics%ensxcorr%selectOpt.eq.1.or.&
            LVT_metrics%ensxcorr%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%ensxcorr%timeOpt.eq.1.and.&
                  LVT_metrics%ensxcorr%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ensxcorr%ftn_ts_loc(i),200,advance='no') &
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
             
             call computeSingleensXcorr(alarm,model, obs, stats, &
                  LVT_metrics%ensxcorr)
             
             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
          
          if(alarm) then 
             if(LVT_metrics%ensxcorr%timeOpt.eq.1.and.&
                  LVT_metrics%ensxcorr%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%ensxcorr%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_computeensXcorr
  

!BOP
! 
! !ROUTINE: computeSingleensXcorr
! \label{computeSingleensXcorr}
!
! !INTERFACE: 
  subroutine computeSingleensXcorr(alarm,model,obs,stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the mean values
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

    integer  :: t,l,k,dummy

    if(metric%timeOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%nparam
                   if(stats%ensxcorr%count_value_ts(t,k,l).gt.0) then 
                      stats%ensxcorr%value_ts(t,k,l) = stats%ensxcorr%value_ts(t,k,l)/&
                           stats%ensxcorr%count_value_ts(t,k,l)                   

                      stats%ensxcorr%tavg_value_ts(t,k,l) = &
                           stats%ensxcorr%tavg_value_ts(t,k,l) + & 
                           stats%ensxcorr%value_ts(t,k,l)
                      stats%ensxcorr%tavg_count_value_ts(t,k,l) = &
                           stats%ensxcorr%tavg_count_value_ts(t,k,l) + 1
                   else
                      stats%ensxcorr%value_ts(t,k,l) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
          
          if(alarm) then 

             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nparam
                      if(stats%ensxcorr%tavg_count_value_ts(t,k,l).gt.0) then
                         stats%ensxcorr%tavg_value_ts(t,k,l) = &
                              stats%ensxcorr%tavg_value_ts(t,k,l) / & 
                              stats%ensxcorr%tavg_count_value_ts(t,k,l) 
                      else
                         stats%ensxcorr%tavg_value_ts(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                enddo
             enddo

             if(metric%extractTS.eq.1) then 
                call LVT_writeTSinfo(metric%ftn_ts_loc,&
                     model,&
                     LVT_rc%ngrid,&
                     stats%ensxcorr%tavg_value_ts,&
                     stats%ensxcorr%tavg_count_value_ts,dummy)
             endif
          endif
       endif
    endif
       
    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%nparam
                   if(stats%ensxcorr%count_value_total(t,k,l).gt.&
                        LVT_rc%obsCountThreshold) then 
                      stats%ensxcorr%value_total(t,k,l) = &
                           stats%ensxcorr%value_total(t,k,l)/&
                              stats%ensxcorr%count_value_total(t,k,l)           
                   else
                      stats%ensxcorr%value_total(t,k,l) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
             
          do k=1,model%selectNlevs
             do l=1, LVT_rc%nparam
                call LVT_computeCI(stats%ensxcorr%value_total(:,k,l),&
                     LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%ensxcorr%value_ci(k,l))
             enddo
          enddo
          
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid,stats%ensxcorr%value_total)
          
       endif
    endif
  end subroutine computeSingleensXcorr


  subroutine LVT_writeMetric_ensXcorr(pass,final,vlevels,stats,obs)
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

    integer                 :: k,l,tind

    if(pass.eq.LVT_metrics%ensxcorr%npass) then 
       if(final.ne.1) then
           if(stats%selectOpt.eq.1) then 
             do k=1,vlevels
                if(LVT_metrics%ensxcorr%timeOpt.eq.1) then 
                   do l=1,LVT_rc%nparam
                      call LVT_writevar_gridded(LVT_metrics%ensxcorr%ftn_ts, &
                           stats%ensxcorr%value_ts(:,k,l),&
                           stats%vid_ts(LVT_ensXcorrid,1),1)
                      call LVT_writevar_gridded(LVT_metrics%ensxcorr%ftn_ts, &
                           real(stats%ensxcorr%count_value_ts(:,k,l)),&
                           stats%vid_count_ts(LVT_ensXcorrid,1),1)
                   enddo
                endif
             enddo
          endif
       else
          if(pass.eq.1) then 
             if(stats%selectOpt.eq.1) then 
                do k=1,vlevels
                   if(LVT_metrics%ensxcorr%selectOpt.eq.1) then 
                      do l=1,LVT_rc%nparam                      
                         call LVT_writevar_gridded(&
                              LVT_metrics%ensxcorr%ftn_total, &
                              stats%ensxcorr%value_total(:,k,l),&
                              stats%vid_total(LVT_ensXcorrid,1))
                         call LVT_writevar_gridded(LVT_metrics%ensxcorr%ftn_total, &
                              real(stats%ensxcorr%count_value_total(:,k,l)),&
                              stats%vid_count_total(LVT_ensXcorrid,1))
                         
                      enddo
                      call LVT_writeSummaryStats(&
                           LVT_metrics%ensxcorr%ftn_summ,&
                           LVT_metrics%ensxcorr%short_name,&
                           LVT_rc%ngrid,&
                           stats%ensxcorr%value_total(:,k,:), &
                           stats%ensxcorr%count_value_total(:,k,:),&
                           stats%standard_name,&
                           stats%ensxcorr%value_ci(k,:))
                   endif
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_writeMetric_ensXcorr
!BOP
! 
! !ROUTINE: LVT_resetMetric_ensXcorr
! \label(LVT_resetMetric_ensXcorr)
!
! !INTERFACE:
  subroutine LVT_resetMetric_ensXcorr(alarm)
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
             if(LVT_metrics%ensxcorr%timeOpt.eq.1) then 
                do l=1,LVT_rc%nparam
                   stats%ensxcorr%value_ts(:,k,l) = 0.0
                   stats%ensxcorr%count_value_ts(:,k,l)=0 
                   if(alarm) then 
                      stats%ensxcorr%tavg_value_ts(:,k,l) = 0.0
                      stats%ensxcorr%tavg_count_value_ts(:,k,l)=0 
                   endif
                enddo
             endif
          enddo
       endif
       
       model => model%next
       stats => stats%next
    enddo

  end subroutine LVT_resetMetric_ensXcorr


!BOP
! 
! !ROUTINE: LVT_writerestart_ensXcorr
! 
! !INTERFACE:
  subroutine LVT_writerestart_ensXcorr(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for ensXcorr metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ensxcorr%selectOpt.eq.1) then 
       
       print*, 'The writerestart method is not implemented for ensXcorr'

    end if
    
  end subroutine LVT_writerestart_ensXcorr

!BOP
! 
! !ROUTINE: LVT_readrestart_ensXcorr
! 
! !INTERFACE:
  subroutine LVT_readrestart_ensXcorr(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for ensXcorr metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ensxcorr%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for ensXcorr'

    end if
    
  end subroutine LVT_readrestart_ensXcorr

end module LVT_ensXcorrMod
