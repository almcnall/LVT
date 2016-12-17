!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
module LVT_DataStreamsMod
!BOP
! 
! !MODULE: LVT_DataStreamMod
! \label(LVT_DataStreamMod)
!
! !INTERFACE:
! 
! !USES:   
  use LVT_histDataMod
  use LVT_coreMod
  use LVT_logMod
  use LVT_LISoutputHandlerMod

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  The code in this file contains the basic datastructures and 
!  control routines for handling the operations associated 
!  with the datastreams. The invocation to read, perform 
!  temporal averaging and resetting of the datastreams are 
!  performed from this module. The calculations of derived
!  variables are also performed in this module.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_DataStreamsInit
  public :: LVT_readDataStreams
  public :: LVT_tavgDataStreams
  public :: LVT_resetDataStreams
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP

contains

!BOP
! 
! !ROUTINE: LVT_DataStreamsInit
! \label{LVT_DataStreamsInit}
!
! !INTERFACE:   
  subroutine LVT_DataStreamsInit
! 
! !USES:   
    use LVT_datastream_pluginMod, only : LVT_datastream_plugin
    implicit none
!
! !DESCRIPTION:
! 
!  This subroutine invokes the call to initialize each datastream.
!
!   The routines invoked are: 
!   \begin{description}
!    \item[LVT\_datastream\_plugin] (\ref{LVT_datastream_plugin}) \newline
!      routine to register all the supported datastream plugin implementations
!    \item[LVT\_LISoutputInit] (\ref{LVT_LISoutputInit}) \newline
!      routine to initialize the handling of LIS output data (if LIS output
!      data is one of the datastreams)
!   \end{description} 
!EOP
    integer      :: kk
    type(LVT_metadataEntry), pointer :: ds1, ds2

    call LVT_datastream_plugin
    
    call observationsetup(trim(LVT_rc%obssource(1))//char(0),1)
    call observationsetup(trim(LVT_rc%obssource(2))//char(0),2)

    call LVT_LISoutputInit


!checking for duplicate entries in a given datastream
    LVT_rc%ds1_dup = .false. 
    ds1 => LVT_histData%head_ds1_list       
    do while(associated(ds1))
       ds2 => ds1%next
       do while(associated(ds2))
          if(ds2%index.ne.ds1%index.and.&
               ds1%short_name.eq.ds2%short_name) then 
             LVT_rc%ds1_dup = .true. 
          endif
          ds2 => ds2%next
       enddo
       ds1 => ds1%next
    enddo

    LVT_rc%ds2_dup = .false. 
    ds1 => LVT_histData%head_ds2_list       
    do while(associated(ds1))
       ds2 => ds1%next
       do while(associated(ds2))
          if(ds2%index.ne.ds1%index.and.&
               ds1%short_name.eq.ds2%short_name) then 
             LVT_rc%ds2_dup = .true. 
          endif
          ds2 => ds2%next
       enddo
       ds1 => ds1%next
    enddo
  end subroutine LVT_DataStreamsInit

!BOP
! 
! !ROUTINE: LVT_readDataStreams
! \label{LVT_readDataStreams}
!
! !INTERFACE: 
  subroutine LVT_readDataStreams
! 
! !USES:   
    implicit none
!
!
! !DESCRIPTION: 
!  This subroutine invokes the routines that read the datastreams
! 
!EOP

    call readObservationSource(trim(LVT_rc%obssource(1))//char(0),1)      
    call readObservationSource(trim(LVT_rc%obssource(2))//char(0),2)      
    
  end subroutine LVT_readDataStreams

!BOP
! 
! !ROUTINE: LVT_tavgDataStreams
! \label{LVT_tavgDataStreams}
!
! !INTERFACE: 
  subroutine LVT_tavgDataStreams
! 
! !USES:   
    use LVT_statsDataMod

    implicit none
!
!
! !DESCRIPTION: 
!   This routine invokes the calls to compute temporal averages of 
!   desired set of variables, based on the specified 
!   temporal averaging frequency. 
!  
!   The routines invoked are: 
!   \begin{description}
!    \item[tavgSingleDataStream](\ref{tavgSingleDataStream})
!     computes the temporal average for a single variable
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer      :: kk
    type(LVT_metadataEntry), pointer :: dataEntry
    type(LVT_metadataEntry), pointer :: ds1, ds2

    if(LVT_rc%computeFlag) then            

!data stream 1
       do kk=1,2
          if(kk.eq.1) then 
             dataEntry => LVT_histData%head_ds1_list
          elseif(kk.eq.2) then 
             dataEntry => LVT_histData%head_ds2_list
          endif
       
          do while(associated(dataEntry))
             call tavgSingleDataStream(dataEntry)
             dataEntry => dataEntry%next
          enddo

!copy duplicate entries
          if(LVT_rc%ds1_dup) then 
             ds1 => LVT_histData%head_ds1_list       
             do while(associated(ds1))
                ds2 => ds1%next
                do while(associated(ds2))
                   if(ds2%index.ne.ds1%index.and.&
                        ds1%short_name.eq.ds2%short_name) then 
                      ds2%value = ds1%value
                      ds2%count = ds1%count
                   endif
                   ds2 => ds2%next
                enddo
                ds1 => ds1%next
             enddo
             
          endif

          if(LVT_rc%ds2_dup) then 
             ds1 => LVT_histData%head_ds2_list       
             do while(associated(ds1))
                ds2 => ds1%next
                do while(associated(ds2))
                   if(ds2%index.ne.ds1%index.and.&
                        ds1%short_name.eq.ds2%short_name) then 
                      ds2%value = ds1%value
                      ds2%count = ds1%count
                   endif
                   ds2 => ds2%next
                enddo
                ds1 => ds1%next
             enddo
             
          endif

          if(LVT_rc%var_based_strat .gt. 0) then
             call tavgSingleDataStream(LVT_histData%strat_varEntry)
             LVT_stats%strat_var(:,:) = LVT_histData%strat_varEntry%value(:,:)
          endif
       enddo
    endif
  end subroutine LVT_tavgDataStreams

!BOP
! 
! !ROUTINE: tavgSingleDataStream
! \label{tavgSingleDataStream}
!
! !INTERFACE:
  subroutine tavgSingleDataStream( dataEntry)
! 
! !USES: 

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine temporally averages the accumulated data in a
!   given datastream
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    type(LVT_metadataEntry) :: dataEntry
!EOP
    integer :: k,t,c,r, gid

    if(dataEntry%selectNlevs.ge.1) then 
       if(LVT_rc%computeEnsMetrics.eq.1) then 
          do t=1,LVT_LIS_rc(1)%ntiles
             do k=1,dataEntry%vlevels
                c = LVT_LIS_domain(1)%tile(t)%col
                r = LVT_LIS_domain(1)%tile(t)%row
                if(LVT_LIS_domain(1)%gindex(c,r).ne.-1) then 
                   gid = LVT_LIS_domain(1)%gindex(c,r)
                   if(dataEntry%count(t,k).ne.0) then 
                      dataEntry%value(t,k) = &
                           dataEntry%value(t,k)/dataEntry%count(t,k)

                   endif
                endif
             enddo
          enddo         
       else
          do r=1,LVT_rc%lnr
             do c=1,LVT_rc%lnc
                do k=1,dataEntry%vlevels
                   if(LVT_domain%gindex(c,r).ne.-1) then 
                      gid = LVT_domain%gindex(c,r)
                      if(dataEntry%count(gid,k).ne.0) then 
                         dataEntry%value(gid,k) = &
                              dataEntry%value(gid,k)/dataEntry%count(gid,k)

                      endif
                   endif
                enddo
             enddo
          enddo
       endif
    endif

  end subroutine tavgSingleDataStream


!BOP
! 
! !ROUTINE: LVT_resetDataStreams
! \label{LVT_resetDataStreams}
!
! !INTERFACE: 
  subroutine LVT_resetDataStreams
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine reinitializes the data structures that hold the observational
!   data
! 
!   The routines invoked are: 
!   \begin{description}
!    \item[resetSingleDataStream2](\ref{resetSingleDataStream2})
!     resets the datastructures for a single variable
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    type(LVT_metadataEntry), pointer :: ds1
    type(LVT_metadataEntry), pointer :: ds2

    if(LVT_rc%computeFlag) then            
!data stream 1
       ds1 => LVT_histData%head_ds1_list
       
       do while(associated(ds1))
          call resetSingleDataStream(ds1)
          ds1 => ds1%next
       enddo
       
!data stream 2
       ds2 => LVT_histData%head_ds2_list
       
       do while(associated(ds2))
          call resetSingleDataStream(ds2)
          ds2 => ds2%next
       enddo
       
!need special handler for LIS output
       if(LVT_rc%lis_output_obs) then 
          if(LVT_rc%obssource(1).eq."LIS output") then 
             call LVT_resetLISoutputContainers(1)
          endif
          if(LVT_rc%obssource(2).eq."LIS output") then 
             call LVT_resetLISoutputContainers(2)
          endif
       endif
    endif
  end subroutine LVT_resetDataStreams


!BOP
! 
! !ROUTINE: resetSingleDataStream
! \label{resetSingleDataStream}
!
! !INTERFACE: 
  subroutine resetSingleDataStream(dataEntry)
! 
! !USES:   
    implicit none 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine resets the data structures that hold the observational 
!  data and the temporal averaging counters
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS: 
    type(LVT_metadataEntry) :: dataEntry
! 
!EOP

    integer                 :: k 

    if(dataEntry%selectNlevs.ge.1) then 
       do k=1,dataEntry%vlevels
          dataEntry%value(:,k) = 0 
          dataEntry%count(:,k) = 0 
          if(dataEntry%stdev_flag) then 
             dataEntry%count_stdev(:,k)= 0 
             dataEntry%stdev(:,k) = 0
          endif
       enddo      
    endif
  end subroutine resetSingleDataStream
end module LVT_DataStreamsMod
