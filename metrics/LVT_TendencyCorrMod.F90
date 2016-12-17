!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_TendencyCorrMod
! \label(LVT_TendencyCorrMod)
!
! !INTERFACE:
module LVT_TendencyCorrMod
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
!  This module handles the computations required to compute mean values
!  of desired variables from the LIS output
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initTendencyCorr
  public :: LVT_diagnoseTendencyCorr
  public :: LVT_computeTendencyCorr
  public :: LVT_writeMetric_TendencyCorr
  public :: LVT_resetMetric_TendencyCorr
  public :: LVT_writerestart_TendencyCorr
  public :: LVT_readrestart_TendencyCorr
!EOP
  
  private

contains
  subroutine LVT_initTendencyCorr(model, obs, stats,metric)
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
       allocate(stats%tendencycorr%sxy_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
       allocate(stats%tendencycorr%sxx_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
       allocate(stats%tendencycorr%syy_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
       allocate(stats%tendencycorr%sx_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
       allocate(stats%tendencycorr%sy_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
       allocate(stats%tendencycorr%rval_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
       allocate(stats%tendencycorr%count_value(LVT_rc%ngrid, model%selectNlevs, LVT_rc%strat_nlevels))       
       stats%tendencycorr%sxy_r = 0 
       stats%tendencycorr%sxx_r = 0 
       stats%tendencycorr%syy_r = 0 
       stats%tendencycorr%sx_r = 0 
       stats%tendencycorr%sy_r = 0 
       stats%tendencycorr%rval_r = 0
       stats%tendencycorr%count_value = 0 

       allocate(stats%tendencycorr%value_ci(model%selectNlevs,LVT_rc%strat_nlevels))
       stats%tendencycorr%value_ci = LVT_rc%udef

       if(metric%timeOpt.eq.1) then 
          allocate(stats%tendencycorr%sxy_ts_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%tendencycorr%sxx_ts_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%tendencycorr%syy_ts_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%tendencycorr%sx_ts_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%tendencycorr%sy_ts_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%tendencycorr%rval_ts_r(LVT_rc%ngrid, model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%tendencycorr%count_value_ts(LVT_rc%ngrid, model%selectNlevs, LVT_rc%strat_nlevels))       
          stats%tendencycorr%sxy_ts_r = 0 
          stats%tendencycorr%sxx_ts_r = 0 
          stats%tendencycorr%syy_ts_r = 0 
          stats%tendencycorr%sx_ts_r = 0 
          stats%tendencycorr%sy_ts_r = 0 
          stats%tendencycorr%rval_ts_r = 0
          stats%tendencycorr%count_value_ts = 0 

          allocate(stats%tendencycorr%tavg_value_ts(LVT_rc%ngrid, &
               model%selectNlevs,LVT_rc%strat_nlevels))
          allocate(stats%tendencycorr%tavg_count_value_ts(LVT_rc%ngrid, &
               model%selectNlevs,LVT_rc%strat_nlevels))
          stats%tendencycorr%tavg_value_ts = 0.0
          stats%tendencycorr%tavg_count_value_ts=0 
          
       endif
       
       allocate(stats%tendencycorr%model_value_cval_ts(LVT_rc%ngrid,model%selectNlevs, &
            LVT_rc%strat_nlevels))
       allocate(stats%tendencycorr%model_value_pval_ts(LVT_rc%ngrid,model%selectNlevs, &
            LVT_rc%strat_nlevels))
       
       stats%tendencycorr%model_value_cval_ts = LVT_rc%udef
       stats%tendencycorr%model_value_pval_ts = LVT_rc%udef
       
       allocate(stats%tendencycorr%obs_value_cval_ts(LVT_rc%ngrid,obs%selectNlevs, &
            LVT_rc%strat_nlevels))
       allocate(stats%tendencycorr%obs_value_pval_ts(LVT_rc%ngrid,obs%selectNlevs, &
            LVT_rc%strat_nlevels))
       
       stats%tendencycorr%obs_value_cval_ts = LVT_rc%udef
       stats%tendencycorr%obs_value_pval_ts = LVT_rc%udef
       
    endif

!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1    
    if(LVT_rc%obssource(2).ne."none") then 
       metric%obsData = .true. 
    else
       metric%obsData = .false. 
    endif
    metric%stdevFlag = .false. 

  end subroutine LVT_initTendencyCorr

!BOP
! 
! !ROUTINE: LVT_diagnoseTendencyCorr
! \label{LVT_diagnoseTendencyCorr}
!
! !INTERFACE: 
  subroutine LVT_diagnoseTendencyCorr(pass)
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
!   calculating the tendency of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelTendencyCorr](\ref{diagnoseSingleModelTendencyCorr})
!     updates the tendency computation for a single variable 
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
       if(LVT_metrics%tendencycorr%selectOpt.eq.1.or.&
            LVT_metrics%tendencycorr%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleTendencyCorr(model,obs,stats,&
                  LVT_metrics%tendencycorr)
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseTendencyCorr

!BOP
! 
! !ROUTINE: diagnoseSingleTendencyCorr
! \label{diagnoseSingleTendencyCorr}
!
! !INTERFACE: 
  subroutine diagnoseSingleTendencyCorr(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the tendency computation of the 
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
    integer    :: t,k,tind
    real       :: m_diffval, o_diffval

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1.and.obs%selectNlevs.ge.1) then        
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             if(metric%selectOpt.eq.1) then
                stats%tendencycorr%model_value_pval_ts(t,k,1) = &
                     stats%tendencycorr%model_value_cval_ts(t,k,1) 
                stats%tendencycorr%model_value_cval_ts(t,k,1) = &
                     model%value(t,k)
                
                stats%tendencycorr%obs_value_pval_ts(t,k,1) = &
                     stats%tendencycorr%obs_value_cval_ts(t,k,1) 
                stats%tendencycorr%obs_value_cval_ts(t,k,1) = &
                     obs%value(t,k)
                
                if(stats%tendencycorr%model_value_cval_ts(t,k,1).ne.LVT_rc%udef.and.&
                     stats%tendencycorr%model_value_pval_ts(t,k,1).ne.LVT_rc%udef.and.&
                     stats%tendencycorr%obs_value_cval_ts(t,k,1).ne.LVT_rc%udef.and.&
                     stats%tendencycorr%obs_value_pval_ts(t,k,1).ne.LVT_rc%udef) then
                   
                   m_diffval = & 
                        (stats%tendencycorr%model_value_cval_ts(t,k,1) - &
                        stats%tendencycorr%model_value_pval_ts(t,k,1))
                   
                   o_diffval = & 
                        (stats%tendencycorr%obs_value_cval_ts(t,k,1) - &
                        stats%tendencycorr%obs_value_pval_ts(t,k,1))
                   
                   stats%tendencycorr%sxy_r(t,k,1) =  stats%tendencycorr%sxy_r(t,k,1) + & 
                        m_diffval*o_diffval
                   
                   stats%tendencycorr%sxx_r(t,k,1) =  stats%tendencycorr%sxx_r(t,k,1) + & 
                        m_diffval*m_diffval
                   
                   stats%tendencycorr%syy_r(t,k,1) =  stats%tendencycorr%syy_r(t,k,1) + & 
                        o_diffval*o_diffval
                   
                   stats%tendencycorr%sx_r(t,k,1) =  stats%tendencycorr%sx_r(t,k,1) + & 
                        m_diffval
                   
                   stats%tendencycorr%sy_r(t,k,1) =  stats%tendencycorr%sy_r(t,k,1) + & 
                        o_diffval
                   
                   stats%tendencycorr%count_value(t,k,1) = & 
                        stats%tendencycorr%count_value(t,k,1)+1

                   if(LVT_rc%strat_nlevels.gt.1) then
                      if(LVT_stats%strat_var(t,k).gt.&
                           LVT_rc%strat_var_threshold) then
                         m_diffval = & 
                              (stats%tendencycorr%model_value_cval_ts(t,k,2) - &
                              stats%tendencycorr%model_value_pval_ts(t,k,2))
                         
                         o_diffval = & 
                              (stats%tendencycorr%obs_value_cval_ts(t,k,2) - &
                              stats%tendencycorr%obs_value_pval_ts(t,k,2))
                         
                         stats%tendencycorr%sxy_r(t,k,2) =  stats%tendencycorr%sxy_r(t,k,2) + & 
                              m_diffval*o_diffval
                         
                         stats%tendencycorr%sxx_r(t,k,2) =  stats%tendencycorr%sxx_r(t,k,2) + & 
                              m_diffval*m_diffval
                         
                         stats%tendencycorr%syy_r(t,k,2) =  stats%tendencycorr%syy_r(t,k,2) + & 
                              o_diffval*o_diffval
                         
                         stats%tendencycorr%sx_r(t,k,2) =  stats%tendencycorr%sx_r(t,k,2) + & 
                              m_diffval
                         
                         stats%tendencycorr%sy_r(t,k,2) =  stats%tendencycorr%sy_r(t,k,2) + & 
                              o_diffval
                         
                         stats%tendencycorr%count_value(t,k,2) = & 
                              stats%tendencycorr%count_value(t,k,2)+1
                         
                      elseif(LVT_stats%strat_var(t,k).le.&
                           LVT_rc%strat_var_threshold) then
                         
                         m_diffval = & 
                              (stats%tendencycorr%model_value_cval_ts(t,k,3) - &
                              stats%tendencycorr%model_value_pval_ts(t,k,3))
                         
                         o_diffval = & 
                              (stats%tendencycorr%obs_value_cval_ts(t,k,3) - &
                              stats%tendencycorr%obs_value_pval_ts(t,k,3))
                         
                         stats%tendencycorr%sxy_r(t,k,3) =  stats%tendencycorr%sxy_r(t,k,3) + & 
                              m_diffval*o_diffval
                         
                         stats%tendencycorr%sxx_r(t,k,3) =  stats%tendencycorr%sxx_r(t,k,3) + & 
                              m_diffval*m_diffval
                         
                         stats%tendencycorr%syy_r(t,k,3) =  stats%tendencycorr%syy_r(t,k,3) + & 
                              o_diffval*o_diffval
                         
                         stats%tendencycorr%sx_r(t,k,3) =  stats%tendencycorr%sx_r(t,k,3) + & 
                              m_diffval
                         
                         stats%tendencycorr%sy_r(t,k,3) =  stats%tendencycorr%sy_r(t,k,3) + & 
                              o_diffval
                         
                         stats%tendencycorr%count_value(t,k,3) = & 
                              stats%tendencycorr%count_value(t,k,3)+1
                      endif
                   end if
                   if(metric%timeOpt.eq.1) then 
                      m_diffval = & 
                           (stats%tendencycorr%model_value_cval_ts(t,k,1) - &
                           stats%tendencycorr%model_value_pval_ts(t,k,1))
                      
                      o_diffval = & 
                           (stats%tendencycorr%obs_value_cval_ts(t,k,1) - &
                           stats%tendencycorr%obs_value_pval_ts(t,k,1))
                      
                      stats%tendencycorr%sxy_ts_r(t,k,1) =  stats%tendencycorr%sxy_ts_r(t,k,1) + & 
                           m_diffval*o_diffval
                      
                      stats%tendencycorr%sxx_ts_r(t,k,1) =  stats%tendencycorr%sxx_ts_r(t,k,1) + & 
                           m_diffval*m_diffval
                      
                      stats%tendencycorr%syy_ts_r(t,k,1) =  stats%tendencycorr%syy_ts_r(t,k,1) + & 
                           o_diffval*o_diffval
                      
                      stats%tendencycorr%sx_ts_r(t,k,1) =  stats%tendencycorr%sx_ts_r(t,k,1) + & 
                           m_diffval
                      
                      stats%tendencycorr%sy_ts_r(t,k,1) =  stats%tendencycorr%sy_ts_r(t,k,1) + & 
                           o_diffval
                      
                      stats%tendencycorr%count_value_ts(t,k,1) = & 
                           stats%tendencycorr%count_value_ts(t,k,1)+1
                      
                      if(LVT_rc%strat_nlevels.gt.1) then
                         if(LVT_stats%strat_var(t,k).gt.&
                              LVT_rc%strat_var_threshold) then
                            m_diffval = & 
                                 (stats%tendencycorr%model_value_cval_ts(t,k,2) - &
                                 stats%tendencycorr%model_value_pval_ts(t,k,2))
                            
                            o_diffval = & 
                                 (stats%tendencycorr%obs_value_cval_ts(t,k,2) - &
                                 stats%tendencycorr%obs_value_pval_ts(t,k,2))
                            
                            stats%tendencycorr%sxy_ts_r(t,k,2) =  stats%tendencycorr%sxy_ts_r(t,k,2) + & 
                                 m_diffval*o_diffval
                            
                            stats%tendencycorr%sxx_ts_r(t,k,2) =  stats%tendencycorr%sxx_ts_r(t,k,2) + & 
                                 m_diffval*m_diffval
                            
                            stats%tendencycorr%syy_ts_r(t,k,2) =  stats%tendencycorr%syy_ts_r(t,k,2) + & 
                                 o_diffval*o_diffval
                            
                            stats%tendencycorr%sx_ts_r(t,k,2) =  stats%tendencycorr%sx_ts_r(t,k,2) + & 
                                 m_diffval
                            
                            stats%tendencycorr%sy_ts_r(t,k,2) =  stats%tendencycorr%sy_ts_r(t,k,2) + & 
                                 o_diffval
                            
                            stats%tendencycorr%count_value_ts(t,k,2) = & 
                                 stats%tendencycorr%count_value_ts(t,k,2)+1
                            
                         elseif(LVT_stats%strat_var(t,k).le.&
                              LVT_rc%strat_var_threshold) then
                            
                            m_diffval = & 
                                 (stats%tendencycorr%model_value_cval_ts(t,k,3) - &
                                 stats%tendencycorr%model_value_pval_ts(t,k,3))
                            
                            o_diffval = & 
                                 (stats%tendencycorr%obs_value_cval_ts(t,k,3) - &
                                 stats%tendencycorr%obs_value_pval_ts(t,k,3))
                            
                            stats%tendencycorr%sxy_ts_r(t,k,3) =  stats%tendencycorr%sxy_ts_r(t,k,3) + & 
                                 m_diffval*o_diffval
                            
                            stats%tendencycorr%sxx_ts_r(t,k,3) =  stats%tendencycorr%sxx_ts_r(t,k,3) + & 
                                 m_diffval*m_diffval
                            
                            stats%tendencycorr%syy_ts_r(t,k,3) =  stats%tendencycorr%syy_ts_r(t,k,3) + & 
                                 o_diffval*o_diffval
                            
                            stats%tendencycorr%sx_ts_r(t,k,3) =  stats%tendencycorr%sx_ts_r(t,k,3) + & 
                                 m_diffval
                            
                            stats%tendencycorr%sy_ts_r(t,k,3) =  stats%tendencycorr%sy_ts_r(t,k,3) + & 
                                 o_diffval
                            
                            stats%tendencycorr%count_value_ts(t,k,3) = & 
                                 stats%tendencycorr%count_value_ts(t,k,3)+1
                            
                         endif
                      endif
                   endif
                end if
             end if
          end do
       end do
    end if
  end subroutine diagnoseSingleTendencyCorr

!BOP
! 
! !ROUTINE: LVT_computeTendencyCorr
! \label{LVT_computeTendencyCorr}
!
! !INTERFACE: 
  subroutine LVT_computeTendencyCorr(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the tendency values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelTendencyCorr](\ref{computeSingleModelTendencyCorr})
!     computes the tendency values for a single variable
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
       if(LVT_metrics%tendencycorr%selectOpt.eq.1.or.&
            LVT_metrics%tendencycorr%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%tendencycorr%timeOpt.eq.1.and.&
                  LVT_metrics%tendencycorr%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%tendencycorr%ftn_ts_loc(i),200,advance='no') &
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
             
             call computeSingleTendencyCorr(alarm,model,obs,stats,&
                  LVT_metrics%tendencycorr)
             
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%tendencycorr%timeOpt.eq.1.and.&
                  LVT_metrics%tendencycorr%extractTS.eq.1) then 
                do i=1,LVT_rc%ntslocs
                   write(LVT_metrics%tendencycorr%ftn_ts_loc(i),fmt='(a1)') ''
                enddo
             endif
          endif
       endif
    endif
  end subroutine LVT_computeTendencyCorr
  

!BOP
! 
! !ROUTINE: computeSingleTendencyCorr
! \label{computeSingleTendencyCorr}
!
! !INTERFACE: 
  subroutine computeSingleTendencyCorr(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the tendency values
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

    integer  :: t,l,k
    real     :: numer, denom, denom_sq

    if(metric%timeOpt.eq.1) then 
       if(alarm) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%tendencycorr%count_value_ts(t,k,l).ne.0.and.&
                        stats%tendencycorr%count_value_ts(t,k,l) &
                        .gt.LVT_rc%obsCountThreshold) then 
                      
                      numer = (float(stats%tendencycorr%count_value_ts(t,k,l))* &
                              stats%tendencycorr%sxy_ts_r(t,k,l) - &
                              stats%tendencycorr%sx_ts_r(t,k,l)*stats%tendencycorr%sy_ts_r(t,k,l))
                      denom_sq =  (float(stats%tendencycorr%count_value_ts(t,k,l))* &
                           stats%tendencycorr%sxx_ts_r(t,k,l)-&
                           stats%tendencycorr%sx_ts_r(t,k,l)**2)* &
                           (float(stats%tendencycorr%count_value_ts(t,k,l))*&
                           stats%tendencycorr%syy_ts_r(t,k,l)-&
                           stats%tendencycorr%sy_ts_r(t,k,l)**2)
                                               
                      if(denom_sq.gt.0) then 
                         denom = sqrt(denom_sq)
                      else
                         denom = 0.0
                      endif
                      
                      if(denom.ne.0) then 
                         stats%tendencycorr%rval_ts_r(t,k,l) = numer/denom
                         stats%tendencycorr%tavg_value_ts(t,k,l) = & 
                              stats%tendencycorr%tavg_value_ts(t,k,l) + & 
                              stats%tendencycorr%rval_ts_r(t,k,l)
                         stats%tendencycorr%tavg_count_value_ts(t,k,l) = & 
                              stats%tendencycorr%tavg_count_value_ts(t,k,l) + 1
                      else
                         stats%tendencycorr%rval_ts_r(t,k,l) = LVT_rc%udef
                      endif
                      
                   endif
                enddo
             enddo
          enddo

          if(metric%extractTS.eq.1) then 
             call LVT_writeTSinfo(metric%ftn_ts_loc,&
                  model,&
                  LVT_rc%ngrid,&
                  stats%tendencycorr%tavg_value_ts,&
                  stats%tendencycorr%tavg_count_value_ts)
          endif
          
       endif
    endif
    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%tendencycorr%count_value(t,k,l).ne.0.and.&
                        stats%tendencycorr%count_value(t,k,l) &
                        .gt.LVT_rc%obsCountThreshold) then 

                      numer = (float(stats%tendencycorr%count_value(t,k,l))* &
                           stats%tendencycorr%sxy_r(t,k,l) - &
                           stats%tendencycorr%sx_r(t,k,l)*stats%tendencycorr%sy_r(t,k,l))
                      denom_sq =  (float(stats%tendencycorr%count_value(t,k,l))* &
                           stats%tendencycorr%sxx_r(t,k,l)-&
                           stats%tendencycorr%sx_r(t,k,l)**2)* &
                           (float(stats%tendencycorr%count_value(t,k,l))*&
                           stats%tendencycorr%syy_r(t,k,l)-&
                           stats%tendencycorr%sy_r(t,k,l)**2)
                      if(denom_sq.gt.0) then 
                         denom = sqrt(denom_sq)
                      else
                         denom = 0.0
                      endif

                      if(denom.ne.0) then 
                         stats%tendencycorr%rval_r(t,k,l) = numer/denom
                      else
                         stats%tendencycorr%rval_r(t,k,l) = LVT_rc%udef
                         stats%tendencycorr%count_value(t,k,l) = LVT_rc%udef
                      endif
                   else
                      stats%tendencycorr%rval_r(t,k,l) = LVT_rc%udef
                      stats%tendencycorr%count_value(t,k,l) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
          
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%tendencycorr%rval_r(:,k,l),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%tendencycorr%value_ci(k,l))
             enddo
          enddo
       
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid,stats%tendencycorr%rval_r)      
       endif
    endif

  end subroutine computeSingleTendencyCorr

!BOP
! 
! !ROUTINE: LVT_writeMetric_TendencyCorr
! \label(LVT_writeMetric_TendencyCorr)
!
! !INTERFACE:
  subroutine LVT_writeMetric_TendencyCorr(pass,final,vlevels,stats,obs)
! 
! !USES:   
    use LVT_statsMod, only : LVT_writeSummaryStats
    use LVT_pluginIndices
    implicit none
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

    integer                 :: pass
    integer                 :: final
    integer                 :: vlevels
    type(LVT_statsEntry)    :: stats
    type(LVT_metaDataEntry) :: obs

    integer                 :: k,l,tind

    if(pass.eq.LVT_metrics%tendencycorr%npass) then 
       if(final.ne.1) then 
          if(stats%selectOpt.eq.1) then 
             do k=1,vlevels
                do l=1,LVT_rc%strat_nlevels
                   
                   call LVT_writevar_gridded(LVT_metrics%tendencycorr%ftn_ts, &
                        stats%tendencycorr%tavg_value_ts(:,k,l),stats%vid_ts(LVT_Tendencycorrid,1))
                   call LVT_writevar_gridded(LVT_metrics%tendencycorr%ftn_ts, &
                        real(stats%tendencycorr%tavg_count_value_ts(:,k,l)),&
                        stats%vid_count_ts(LVT_Tendencycorrid,1))                
                enddo
             enddo
          endif
       else
          if(LVT_metrics%tendencycorr%selectOpt.eq.1) then
             if(stats%selectOpt.eq.1) then 
                do k=1,vlevels
                   do l=1,LVT_rc%strat_nlevels
                      
                      call LVT_writevar_gridded(LVT_metrics%tendencycorr%ftn_total, &
                           stats%tendencycorr%rval_r(:,k,l),stats%vid_total(LVT_Tendencycorrid,1))
                      call LVT_writevar_gridded(LVT_metrics%tendencycorr%ftn_total, &
                           real(stats%tendencycorr%count_value(:,k,l)),&
                           stats%vid_count_total(LVT_Tendencycorrid,1))                
                   enddo
                   call LVT_writeSummaryStats(&
                        LVT_metrics%tendencycorr%ftn_summ,&
                        LVT_metrics%tendencycorr%short_name,&
                        LVT_rc%ngrid,&
                        stats%tendencycorr%rval_r(:,k,:), &
                        stats%tendencycorr%count_value(:,k,:),stats%standard_name,&
                        stats%tendencycorr%value_ci(k,:))
                enddo
             endif
          endif
       endif
    endif

  end subroutine LVT_writeMetric_TendencyCorr

!BOP
! 
! !ROUTINE: LVT_resetMetric_TendencyCorr
! \label(LVT_resetMetric_TendencyCorr)
!
! !INTERFACE:
  subroutine LVT_resetMetric_TendencyCorr(alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
    logical                :: alarm
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
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats


   call LVT_getDataStream1Ptr(model)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       
       if(stats%selectOpt.eq.1) then 
          do k=1,model%selectNlevs
             if(LVT_metrics%tendencycorr%timeOpt.eq.1) then 
                do l=1,LVT_rc%strat_nlevels
                   if(alarm) then
                      stats%tendencycorr%sxy_ts_r(:,k,l) = 0.0
                      stats%tendencycorr%sxx_ts_r(:,k,l) = 0.0
                      stats%tendencycorr%syy_ts_r(:,k,l) = 0.0
                      stats%tendencycorr%sx_ts_r(:,k,l) = 0.0
                      stats%tendencycorr%sy_ts_r(:,k,l) = 0.0
                      stats%tendencycorr%rval_ts_r(:,k,l) = 0.0
                      stats%tendencycorr%count_value_ts(:,k,l)=0 
                      
                      stats%tendencycorr%tavg_value_ts(:,k,l) = 0.0
                      stats%tendencycorr%tavg_count_value_ts(:,k,l)=0 
                   endif

                enddo
             endif
             
          enddo
          
       endif
       model => model%next
       stats => stats%next

    enddo
  end subroutine LVT_resetMetric_TendencyCorr

!BOP
! 
! !ROUTINE: LVT_writerestart_TendencyCorr
! 
! !INTERFACE:
  subroutine LVT_writerestart_TendencyCorr(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for TendencyCorr metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    integer              :: k,l
    type(LVT_metaDataEntry), pointer :: model
    type(LVT_metaDataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats

    if(LVT_metrics%tendencycorr%selectOpt.eq.1) then 
       
       call LVT_getDataStream1Ptr(model)
       call LVT_getDataStream2Ptr(obs)
       call LVT_getstatsEntryPtr(stats)

       do while(associated(model))
          
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels     
                   call LVT_writevar_restart(ftn,stats%tendencycorr%sxy_r(:,k,l))
                   call LVT_writevar_restart(ftn,stats%tendencycorr%sxx_r(:,k,l))
                   call LVT_writevar_restart(ftn,stats%tendencycorr%syy_r(:,k,l))
                   call LVT_writevar_restart(ftn,stats%tendencycorr%sx_r(:,k,l))
                   call LVT_writevar_restart(ftn,stats%tendencycorr%sy_r(:,k,l))
                   call LVT_writevar_restart(ftn,stats%tendencycorr%rval_r(:,k,l))
                   call LVT_writevar_restart(ftn,stats%tendencycorr%count_value(:,k,l))
                enddo
             enddo
          endif

          model => model%next
          obs   => obs%next
          stats => stats%next
          
       enddo
    end if
    
  end subroutine LVT_writerestart_TendencyCorr


!BOP
! 
! !ROUTINE: LVT_readrestart_TendencyCorr
! 
! !INTERFACE:
  subroutine LVT_readrestart_TendencyCorr(ftn)
! !USES: 
! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for TendencyCorr metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP

    integer              :: k,l,index
    type(LVT_metaDataEntry), pointer :: model
    type(LVT_metaDataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats

    if(LVT_metrics%tendencycorr%selectOpt.eq.1) then 
       call LVT_getDataStream1Ptr(model)
       call LVT_getDataStream2Ptr(obs)
       call LVT_getstatsEntryPtr(stats)
       
       do while(associated(model))
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   call LVT_readvar_restart(ftn,stats%tendencycorr%sxy_r(:,k,l))
                   call LVT_readvar_restart(ftn,stats%tendencycorr%sxx_r(:,k,l))
                   call LVT_readvar_restart(ftn,stats%tendencycorr%syy_r(:,k,l))
                   call LVT_readvar_restart(ftn,stats%tendencycorr%sx_r(:,k,l))
                   call LVT_readvar_restart(ftn,stats%tendencycorr%sy_r(:,k,l))
                   call LVT_readvar_restart(ftn,stats%tendencycorr%rval_r(:,k,l))
                   call LVT_readvar_restart(ftn,stats%tendencycorr%count_value(:,k,l))
                enddo
             enddo
          endif
          
          model => model%next
          obs   => obs%next
          stats => stats%next

       enddo
    end if

  end subroutine LVT_readrestart_TendencyCorr


end module LVT_TendencyCorrMod
