!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: LVT_statsDataMod
! \label(LVT_statsDataMod)
!
! !INTERFACE:
module LVT_statsDataMod
! 
! !USES:   
  use ESMF
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This module defines the data structures for storing computed
!  metrics. 
!  
!  The 'LVT_stats' object contains the list of supported
!  variables of data type 'LVT_statsEntry'. Each 'LVT_statsEntry'
!  object contains the variables for all supported variables. 
!  The memory allocation and initialization are done only for 
!  selected variables and selected metrics. 
! 
!  The 'LVT_metrics' variable contains objects to store
!  the configuration options for each metri (stored through
!  the data strcture called 'LVT_metricEntry')
! 
! !FILES USED:
!
! !REVISION HISTORY:
!  02 Oct 2008: Sujay Kumar; Initial version
! 
!EOP

  private

  type min_metric_spec 
     real,    allocatable :: model_value_total(:,:,:)     
     real,    allocatable :: model_value_ci(:,:)
     real,    allocatable :: model_value_ts(:,:,:)        
     real,    allocatable :: tavg_model_value_ts(:,:,:)        
     real,    allocatable :: model_value_asc(:,:,:)
     real,    allocatable :: model_value_adc(:,:,:)

     real,    allocatable :: obs_value_total(:,:,:)     
     real,    allocatable :: obs_value_ci(:,:)
     real,    allocatable :: obs_value_ts(:,:,:)        
     real,    allocatable :: tavg_obs_value_ts(:,:,:)        
     real,    allocatable :: obs_value_asc(:,:,:)
     real,    allocatable :: obs_value_adc(:,:,:)

     integer, allocatable :: tavg_count_model_value_ts(:,:,:)
     integer, allocatable :: tavg_count_obs_value_ts(:,:,:)
  end type min_metric_spec

  type max_metric_spec
     real,    allocatable :: model_value_total(:,:,:)     
     real,    allocatable :: model_value_ci(:,:)
     real,    allocatable :: model_value_ts(:,:,:)        
     real,    allocatable :: tavg_model_value_ts(:,:,:)        
     real,    allocatable :: model_value_asc(:,:,:)
     real,    allocatable :: model_value_adc(:,:,:)

     real,    allocatable :: obs_value_total(:,:,:)     
     real,    allocatable :: obs_value_ci(:,:)
     real,    allocatable :: obs_value_ts(:,:,:)        
     real,    allocatable :: tavg_obs_value_ts(:,:,:)        
     real,    allocatable :: obs_value_asc(:,:,:)
     real,    allocatable :: obs_value_adc(:,:,:)

     integer, allocatable :: tavg_count_model_value_ts(:,:,:)
     integer, allocatable :: tavg_count_obs_value_ts(:,:,:)
  end type max_metric_spec

  type mintime_metric_spec
     real,    allocatable :: model_value_total(:,:,:)     
     real,    allocatable :: model_value_min_total(:,:,:)     
     real,    allocatable :: model_value_ci(:,:)
     real,    allocatable :: model_value_ts(:,:,:)        
     real,    allocatable :: model_value_min_ts(:,:,:)        
     real,    allocatable :: model_value_asc(:,:,:)
     real,    allocatable :: model_value_adc(:,:,:)
     real,    allocatable :: model_value_min_asc(:,:,:)
     real,    allocatable :: model_value_min_adc(:,:,:)

     real,    allocatable :: obs_value_total(:,:,:)     
     real,    allocatable :: obs_value_min_total(:,:,:)     
     real,    allocatable :: obs_value_ci(:,:)
     real,    allocatable :: obs_value_ts(:,:,:)        
     real,    allocatable :: obs_value_min_ts(:,:,:)        
     real,    allocatable :: obs_value_asc(:,:,:)
     real,    allocatable :: obs_value_adc(:,:,:)
     real,    allocatable :: obs_value_min_asc(:,:,:)
     real,    allocatable :: obs_value_min_adc(:,:,:)

  end type mintime_metric_spec

  type maxtime_metric_spec
     real,    allocatable :: model_value_total(:,:,:)     
     real,    allocatable :: model_value_max_total(:,:,:)     
     real,    allocatable :: model_value_ci(:,:)
     real,    allocatable :: model_value_ts(:,:,:)        
     real,    allocatable :: model_value_max_ts(:,:,:)        
     real,    allocatable :: model_value_asc(:,:,:)
     real,    allocatable :: model_value_adc(:,:,:)
     real,    allocatable :: model_value_max_asc(:,:,:)
     real,    allocatable :: model_value_max_adc(:,:,:)

     real,    allocatable :: obs_value_total(:,:,:)     
     real,    allocatable :: obs_value_max_total(:,:,:)     
     real,    allocatable :: obs_value_ci(:,:)
     real,    allocatable :: obs_value_ts(:,:,:)        
     real,    allocatable :: obs_value_max_ts(:,:,:)        
     real,    allocatable :: obs_value_asc(:,:,:)
     real,    allocatable :: obs_value_adc(:,:,:)
     real,    allocatable :: obs_value_max_asc(:,:,:)
     real,    allocatable :: obs_value_max_adc(:,:,:)

  end type maxtime_metric_spec

  type sum_metric_spec
     real,    allocatable :: model_value_total(:,:,:)     
     integer, allocatable :: count_model_value_total(:,:,:)
     real,    allocatable :: model_value_ci(:,:)
     real,    allocatable :: model_value_ts(:,:,:)        
     integer, allocatable :: count_model_value_ts(:,:,:)
     real,    allocatable :: tavg_model_value_ts(:,:,:)        
     integer, allocatable :: tavg_count_model_value_ts(:,:,:)
     real,    allocatable :: model_value_asc(:,:,:)
     integer, allocatable :: count_model_value_asc(:,:,:)     
     real,    allocatable :: model_value_adc(:,:,:)
     integer, allocatable :: count_model_value_adc(:,:,:)     

     real,    allocatable :: obs_value_total(:,:,:)     
     integer, allocatable :: count_obs_value_total(:,:,:)
     real,    allocatable :: obs_value_ci(:,:)
     real,    allocatable :: obs_value_ts(:,:,:)        
     integer, allocatable :: count_obs_value_ts(:,:,:)
     real,    allocatable :: tavg_obs_value_ts(:,:,:)        
     integer, allocatable :: tavg_count_obs_value_ts(:,:,:)
     real,    allocatable :: obs_value_asc(:,:,:)
     integer, allocatable :: count_obs_value_asc(:,:,:)     
     real,    allocatable :: obs_value_adc(:,:,:)
     integer, allocatable :: count_obs_value_adc(:,:,:)     
  end type sum_metric_spec

  type mean_metric_spec
     real,    allocatable :: model_value_total(:,:,:)     
     integer, allocatable :: count_model_value_total(:,:,:)
     real,    allocatable :: model_value_ci(:,:)
     real,    allocatable :: model_value_ts(:,:,:)        
     integer, allocatable :: count_model_value_ts(:,:,:)
     real,    allocatable :: tavg_model_value_ts(:,:,:)        
     integer, allocatable :: tavg_count_model_value_ts(:,:,:)
     real,    allocatable :: model_value_asc(:,:,:)
     integer, allocatable :: count_model_value_asc(:,:,:)     
     real,    allocatable :: model_value_adc(:,:,:)
     integer, allocatable :: count_model_value_adc(:,:,:)     

     real,    allocatable :: obs_value_total(:,:,:)     
     integer, allocatable :: count_obs_value_total(:,:,:)
     real,    allocatable :: obs_value_ci(:,:)
     real,    allocatable :: obs_value_ts(:,:,:)        
     integer, allocatable :: count_obs_value_ts(:,:,:)     
     real,    allocatable :: tavg_obs_value_ts(:,:,:)        
     integer, allocatable :: tavg_count_obs_value_ts(:,:,:)     
     real,    allocatable :: obs_value_stdev_ts(:,:,:) !specified stdev of obs
     integer, allocatable :: count_obs_value_stdev_ts(:,:,:)
     real,    allocatable :: tavg_obs_value_stdev_ts(:,:,:) !specified stdev of obs
     integer, allocatable :: tavg_count_obs_value_stdev_ts(:,:,:)

     real,    allocatable :: obs_value_asc(:,:,:)
     integer, allocatable :: count_obs_value_asc(:,:,:)     
     real,    allocatable :: obs_value_adc(:,:,:)
     integer, allocatable :: count_obs_value_adc(:,:,:)     
  end type mean_metric_spec

  type tendency_metric_spec
     real,    allocatable :: model_value_total(:,:,:)     
     integer, allocatable :: count_model_value_total(:,:,:)
     real,    allocatable :: model_value_ci(:,:)
     real,    allocatable :: model_value_ts(:,:,:)        
     integer, allocatable :: count_model_value_ts(:,:,:)
     real,    allocatable :: model_value_pval_ts(:,:,:)        
     real,    allocatable :: model_value_cval_ts(:,:,:)        
     real,    allocatable :: tavg_model_value_ts(:,:,:)        
     integer, allocatable :: tavg_count_model_value_ts(:,:,:)
     real,    allocatable :: model_value_asc(:,:,:)
     integer, allocatable :: count_model_value_asc(:,:,:)     
     real,    allocatable :: model_value_adc(:,:,:)
     integer, allocatable :: count_model_value_adc(:,:,:)     

     real,    allocatable :: obs_value_total(:,:,:)     
     integer, allocatable :: count_obs_value_total(:,:,:)
     real,    allocatable :: obs_value_ci(:,:)
     real,    allocatable :: obs_value_ts(:,:,:)        
     integer, allocatable :: count_obs_value_ts(:,:,:)  
     real,    allocatable :: obs_value_pval_ts(:,:,:)        
     real,    allocatable :: obs_value_cval_ts(:,:,:)           
     real,    allocatable :: tavg_obs_value_ts(:,:,:)        
     integer, allocatable :: tavg_count_obs_value_ts(:,:,:)     

     real,    allocatable :: obs_value_asc(:,:,:)
     integer, allocatable :: count_obs_value_asc(:,:,:)     
     real,    allocatable :: obs_value_adc(:,:,:)
     integer, allocatable :: count_obs_value_adc(:,:,:)     
  end type tendency_metric_spec

  type tendcorr_metric_spec
     real,    allocatable :: model_value_pval_ts(:,:,:)        
     real,    allocatable :: model_value_cval_ts(:,:,:)        

     real,    allocatable :: obs_value_pval_ts(:,:,:)        
     real,    allocatable :: obs_value_cval_ts(:,:,:)           

     real   , allocatable :: sxx_r(:,:,:)
     real   , allocatable :: syy_r(:,:,:)
     real   , allocatable :: sxy_r(:,:,:)
     real   , allocatable :: sx_r(:,:,:)
     real   , allocatable :: sy_r(:,:,:)
     real   , allocatable :: rval_r(:,:,:)
     integer, allocatable :: count_value(:,:,:)

     real   , allocatable :: sxx_ts_r(:,:,:)
     real   , allocatable :: syy_ts_r(:,:,:)
     real   , allocatable :: sxy_ts_r(:,:,:)
     real   , allocatable :: sx_ts_r(:,:,:)
     real   , allocatable :: sy_ts_r(:,:,:)
     real   , allocatable :: rval_ts_r(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)
     
     real   , allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_ci(:,:)

  end type tendcorr_metric_spec

  type ensmean_metric_spec
     real,    allocatable :: model_value_total(:,:,:)     
     integer, allocatable :: count_model_value_total(:,:,:)
     real,    allocatable :: model_value_ci(:,:)
     real,    allocatable :: model_value_ts(:,:,:)        
     integer, allocatable :: count_model_value_ts(:,:,:)
     real,    allocatable :: tavg_model_value_ts(:,:,:)        
     integer, allocatable :: tavg_count_model_value_ts(:,:,:)
     real,    allocatable :: model_value_asc(:,:,:)
     integer, allocatable :: count_model_value_asc(:,:,:)     
     real,    allocatable :: model_value_adc(:,:,:)
     integer, allocatable :: count_model_value_adc(:,:,:)     

     real,    allocatable :: obs_value_total(:,:,:)     
     integer, allocatable :: count_obs_value_total(:,:,:)
     real,    allocatable :: obs_value_ci(:,:)
     real,    allocatable :: obs_value_ts(:,:,:)        
     integer, allocatable :: count_obs_value_ts(:,:,:)     
     real,    allocatable :: tavg_obs_value_ts(:,:,:)        
     integer, allocatable :: tavg_count_obs_value_ts(:,:,:)     
     real,    allocatable :: obs_value_stdev_ts(:,:,:) !specified stdev of obs
     integer, allocatable :: count_obs_value_stdev_ts(:,:,:)
     real,    allocatable :: tavg_obs_value_stdev_ts(:,:,:) !specified stdev of obs
     integer, allocatable :: tavg_count_obs_value_stdev_ts(:,:,:)

     real,    allocatable :: obs_value_asc(:,:,:)
     integer, allocatable :: count_obs_value_asc(:,:,:)     
     real,    allocatable :: obs_value_adc(:,:,:)
     integer, allocatable :: count_obs_value_adc(:,:,:)     
  end type ensmean_metric_spec

  type anomaly_metric_spec
     real,    allocatable :: model_value_total(:,:,:)     
     integer, allocatable :: count_model_value_total(:,:,:)
     real,    allocatable :: model_value_ci(:,:)
     real,    allocatable :: model_value_ts(:,:,:)        
     integer, allocatable :: count_model_value_ts(:,:,:)
     real,    allocatable :: tavg_model_value_ts(:,:,:)        
     integer, allocatable :: tavg_count_model_value_ts(:,:,:)
     real,    allocatable :: model_value_asc(:,:,:)
     integer, allocatable :: count_model_value_asc(:,:,:)     
     real,    allocatable :: model_value_climo(:,:,:,:)
     real,    allocatable :: obs_value_climo(:,:,:,:)
     integer, allocatable :: count_obs_value_climo(:,:,:,:)
     integer, allocatable :: count_model_value_climo(:,:,:,:)

     real,    allocatable :: obs_value_total(:,:,:)     
     integer, allocatable :: count_obs_value_total(:,:,:)
     real,    allocatable :: obs_value_ci(:,:)
     real,    allocatable :: obs_value_ts(:,:,:)        
     integer, allocatable :: count_obs_value_ts(:,:,:)
     real,    allocatable :: tavg_obs_value_ts(:,:,:)        
     integer, allocatable :: tavg_count_obs_value_ts(:,:,:)
     real,    allocatable :: obs_value_asc(:,:,:)
     integer, allocatable :: count_obs_value_asc(:,:,:)     
  end type anomaly_metric_spec

  type stdev_metric_spec
     real,    allocatable :: model_value_total(:,:,:)     
     real,    allocatable :: model_value_total_sxsx(:,:,:)     
     integer, allocatable :: count_model_value_total(:,:,:)
     real,    allocatable :: model_value_ci(:,:)
     real,    allocatable :: model_value_ts_sxsx(:,:,:) 
     real,    allocatable :: model_value_ts(:,:,:)               
     integer, allocatable :: count_model_value_ts(:,:,:)
     real,    allocatable :: tavg_model_value_ts(:,:,:)               
     integer, allocatable :: tavg_count_model_value_ts(:,:,:)
     real,    allocatable :: model_value_asc(:,:,:)
     real,    allocatable :: model_value_asc_sxsx(:,:,:)
     integer, allocatable :: count_model_value_asc(:,:,:)     
     real,    allocatable :: model_value_adc(:,:,:)
     real,    allocatable :: model_value_adc_sxsx(:,:,:)
     integer, allocatable :: count_model_value_adc(:,:,:)

     real,    allocatable :: obs_value_total(:,:,:)     
     real,    allocatable :: obs_value_total_sxsx(:,:,:)     
     integer, allocatable :: count_obs_value_total(:,:,:)
     real,    allocatable :: obs_value_ci(:,:)
     real,    allocatable :: obs_value_ts_sxsx(:,:,:)        
     real,    allocatable :: obs_value_ts(:,:,:)        
     integer, allocatable :: count_obs_value_ts(:,:,:)
     real,    allocatable :: tavg_obs_value_ts(:,:,:)        
     integer, allocatable :: tavg_count_obs_value_ts(:,:,:)
     real,    allocatable :: obs_value_asc(:,:,:)
     real,    allocatable :: obs_value_asc_sxsx(:,:,:)
     integer, allocatable :: count_obs_value_asc(:,:,:)     
     real,    allocatable :: obs_value_adc(:,:,:)
     real,    allocatable :: obs_value_adc_sxsx(:,:,:)
     integer, allocatable :: count_obs_value_adc(:,:,:)  
  
  end type stdev_metric_spec

  type ensstdev_metric_spec
     real,    allocatable :: model_value_total(:,:,:)     
     real,    allocatable :: model_value_total_sxsx(:,:,:)     
     integer, allocatable :: count_model_value_total(:,:,:)
     real,    allocatable :: model_value_ci(:,:)
     real,    allocatable :: model_value_ts_sxsx(:,:,:) 
     real,    allocatable :: model_value_ts(:,:,:)               
     integer, allocatable :: count_model_value_ts(:,:,:)
     real,    allocatable :: tavg_model_value_ts(:,:,:)               
     integer, allocatable :: tavg_count_model_value_ts(:,:,:)
     real,    allocatable :: model_value_asc(:,:,:)
     real,    allocatable :: model_value_asc_sxsx(:,:,:)
     integer, allocatable :: count_model_value_asc(:,:,:)     
     real,    allocatable :: model_value_adc(:,:,:)
     real,    allocatable :: model_value_adc_sxsx(:,:,:)
     integer, allocatable :: count_model_value_adc(:,:,:)

     real,    allocatable :: obs_value_total(:,:,:)     
     real,    allocatable :: obs_value_total_sxsx(:,:,:)     
     integer, allocatable :: count_obs_value_total(:,:,:)
     real,    allocatable :: obs_value_ci(:,:)
     real,    allocatable :: obs_value_ts_sxsx(:,:,:)        
     real,    allocatable :: obs_value_ts(:,:,:)        
     integer, allocatable :: count_obs_value_ts(:,:,:)
     real,    allocatable :: tavg_obs_value_ts(:,:,:)        
     integer, allocatable :: tavg_count_obs_value_ts(:,:,:)
     real,    allocatable :: obs_value_asc(:,:,:)
     real,    allocatable :: obs_value_asc_sxsx(:,:,:)
     integer, allocatable :: count_obs_value_asc(:,:,:)     
     real,    allocatable :: obs_value_adc(:,:,:)
     real,    allocatable :: obs_value_adc_sxsx(:,:,:)
     integer, allocatable :: count_obs_value_adc(:,:,:)  
  
  end type ensstdev_metric_spec

  type ensspread_metric_spec

     real,    allocatable :: model_value_total(:,:,:)     
     integer, allocatable :: count_model_value_total(:,:,:)
     real,    allocatable :: model_value_ci(:,:)
     real,    allocatable :: model_value_ts(:,:,:)        
     integer, allocatable :: count_model_value_ts(:,:,:)
     real,    allocatable :: tavg_model_value_ts(:,:,:)        
     integer, allocatable :: tavg_count_model_value_ts(:,:,:)
     real,    allocatable :: model_value_asc(:,:,:)
     integer, allocatable :: count_model_value_asc(:,:,:)     
     real,    allocatable :: model_value_adc(:,:,:)
     integer, allocatable :: count_model_value_adc(:,:,:)  

  end type ensspread_metric_spec
  
  type rmse_metric_spec
     real,    allocatable :: value_total(:,:,:)
     integer, allocatable :: count_value_total(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real,    allocatable :: value_ts(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)
     real,    allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)
     real,    allocatable :: value_stdev_total(:,:,:)
     real,    allocatable :: value_stdev_total_sxx(:,:,:)
     integer, allocatable :: count_value_stdev_total(:,:,:)
  end type rmse_metric_spec

  type ubrmse_metric_spec
     real,    allocatable :: value_total(:,:,:)
     real,    allocatable :: value_b1_total(:,:,:)
     real,    allocatable :: value_b2_total(:,:,:)
     integer, allocatable :: count_value_total(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real,    allocatable :: value_ts(:,:,:)
     real,    allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_b1_ts(:,:,:)
     real,    allocatable :: value_b2_ts(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)
     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)
  end type ubrmse_metric_spec

  type bias_metric_spec
     real,    allocatable :: value_total(:,:,:)
     integer, allocatable :: count_value_total(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real,    allocatable :: value_ts(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)
     real,    allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)
  end type bias_metric_spec

  type mae_metric_spec
     real,    allocatable :: value_total(:,:,:)
     integer, allocatable :: count_value_total(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real,    allocatable :: value_ts(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)
     real,    allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)
  end type mae_metric_spec

  type pody_metric_spec
     real   , allocatable :: value_total(:,:,:)
     integer, allocatable :: count_value_total(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real   , allocatable :: value_ts(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)
     real   , allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)
  end type pody_metric_spec
  
  type podn_metric_spec
     real   , allocatable :: value_total(:,:,:)
     integer, allocatable :: count_value_total(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real   , allocatable :: value_ts(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)
     real   , allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)
  end type podn_metric_spec

  type far_metric_spec
     real   , allocatable :: value_total(:,:,:)
     integer, allocatable :: count_value_total(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real   , allocatable :: value_ts(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)
     real   , allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)
  end type far_metric_spec

  type pofd_metric_spec
     real   , allocatable :: value_total(:,:,:)
     integer, allocatable :: count_value_total(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real   , allocatable :: value_ts(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)
     real   , allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)
  end type pofd_metric_spec

  type csi_metric_spec
     real   , allocatable :: value_total(:,:,:)
     real   , allocatable :: value_total_a(:,:,:)
     real   , allocatable :: value_total_b(:,:,:)
     real   , allocatable :: value_total_c(:,:,:)
     integer, allocatable :: count_value_total(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real   , allocatable :: value_ts(:,:,:)
     real   , allocatable :: value_ts_a(:,:,:)
     real   , allocatable :: value_ts_b(:,:,:)
     real   , allocatable :: value_ts_c(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)
     real   , allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)
  end type csi_metric_spec

  type acc_metric_spec
     real   , allocatable :: value_total(:,:,:)
     integer, allocatable :: count_value_total(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real   , allocatable :: value_ts(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)
     real   , allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)
  end type acc_metric_spec

  type fbias_metric_spec
     real   , allocatable :: value_total(:,:,:)
     integer, allocatable :: count_value_total(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real   , allocatable :: value_ts(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)
     real   , allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)
  end type fbias_metric_spec

  type ets_metric_spec
     real   , allocatable :: value_total(:,:,:)
     real   , allocatable :: value_total_a(:,:,:)
     real   , allocatable :: value_total_b(:,:,:)
     real   , allocatable :: value_total_c(:,:,:)
     real   , allocatable :: value_total_d(:,:,:)
     integer, allocatable :: count_value_total(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real   , allocatable :: value_ts(:,:,:)
     real   , allocatable :: value_ts_a(:,:,:)
     real   , allocatable :: value_ts_b(:,:,:)
     real   , allocatable :: value_ts_c(:,:,:)
     real   , allocatable :: value_ts_d(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)
     real   , allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)
  end type ets_metric_spec

  type rcorr_metric_spec
     real   , allocatable :: sxx_r(:,:,:)
     real   , allocatable :: syy_r(:,:,:)
     real   , allocatable :: sxy_r(:,:,:)
     real   , allocatable :: sx_r(:,:,:)
     real   , allocatable :: sy_r(:,:,:)
     real   , allocatable :: rval_r(:,:,:)
     integer, allocatable :: count_value(:,:,:)

     real   , allocatable :: sxx_ts_r(:,:,:)
     real   , allocatable :: syy_ts_r(:,:,:)
     real   , allocatable :: sxy_ts_r(:,:,:)
     real   , allocatable :: sx_ts_r(:,:,:)
     real   , allocatable :: sy_ts_r(:,:,:)
     real   , allocatable :: rval_ts_r(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)
     
     real   , allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_ci(:,:)

     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)
  end type rcorr_metric_spec
  
  type acorr_metric_spec
     real,    allocatable :: model_value_climo(:,:,:,:)
     real,    allocatable :: obs_value_climo(:,:,:,:)
     integer, allocatable :: count_obs_value_climo(:,:,:,:)
     integer, allocatable :: count_model_value_climo(:,:,:,:)
     real   , allocatable :: sxx_a(:,:,:)
     real   , allocatable :: syy_a(:,:,:)
     real   , allocatable :: sxy_a(:,:,:)
     real   , allocatable :: sx_a(:,:,:)
     real   , allocatable :: sy_a(:,:,:)
     real   , allocatable :: rval_a(:,:,:)
     real,    allocatable :: rval_a_ci(:,:)
     integer, allocatable :: count_value(:,:,:)

     real   , allocatable :: sxx_ts_a(:,:,:)
     real   , allocatable :: syy_ts_a(:,:,:)
     real   , allocatable :: sxy_ts_a(:,:,:)
     real   , allocatable :: sx_ts_a(:,:,:)
     real   , allocatable :: sy_ts_a(:,:,:)
     real   , allocatable :: rval_ts_a(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)

     real   , allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)

     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)
  end type acorr_metric_spec

  type armse_metric_spec
     real,    allocatable :: model_value_climo(:,:,:,:)
     real,    allocatable :: obs_value_climo(:,:,:,:)
     integer, allocatable :: count_obs_value_climo(:,:,:,:)
     integer, allocatable :: count_model_value_climo(:,:,:,:)
     real,    allocatable :: value_total(:,:,:)
     integer, allocatable :: count_value_total(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real,    allocatable :: value_ts(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)
     real,    allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)
  end type armse_metric_spec

  type nse_metric_spec
     real,    allocatable :: value_obs_mean(:,:,:)
     integer, allocatable :: count_value_obs_mean(:,:,:)
     real,    allocatable :: value_numer(:,:,:)
     real,    allocatable :: value_denom(:,:,:)
     real,    allocatable :: value(:,:,:)
     real,    allocatable :: value_obs_mean_ts(:,:,:)
     integer, allocatable :: count_value_obs_mean_ts(:,:,:)
     real,    allocatable :: value_numer_ts(:,:,:)
     real,    allocatable :: value_denom_ts(:,:,:)
     real,    allocatable :: value_ts(:,:,:)
     real,    allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_ci(:,:)
  end type nse_metric_spec
  
  type area_metric_spec
     real,    allocatable :: model_area(:)
     real,    allocatable :: obs_area(:)
  end type area_metric_spec

  type ensll_metric_spec
     real,    allocatable :: value_total(:,:,:)
     integer, allocatable :: count_value_total(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real,    allocatable :: value_ts(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)
     real,    allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)
  end type ensll_metric_spec

  type ensxcorr_metric_spec
     real,    allocatable :: value_param(:,:)
     real,    allocatable :: value_total(:,:,:)
     integer, allocatable :: count_value_total(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real,    allocatable :: value_ts(:,:,:)
     integer, allocatable :: count_value_ts(:,:,:)
     real,    allocatable :: tavg_value_ts(:,:,:)
     integer, allocatable :: tavg_count_value_ts(:,:,:)
  end type ensxcorr_metric_spec

  type mentropy_metric_spec
     real,    allocatable :: value_model_ts(:,:,:)
     real,    allocatable :: value_final(:,:)
     real,    allocatable :: value_ci(:)
  end type mentropy_metric_spec

  type igain_metric_spec
     real,    allocatable :: value_model_ts(:,:,:)
     real,    allocatable :: value_final(:,:)
     real,    allocatable :: value_ci(:)
  end type igain_metric_spec

  type fcomplexity_metric_spec
     real,    allocatable :: value_model_ts(:,:,:)
     real,    allocatable :: value_final(:,:)
     real,    allocatable :: value_ci(:)
  end type fcomplexity_metric_spec

  type ecomplexity_metric_spec
     real,    allocatable :: value_model_ts(:,:,:)
     real,    allocatable :: value_final(:,:)
     real,    allocatable :: value_ci(:)
  end type ecomplexity_metric_spec

  type wvt_metric_spec
     real,    allocatable :: value_model_mean_total(:,:,:)     
     integer, allocatable :: count_value_model_mean_total(:,:,:)
     real,    allocatable :: value_model_mean_ts(:,:,:)        
     integer, allocatable :: count_value_model_mean_ts(:,:,:)
     real,    allocatable :: value_obs_mean_total(:,:,:)     
     integer, allocatable :: count_value_obs_mean_total(:,:,:)
     real,    allocatable :: value_obs_mean_ts(:,:,:)        
     integer, allocatable :: count_value_obs_mean_ts(:,:,:)
     real,    allocatable :: value_mse_ts(:)
     real,    allocatable :: value_mse_total(:)
     real,    allocatable :: value_mse_pct_ts(:)
     real,    allocatable :: value_mse_pct_total(:)
  end type wvt_metric_spec

  type hn_metric_spec
     real,    allocatable :: value_total(:,:)
     real,    allocatable :: value_ts(:,:)
     integer, allocatable :: count_value_ts(:,:)     
     real,    allocatable :: tavg_value_ts(:,:)
     integer, allocatable :: tavg_count_value_ts(:,:)     
     real,    allocatable :: value_model_total(:,:,:)
     real,    allocatable :: value_obs_total(:,:,:)
     integer, allocatable :: count_value_model_total(:,:,:)
     integer, allocatable :: count_value_obs_total(:,:,:)
  end type hn_metric_spec
  
  type ensskill_metric_spec
     real,    allocatable :: value_final(:,:,:)     
     integer, allocatable :: count_value_final(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real,    allocatable :: value_ts(:,:,:)        
     integer, allocatable :: count_value_ts(:,:,:)
     real,    allocatable :: tavg_value_ts(:,:,:)        
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)     
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)     
  end type ensskill_metric_spec

  type ensme_metric_spec
     real,    allocatable :: value_final(:,:,:)     
     integer, allocatable :: count_value_final(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real,    allocatable :: value_ts(:,:,:)        
     integer, allocatable :: count_value_ts(:,:,:)
     real,    allocatable :: tavg_value_ts(:,:,:)        
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)     
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)     
  end type ensme_metric_spec

  type ensmbias_metric_spec
     real,    allocatable :: value_final(:,:,:)      
     integer, allocatable :: count_value_final(:,:,:)
     real,    allocatable :: value_ci(:,:)
     real,    allocatable :: value_ts(:,:,:)        
     integer, allocatable :: count_value_ts(:,:,:)
     real,    allocatable :: tavg_value_ts(:,:,:)        
     integer, allocatable :: tavg_count_value_ts(:,:,:)
     real,    allocatable :: value_asc(:,:,:)
     integer, allocatable :: count_value_asc(:,:,:)     
     real,    allocatable :: value_adc(:,:,:)
     integer, allocatable :: count_value_adc(:,:,:)     
  end type ensmbias_metric_spec

  type pctile_metric_spec
     integer, allocatable :: value_model_nsize(:,:)
     real,    allocatable :: value_model(:,:,:)
     integer, allocatable :: value_count_model(:,:,:)
     real,    allocatable :: tavg_value_model(:,:,:)
     integer, allocatable :: tavg_count_value_model(:,:,:)
     real,    allocatable :: value_model_ci(:,:)
  end type pctile_metric_spec

  type enspctile_metric_spec
     integer, allocatable :: value_model_nsize(:,:)
     real,    allocatable :: value_model(:,:,:)
     integer, allocatable :: value_count_model(:,:,:)
     real,    allocatable :: value_model_ci(:,:)
  end type enspctile_metric_spec

  type kstest_metric_spec
     real,    allocatable :: value_model_max(:,:)
     real,    allocatable :: value_model_min(:,:)
     real,    allocatable :: value_model_bincounts(:,:,:)
     real,    allocatable :: value_model_xrange(:,:,:)
     real,    allocatable :: value_model_cdf(:,:,:)

     real,    allocatable :: value_obs_max(:,:)
     real,    allocatable :: value_obs_min(:,:)
     real,    allocatable :: value_obs_bincounts(:,:,:)
     real,    allocatable :: value_obs_xrange(:,:,:)
     real,    allocatable :: value_obs_cdf(:,:,:)
     real,    allocatable :: value_d(:,:)
     real,    allocatable :: value_prob(:,:)
  end type kstest_metric_spec

  type, public :: LVT_statsEntry
     
     type(min_metric_spec)       :: min   
     type(max_metric_spec)       :: max   
     type(maxtime_metric_spec)   :: maxtime
     type(mintime_metric_spec)   :: mintime
     type(sum_metric_spec)       :: sum   
     type(mean_metric_spec)      :: mean  
     type(tendency_metric_spec)  :: tendency
     type(tendcorr_metric_spec)  :: tendencycorr
     type(anomaly_metric_spec)   :: anomaly
     type(stdev_metric_spec)     :: stdev   
     type(area_metric_spec)      :: area  
     
     type(rmse_metric_spec)      :: rmse  
     type(ubrmse_metric_spec)    :: ubrmse
     type(bias_metric_spec)      :: bias  
     type(mae_metric_spec)       :: mae   

     type(pody_metric_spec)      :: pody 
     type(podn_metric_spec)      :: podn
     type(far_metric_spec)       :: far
     type(pofd_metric_spec)      :: pofd
     type(csi_metric_spec)       :: csi
     type(acc_metric_spec)       :: acc
     type(fbias_metric_spec)     :: fbias
     type(ets_metric_spec)       :: ets
     type(rcorr_metric_spec)     :: rcorr
     type(acorr_metric_spec)     :: acorr
     type(armse_metric_spec)     :: armse
     type(nse_metric_spec)       :: nse

     type(ensll_metric_spec)     :: ensll
     type(ensmean_metric_spec)   :: ensmean
     type(ensmean_metric_spec)   :: ensstdev
     type(ensxcorr_metric_spec)  :: ensxcorr
     type(ensskill_metric_spec)  :: ensskill
     type(ensme_metric_spec)     :: ensme
     type(ensmbias_metric_spec)  :: ensmbias
     type(ensspread_metric_spec) :: ensspread
     
     type(mentropy_metric_spec)  :: mentropy
     type(igain_metric_spec)     :: igain
     type(fcomplexity_metric_spec)     :: fcomplexity
     type(ecomplexity_metric_spec)     :: ecomplexity
     type(wvt_metric_spec)             :: wvt
     type(hn_metric_spec)             :: hn

     type(pctile_metric_spec)         :: pctile
     type(enspctile_metric_spec)      :: enspctile
     type(kstest_metric_spec)         :: kstest


     integer          :: selectOpt    
     integer          :: computeVar
     character*100    :: short_name
     character*100    :: standard_name
     character*100    :: long_name
     character*20     :: units
     integer, allocatable :: vid_total(:,:)
     integer, allocatable :: vid_count_total(:,:)
     integer, allocatable :: vid_stdev_total(:,:)
     integer, allocatable :: vid_count_stdev_total(:,:)
     integer, allocatable :: vid_ts(:,:)
     integer, allocatable :: vid_count_ts(:,:)
     integer, allocatable :: vid_sc_total(:,:,:)
     integer, allocatable :: vid_adc_total(:,:,:)

     type(LVT_statsEntry), pointer :: next

  end type LVT_statsEntry

  type, public :: stats_struc
     
     logical          :: computeFlag
     integer, allocatable :: datamask(:,:)
     real,    allocatable :: strat_var(:,:)

     integer          :: prev_mo_tavg
     integer          :: prev_yr_tavg

     type(LVT_statsEntry) :: swnet        ! Net shortwave radiation (surface) (W/m2)
     type(LVT_statsEntry) :: lwnet        ! Net longwave radiation (surface) (W/m2)
     type(LVT_statsEntry) :: rnet         ! Net absorbed radiation (surface) (W/m2)
     type(LVT_statsEntry) :: qle          ! Latent Heat Flux (W/m2)
     type(LVT_statsEntry) :: qh           ! Sensible Heat Flux (W/m2)
     type(LVT_statsEntry) :: qg           ! Ground Heat Flux (W/m2)
     type(LVT_statsEntry) :: qf           ! Energy of fusion (W/m2)
     type(LVT_statsEntry) :: qv           ! Energy of sublimation (W/m2)
     type(LVT_statsEntry) :: qtau         ! Momentum flux (N/m2)
     type(LVT_statsEntry) :: qa           ! Advective flux (W/m2)
     type(LVT_statsEntry) :: delsurfheat  ! Change in surface heat storage (J/m2)
     type(LVT_statsEntry) :: delcoldcont  ! Change in snow water content (J/m2)
     type(LVT_statsEntry) :: br           ! Bowen Ratio
     type(LVT_statsEntry) :: ef           ! Evaporative Fraction

     type(LVT_statsEntry) :: snowf        ! Snowfall rate (kg/m2s)
     type(LVT_statsEntry) :: rainf        ! Rainfall rate (kg/m2s)
     type(LVT_statsEntry) :: evap         ! Evapotranspiration (kg/m2s)
     type(LVT_statsEntry) :: qs           ! Surface Runoff(kg/m2s)
     type(LVT_statsEntry) :: qrec         ! Recharge from river to the floodplain (kg/m2s)
     type(LVT_statsEntry) :: qsb          ! Subsurface Runoff (kg/m2s)
     type(LVT_statsEntry) :: qsm          ! Snowmelt (kg/m2s)
     type(LVT_statsEntry) :: qfz          ! Refreezing of water in the snowpack (kg/m2s)
     type(LVT_statsEntry) :: qst          ! Snow throughfall (kg/m2s)
     type(LVT_statsEntry) :: delsoilmoist ! DelSoilMoist
     type(LVT_statsEntry) :: delswe       ! DelSWE
     type(LVT_statsEntry) :: delsurfstor  ! Change in surface water storage (kg/m2)
     type(LVT_statsEntry) :: delintercept ! Change in interception storage (kg/m2)
     
     type(LVT_statsEntry) :: snowt        ! Snow surface temperature (K)
     type(LVT_statsEntry) :: vegt         ! Vegetation canopy temperature (K)
     type(LVT_statsEntry) :: baresoilt    ! Temperature of bare soil (K)
     type(LVT_statsEntry) :: avgsurft     ! Average Surface Temperature (K)
     type(LVT_statsEntry) :: radt         ! Surface Radiative Tempearture (K)
     type(LVT_statsEntry) :: albedo       ! Surface Albedo (-)
     type(LVT_statsEntry) :: swe          ! Snow water equivalent (kg/m2)
     type(LVT_statsEntry) :: sweveg       ! SWE intercepted by vegetation (kg/m2)
     type(LVT_statsEntry) :: snowfrac     ! Grid cell snow covered fraction
     type(LVT_statsEntry) :: snowdepth    ! Snow Depth(m)
     type(LVT_statsEntry) :: snowcover    ! Snow cover
     type(LVT_statsEntry) :: surfstor     ! Surface water storage (kg/m2)

     type(LVT_statsEntry) :: soilmoist
     type(LVT_statsEntry) :: soiltemp
     type(LVT_statsEntry) :: sliqfrac   ! fraction of SWE which is in the liquid phase
     type(LVT_statsEntry) :: smliqfrac  ! Average layer fraction of liquid
                                          ! moisture
     type(LVT_statsEntry) :: smfrozfrac ! Average layer fraction of liquid
                                          ! moisture
          
     type(LVT_statsEntry) :: soilwet      ! Total Soil Wetness (-)
     type(LVT_statsEntry) :: soilet       ! Plant transpiration from a particular root layer (W/m2)
     type(LVT_statsEntry) :: z0brd        ! Background (i.e., snow-free) roughness length (m)
     type(LVT_statsEntry) :: roughness    ! Roughness length (m)

     type(LVT_statsEntry) :: potevap      ! Potential Evapotranspiration (kg/m2s)
     type(LVT_statsEntry) :: ecanop       ! Interception evaporation (kg/m2s)
     type(LVT_statsEntry) :: tveg         ! Vegetation transpiration (kg/m2s)
     type(LVT_statsEntry) :: esoil        ! Bare soil evaporation (kg/m2s)
     type(LVT_statsEntry) :: ewater       ! Open water evaporation (kg/m2s)
     type(LVT_statsEntry) :: rootmoist    ! Root zone soil moisture (kg/m2)
     type(LVT_statsEntry) :: canopint     ! Total canopy water storage (kg/m2s)
     type(LVT_statsEntry) :: evapsnow     ! Snow evaporation (kg/m2s)
     type(LVT_statsEntry) :: subsnow      ! Snow sublimation (kg/m2s)
     type(LVT_statsEntry) :: subsurf      ! Sublimation of the snow free area (kg/m2s)
     type(LVT_statsEntry) :: acond        ! Aerodynamic conductance (m/s)

     type(LVT_statsEntry) :: totalprecip  ! Total precipitation rate (kg/m2/s)
     type(LVT_statsEntry) :: rainfconv  ! Convective Rainfall rate (kg/m2/s)

     type(LVT_statsEntry) :: autoresp   ! Autotrophic Respiration
     type(LVT_statsEntry) :: heteroresp ! Heterotrophic Respiration
     type(LVT_statsEntry) :: leafresp   ! Leaf Respiration
     type(LVT_statsEntry) :: npp        ! Net primary productivity in gridbox
     type(LVT_statsEntry) :: gpp        ! Gross primary productivity in gridbox
     type(LVT_statsEntry) :: nee        ! Net Ecosystem Exchange

     type(LVT_statsEntry) :: t2diag
     type(LVT_statsEntry) :: q2diag
     type(LVT_statsEntry) :: tws

     type(LVT_statsEntry) ::  windforc
     type(LVT_statsEntry) ::  rainfforc
     type(LVT_statsEntry) ::  snowfforc
     type(LVT_statsEntry) ::  tairforc
     type(LVT_statsEntry) ::  qairforc
     type(LVT_statsEntry) ::  psurfforc
     type(LVT_statsEntry) ::  swdownforc
     type(LVT_statsEntry) ::  lwdownforc
     type(LVT_statsEntry) ::  directswforc 
     type(LVT_statsEntry) ::  diffuseswforc 
     type(LVT_statsEntry) ::  nwindforc   
     type(LVT_statsEntry) ::  ewindforc   
     type(LVT_statsEntry) ::  fheightforc   
     type(LVT_statsEntry) ::  chforc   
     type(LVT_statsEntry) ::  cmforc   
     type(LVT_statsEntry) ::  emissforc   
     type(LVT_statsEntry) ::  mixratioforc   
     type(LVT_statsEntry) ::  coszenforc   
     type(LVT_statsEntry) ::  albedoforc   

     type(LVT_statsEntry) :: landmask
     type(LVT_statsEntry) :: landcover
     type(LVT_statsEntry) :: soiltype
     type(LVT_statsEntry) :: sandfrac
     type(LVT_statsEntry) :: clayfrac
     type(LVT_statsEntry) :: siltfrac
     type(LVT_statsEntry) :: porosity
     type(LVT_statsEntry) :: soilcolor
     type(LVT_statsEntry) :: elevation
     type(LVT_statsEntry) :: slope
     type(LVT_statsEntry) :: lai
     type(LVT_statsEntry) :: sai
     type(LVT_statsEntry) :: snfralbedo
     type(LVT_statsEntry) :: mxsnalbedo
     type(LVT_statsEntry) :: greenness
     type(LVT_statsEntry) :: ndvi
     type(LVT_statsEntry) :: tempbot
     type(LVT_statsEntry) :: roottemp
     type(LVT_statsEntry) :: ccond 

     type(LVT_statsEntry) :: relsmc 
     type(LVT_statsEntry) :: rhmin

     !fldas 
     type(LVT_statsEntry) :: petforc
     type(LVT_statsEntry) :: refetforc
     type(LVT_statsEntry) :: sos
     type(LVT_statsEntry) :: wrsi     
     type(LVT_statsEntry) :: kf2     
     type(LVT_statsEntry) :: sumWR     
     type(LVT_statsEntry) :: sumET     
     type(LVT_statsEntry) :: SWI   
     type(LVT_statsEntry) :: SOSa     
     type(LVT_statsEntry) :: TotalSurplusWater     
     type(LVT_statsEntry) :: MaxSurplusWater     
     type(LVT_statsEntry) :: TotalWaterDeficit     
     type(LVT_statsEntry) :: MaxWaterDeficit     
     type(LVT_statsEntry) :: TotalAETInitial     
     type(LVT_statsEntry) :: TotalWRInitial     
     type(LVT_statsEntry) :: TotalSurplusWaterInitial    
     type(LVT_statsEntry) :: TotalWaterDeficitInitial     
     type(LVT_statsEntry) :: TotalAETVeg
     type(LVT_statsEntry) :: TotalWRVeg
     type(LVT_statsEntry) :: TotalSurplusWaterVeg
     type(LVT_statsEntry) :: TotalWaterDeficitVeg
     type(LVT_statsEntry) :: TotalAETFlower
     type(LVT_statsEntry) :: TotalWRFlower
     type(LVT_statsEntry) :: TotalSurplusWaterFlower 
     type(LVT_statsEntry) :: TotalWaterDeficitFlower
     type(LVT_statsEntry) :: TotalAETRipe
     type(LVT_statsEntry) :: TotalWRRipe
     type(LVT_statsEntry) :: TotalSurplusWaterRipe
     type(LVT_statsEntry) :: TotalWaterDeficitRipe
     type(LVT_statsEntry) :: PermWiltDate 
     type(LVT_statsEntry) :: Wilting1
     type(LVT_statsEntry) :: Wilting2
     type(LVT_statsEntry) :: WRSIa
     type(LVT_statsEntry) :: growing_season

     type(LVT_statsEntry) :: ebal
     type(LVT_statsEntry) :: wbal
     type(LVT_statsEntry) :: evapbal
     type(LVT_statsEntry) :: sweoverp
     type(LVT_statsEntry) :: etoverp
     type(LVT_statsEntry) :: qsoverp
     type(LVT_statsEntry) :: qsboverp

     type(LVT_statsEntry) :: tairforc_min
     type(LVT_statsEntry) :: tairforc_max

     type(LVT_statsEntry) :: streamflow
     
     type(LVT_statsEntry) :: rtm_emissivity 
     type(LVT_statsEntry) :: rtm_tb

  end type stats_struc

  type(stats_struc),  save :: LVT_stats

  type, public :: LVT_metricEntry
     !variable to indicate if the metric is selected or not
     integer                   :: selectOpt 
     !variable to indicate if temporal computations are to be performed
     integer                   :: timeOpt   
     !variable to indicate if time series files are to be extracted (ASCII)
     integer                   :: extractTS  
     !variable to indicate if gridded time series files are to be written
     integer                   :: writeTS   
     !
     real                      :: threshold
     !variable to indicate if average seasonal cycle of the metric is to be
     !computed
     integer                   :: computeSC
     !variable to indicate if average diurnal cycle of the metric is to be
     !computed
     integer                   :: computeADC
     !short_name of the metric
     character*100             :: short_name
     !number of passes through the data record needed to compute this metric
     integer                   :: npass
     !if observation data is needed/selected in this metric computation
     logical                   :: obsData
     !if computing std deviation is enabled in this metric computation
     logical                   :: stdevFlag
     !file unit used for writing the metadata file
     integer                   :: ftn_meta_out
     !file units used for writing time series files (ASCII)
     integer, allocatable      :: ftn_ts_loc(:)
     !file unit used for writing gridded time series files
     integer                   :: ftn_ts
     !file unit used for writing the final gridded metric file
     integer                   :: ftn_total
     !file unit used for writing the summary statistics file
     integer                   :: ftn_summ
     !
     logical                   :: customNames
     !
     integer                   :: nfields
     !
     character*100,allocatable :: mName(:)
  end type LVT_metricEntry

  type, public :: metrics_struc
     type(LVT_metricEntry)   :: min
     type(LVT_metricEntry)   :: max
     type(LVT_metricEntry)   :: mintime
     type(LVT_metricEntry)   :: maxtime
     type(LVT_metricEntry)   :: sum
     type(LVT_metricEntry)   :: mean
     type(LVT_metricEntry)   :: tendency
     type(LVT_metricEntry)   :: tendencycorr
     type(LVT_metricEntry)   :: stdev

     type(LVT_metricEntry)   :: rmse
     type(LVT_metricEntry)   :: armse
     type(LVT_metricEntry)   :: bias
     type(LVT_metricEntry)   :: mae
     type(LVT_metricEntry)   :: acorr
     type(LVT_metricEntry)   :: rcorr

     type(LVT_metricEntry)   :: pody
     type(LVT_metricEntry)   :: podn
     type(LVT_metricEntry)   :: far
     type(LVT_metricEntry)   :: pofd
     type(LVT_metricEntry)   :: csi
     type(LVT_metricEntry)   :: acc
     type(LVT_metricEntry)   :: fbias
     type(LVT_metricEntry)   :: ets

     type(LVT_metricEntry)   :: area
     type(LVT_metricEntry)   :: nse
     type(LVT_metricEntry)   :: ubrmse
!ensemble metrics
     type(LVT_metricEntry)   :: ensmean
     type(LVT_metricEntry)   :: ensstdev
     type(LVT_metricEntry)   :: ensspread
     type(LVT_metricEntry)   :: ensll
     type(LVT_metricEntry)   :: ensXcorr
     type(LVT_metricEntry)   :: ensSkill
     type(LVT_metricEntry)   :: ensME
     type(LVT_metricEntry)   :: ensMeanBias
     type(LVT_metricEntry)   :: ensPercentile

     type(LVT_metricEntry)   :: mentropy
     type(LVT_metricEntry)   :: igain
     type(LVT_metricEntry)   :: fcomplexity
     type(LVT_metricEntry)   :: ecomplexity

     type(LVT_metricEntry)   :: waveletStat
     
     type(LVT_metricEntry)   :: hn
     type(LVT_metricEntry)   :: spi
     type(LVT_metricEntry)   :: sri
     type(LVT_metricEntry)   :: sswi
     type(LVT_metricEntry)   :: sgwi
     type(LVT_metricEntry)   :: percentile

     type(LVT_metricEntry)   :: rfv
     type(LVT_metricEntry)   :: anomaly
     type(LVT_metricEntry)   :: kstest

  end type metrics_struc

  type(metrics_struc) , save :: LVT_metrics

  type, public :: metricsdep
     type(LVT_metricEntry), pointer :: metricEntryPtr
  end type metricsdep

  type(metricsdep),  pointer :: LVT_metricsPtr(:)

  PUBLIC :: LVT_stats
  PUBLIC :: LVT_metrics
  PUBLIC :: LVT_metricsPtr

end module LVT_statsDataMod
