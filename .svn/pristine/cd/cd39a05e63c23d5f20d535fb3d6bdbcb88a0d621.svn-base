!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !ROUTINE: LVT_readMetricsAttributes
!  \label{LVT_readMetricsAttributes}
!
! !INTERFACE: 
subroutine LVT_readMetricsAttributes(attribFile)
! 
! !USES: 
  use ESMF
  use LVT_statsDataMod, only : LVT_metrics
  use LVT_logMod,       only : LVT_verify

  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in) :: attribFile
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the attributes of different analysis
!  metrics within LVT. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  type(ESMF_Config)            :: attribConfig
  integer                      :: rc

  attribConfig = ESMF_ConfigCreate(rc=rc)
  call ESMF_ConfigLoadFile(attribConfig,trim(attribFile),rc=rc)
  call LVT_verify(rc,'loading file '//trim(attribFile)//' failed')


  call ESMF_ConfigFindLabel(attribConfig,"Mean:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%mean,"MEAN",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Anomaly:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%anomaly,"Anomaly",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Min:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%min,"MIN",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Max:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%max,"MAX",rc)

  call ESMF_ConfigFindLabel(attribConfig,"MinTime:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%mintime,"MINTIME",rc)

  call ESMF_ConfigFindLabel(attribConfig,"MaxTime:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%maxtime,"MAXTIME",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Sum:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%sum,"SUM",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Standard deviation:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%stdev,"STDEV",rc)

  call ESMF_ConfigFindLabel(attribConfig,"RMSE:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%rmse,"RMSE",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Anomaly RMSE:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%armse,"ARMSE",rc)
  
  call ESMF_ConfigFindLabel(attribConfig,"Bias:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%bias,"BIAS",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Mean absolute error:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%mae,"MAE",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Anomaly correlation:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%acorr,"ACORR",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Raw correlation:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%rcorr,"RCORR",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Probability of detection (PODy):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%pody,"PODY",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Probability of detection (PODn):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%podn,"PODN",rc)

  call ESMF_ConfigFindLabel(attribConfig,"False alarm ratio (FAR):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%far,"FAR",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Probability of false detection (POFD):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%pofd,"POFD",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Critical success index (CSI):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%csi,"CSI",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Accuracy measure (ACC):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%acc,"ACC",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Frequency bias (FBIAS):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%fbias,"FBIAS",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Equitable threat score (ETS):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ets,"ETS",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Area metric:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%area,"AREA",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Nash sutcliffe efficiency:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%nse,"NSE",rc)

  call ESMF_ConfigFindLabel(attribConfig,"ubRMSE:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ubrmse,"ubRMSE",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Ensemble mean:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ensmean,"EnsMEAN",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Ensemble standard deviation:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ensstdev,"EnsStdev",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Ensemble likelihood:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ensll,"EnsLL",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Ensemble cross correlation:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ensxcorr,"EnsXcorr",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Ensemble skill:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ensSkill,"EnsSkill",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Ensemble mean error:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ensME,"EnsME",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Ensemble mean bias:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ensMeanBias,"EnsMbias",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Ensemble spread:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ensspread,"EnsSpread",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Ensemble percentile:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%enspercentile,"EnsPercentile",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Metric entropy:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%mentropy,"Mentropy",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Information gain:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%igain,"Igain",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Fluctuation complexity:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%fcomplexity,"Fcomplexity",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Effective complexity:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ecomplexity,"Ecomplexity",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Wavelet stat:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%waveletStat,"Waveletstat",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Hausdorff norm:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%hn,"Hnorm",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Standard precipitation index:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%spi,"SPI",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Standard runoff index:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%sri,"SRI",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Standardized soil water index:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%sswi,"SSWI",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Standardized ground water index:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%sgwi,"SGWI",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Percentile:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%percentile,&
       "Percentile",rc)

  call ESMF_ConfigFindLabel(attribConfig,"River flow variate:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%rfv,&
       "RFV",rc)

  call ESMF_ConfigFindLabel(attribConfig,"K-S test:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%kstest,&
       "KStest",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Tendency:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%tendency,"TENDENCY",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Tendency correlation:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%tendencycorr,"TENDENCYCORR",rc)

end subroutine LVT_readMetricsAttributes

!BOP
! 
! !ROUTINE: get_metric_attributes
! \label{get_metric_attributes}
!
! !INTERFACE:
subroutine get_metric_attributes(attribConfig,attribEntry,&
     short_name, status)
! 
! !USES: 
  use ESMF
  use LVT_statsDataMod
  use LVT_logMod, only : LVT_verify

  implicit none
!
! !INPUT PARAMETERS: 
  type(ESMF_Config),      intent(inout) :: attribConfig
  type(LVT_metricEntry), intent(inout) :: attribEntry
  character(len=*),        intent(in)   :: short_name
  integer,                 intent(in)   :: status
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine reads the run time specification for a specific
!  analysis metric. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  integer                               :: rc

  if(status.eq.0) then 
     call ESMF_ConfigGetAttribute(attribConfig,attribEntry%selectOpt,&
          default=0,rc=rc)
     call LVT_verify(rc,'Error reading select option for '//trim(short_name))
     call ESMF_ConfigGetAttribute(attribConfig,attribEntry%timeOpt,&
          default=0,rc=rc)
     call LVT_verify(rc,'Error reading in-time option for '//trim(short_name))
     call ESMF_ConfigGetAttribute(attribConfig,attribEntry%writeTS,&
          default=0,rc=rc)
     call LVT_verify(rc,'Error reading write time series files option for '//trim(short_name))
     call ESMF_ConfigGetAttribute(attribConfig,attribEntry%extractTS,&
          default=0,rc=rc)
     call LVT_verify(rc,'Error reading extract time series for '//trim(short_name))
     call ESMF_ConfigGetAttribute(attribConfig,attribEntry%threshold,&
          default=0.0,rc=rc)
     call LVT_verify(rc,'Error reading threshold option for '//trim(short_name))
     call ESMF_ConfigGetAttribute(attribConfig,attribEntry%computeSC,&
          default=0,rc=rc)
     call LVT_verify(rc,'Error reading compute seasonal cycle option for '//trim(short_name))
     call ESMF_ConfigGetAttribute(attribConfig,attribEntry%computeADC,&
          default=0,rc=rc)
     call LVT_verify(rc,'Error reading compute average diurnal cycle option for '//trim(short_name))
     attribEntry%short_name  = (short_name)
  else

     attribEntry%selectOpt = 0
     attribEntry%timeOpt = 0 
     attribEntry%short_name = (short_name)

  endif
end subroutine get_metric_attributes
