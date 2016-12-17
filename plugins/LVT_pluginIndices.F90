!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: LVT_pluginIndices
!  \label(LVT_pluginIndices)
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
!  !DESCRIPTION: 
!   The code in this file provides values of indices used to 
!   to register functions in the plugin modules
!
!   The index definitions are simply a convention
!   The user may change these options, and the lis.config 
!   should be changed appropriately to ensure that the correct function
!   is called at run time
! 
!   NOTES: The indices for metrics should be in increasing order, whereas
!   the indices for other plugin sets need not be. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  23 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP
module LVT_pluginIndices

  PRIVATE
   
!BOC
!-------------------------------------------------------------------------
! supported metrics 
!------------------------------------------------------------------------- 
  integer, public,  parameter :: LVT_METRIC_SINDEX = 1
  integer, public,  parameter :: LVT_MEANid        = 1
  integer, public,  parameter :: LVT_Minid        = 2
  integer, public,  parameter :: LVT_Maxid        = 3
  integer, public,  parameter :: LVT_Sumid        = 4
  integer, public,  parameter :: LVT_Stdevid      = 5
  integer, public,  parameter :: LVT_RMSEid        = 6
  integer, public,  parameter :: LVT_BIASid        = 7
  integer, public,  parameter :: LVT_MAEid         = 8
  integer, public,  parameter :: LVT_PODYid        = 9
  integer, public,  parameter :: LVT_PODNid        = 10
  integer, public,  parameter :: LVT_FARid         = 11
  integer, public,  parameter :: LVT_POFDid        = 12
  integer, public,  parameter :: LVT_CSIid         = 13
  integer, public,  parameter :: LVT_ACCid         = 14
  integer, public,  parameter :: LVT_FBIASid       = 15
  integer, public,  parameter :: LVT_ETSid         = 16
  integer, public,  parameter :: LVT_Rcorrid       = 17
  integer, public,  parameter :: LVT_Acorrid       = 18
  integer, public,  parameter :: LVT_ARMSEid       = 19
  integer, public,  parameter :: LVT_NSEid         = 20
  integer, public,  parameter :: LVT_AREAid        = 21
  integer, public,  parameter :: LVT_ubRMSEid      = 22
  integer, public,  parameter :: LVT_waveletStatId = 23
  integer, public,  parameter :: LVT_hnId          = 24
  integer, public,  parameter :: LVT_spiId         = 25
  integer, public,  parameter :: LVT_sriId         = 26
  integer, public,  parameter :: LVT_sswiId        = 27
  integer, public,  parameter :: LVT_sgwiId        = 28
  integer, public,  parameter :: LVT_percentileId  = 29
  integer, public,  parameter :: LVT_RFVid         = 30
  integer, public,  parameter :: LVT_Anomalyid     = 31
  integer, public,  parameter :: LVT_KStestid      = 32
  integer, public,  parameter :: LVT_MinTimeid    = 33
  integer, public,  parameter :: LVT_MaxTimeid    = 34
  integer, public,  parameter :: LVT_Tendencyid    = 35
  integer, public,  parameter :: LVT_TendencyCorrid    = 36
  integer, public,  parameter :: LVT_METRIC_EINDEX = 36

!ensemble metrics
  integer, public,  parameter :: LVT_ENSMETRIC_SINDEX = 37
  integer, public,  parameter :: LVT_EnsMEANid        = 37
  integer, public,  parameter :: LVT_EnsStdevid       = 38
  integer, public,  parameter :: LVT_EnsSpreadid      = 39
  integer, public,  parameter :: LVT_EnsLLid          = 40
  integer, public,  parameter :: LVT_EnsXcorrid       = 41
  integer, public,  parameter :: LVT_EnsSkillid       = 42
  integer, public,  parameter :: LVT_EnsMEid          = 43
  integer, public,  parameter :: LVT_EnsMeanBiasid    = 44
  integer, public,  parameter :: LVT_EnsPercentileid  = 45
  integer, public,  parameter :: LVT_ENSMETRIC_EINDEX = 45
!Information content metrics
  integer, public,  parameter :: LVT_ICMETRIC_SINDEX = 46
  integer, public,  parameter :: LVT_mentropyid      = 46
  integer, public,  parameter :: LVT_igainid         = 47
  integer, public,  parameter :: LVT_fcomplexityid   = 48
  integer, public,  parameter :: LVT_ecomplexityid   = 49
  integer, public,  parameter :: LVT_ICMETRIC_EINDEX = 49

  integer, public,  parameter :: LVT_NMETRICS        = 49

!-------------------------------------------------------------------------
! Run modes
!------------------------------------------------------------------------- 
   character*50, public,  parameter :: LVT_DataCompId = "Data intercomparison"
   character*50, public,  parameter :: LVT_dastatId = "DA statistics processing"
   character*50, public,  parameter :: LVT_benchMarkId = "Benchmarking"
   character*50, public,  parameter :: LVT_daobsId = "DA observation processing"
   character*50, public,  parameter :: LVT_optUEId  = "OPTUE output processing"
   character*50, public,  parameter :: LVT_rtmrunId  = "RTM output processing"
!-------------------------------------------------------------------------
! Domains
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LVT_latlonId    = "latlon"
   character*50, public,  parameter :: LVT_mercId      = "mercator"
   character*50, public,  parameter :: LVT_lambertId   = "lambert"
   character*50, public,  parameter :: LVT_gaussId     = "gaussian"
   character*50, public,  parameter :: LVT_polarId     = "polar"
   character*50, public,  parameter :: LVT_utmId       = "UTM"
!-------------------------------------------------------------------------
! Observations
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LVT_templateobsId      = "none" 
   character*50, public,  parameter :: LVT_LISoutputId        = "LIS output"
   character*50, public,  parameter :: LVT_LISdaobsId         = "LIS DAOBS"
   character*50, public,  parameter :: LVT_ceopobsId          = "CEOP"
   character*50, public,  parameter :: LVT_ISCCP_TskinobsId   = "ISCCP LST"
   character*50, public,  parameter :: LVT_MODIS_LSTobsId     = "MODIS LST"
   character*50, public,  parameter :: LVT_SCANobsId          = "SCAN"
   character*50, public,  parameter :: LVT_NASMDobsId         = "NASMD"
   character*50, public,  parameter :: LVT_ISMNobsId          = "ISMN"
   character*50, public,  parameter :: LVT_SURFRADobsId       = "SURFRAD"
   character*50, public,  parameter :: LVT_wgPBMRobsId        = "WG PBMRsm"
   character*50, public,  parameter :: LVT_SNOTELobsId        = "SNOTEL"
   character*50, public,  parameter :: LVT_LSWG_TbobsId       = "LSWG Tb"
   character*50, public,  parameter :: LVT_FMISWEobsId        = "FMI SWE"
   character*50, public,  parameter :: LVT_CMCSNWDobsId       = "CMC"
   character*50, public,  parameter :: LVT_SNODASobsId        = "SNODAS"
   character*50, public,  parameter :: LVT_NASAAMSREsmobsId   = "AMSR-E NASA soil moisture"
   character*50, public,  parameter :: LVT_LPRMAMSREsmobsId   = "AMSR-E LPRM soil moisture"
   character*50, public,  parameter :: LVT_AMMAobsId          = "AMMA"
   character*50, public,  parameter :: LVT_AmerifluxobsId     = "Ameriflux"
   character*50, public,  parameter :: LVT_ARMobsId           = "ARM"
   character*50, public,  parameter :: LVT_SMOSREXobsId       = "SMOSREX"
   character*50, public,  parameter :: LVT_AGRMETdataId       = "AGRMET"
   character*50, public,  parameter :: LVT_GlobSnowObsId      = "Globsnow"
   character*50, public,  parameter :: LVT_SNODEPmetobsId        = "SNODEP metobs"
   character*50, public,  parameter :: LVT_MOD10A1obsId       = "MOD10A1"
   character*50, public,  parameter :: LVT_ANSASNWDobsId      = "ANSA snow depth"
   character*50, public,  parameter :: LVT_ANSASWEobsId       = "ANSA SWE"
   character*50, public,  parameter :: LVT_CPCPRCPobsId       = "CPC precipitation"
   character*50, public,  parameter :: LVT_USGSSFobsId        = "USGS streamflow"
   character*50, public,  parameter :: LVT_NatSFobsId        = "Naturalized streamflow"
   character*50, public,  parameter :: LVT_FLUXNETobsId       = "FLUXNET"
   character*50, public,  parameter :: LVT_MOD16A2obsId       = "MOD16A2"
   character*50, public,  parameter :: LVT_UWETobsId          = "UW ET"
   character*50, public,  parameter :: LVT_ARSsmobsId         = "USDA ARS soil moisture"  
   character*50, public,  parameter :: LVT_NLDAS2obsId        = "NLDAS2"
   character*50, public,  parameter :: LVT_GHCNobsId          = "GHCN"
   character*50, public,  parameter :: LVT_ALEXIobsId         = "ALEXI"
   character*50, public,  parameter :: LVT_GRACEobsId         = "GRACE"
   character*50, public,  parameter :: LVT_USGSGWwellobsId    = "USGS ground water well data"
   character*50, public,  parameter :: LVT_PBOH2OobsId        = "PBO H2O"
   character*50, public,  parameter :: LVT_SMOSL2smobsId      = "SMOS L2 soil moisture"
   character*50, public,  parameter :: LVT_SMOSL1TBobsId      = "SMOS L1 TB"
   character*50, public,  parameter :: LVT_GCOMW_AMSR2L3smobsId  = "GCOMW AMSR2 L3 soil moisture"
   character*50, public,  parameter :: LVT_SMOPSsmobsId  = "SMOPS soil moisture"
   character*50, public,  parameter :: LVT_ESACCIsmobsId = "ESA CCI soil moisture"
   character*50, public,  parameter :: LVT_GIMMS_NDVIobsId = "GIMMS NDVI"
   character*50, public,  parameter :: LVT_GLDAS2obsId = "GLDAS2"
   character*50, public,  parameter :: LVT_MERRA2obsId = "MERRA2"
   character*50, public,  parameter :: LVT_ERAIlandobsId = "ERA interim land"
   character*50, public,  parameter :: LVT_SSEBopobsId = "SSEB"
   character*50, public,  parameter :: LVT_GRDCobsId = "GRDC"
   character*50, public,  parameter :: LVT_GLERLobsId = "GLERL hydro data"
   character*50, public,  parameter :: LVT_GL6JULESobsId = "GL6 JULES data"
   character*50, public,  parameter :: LVT_LVTbenchmarkobsId = "LVT benchmark"
!-------------------------------------------------------------------------
! Training algorithms
!------------------------------------------------------------------------- 
   character*50, public,  parameter :: LVT_LinearRegressionId = "Linear regression"

!EOC
 end module LVT_pluginIndices
