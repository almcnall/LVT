!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: LVT_datastream_pluginMod
!  \label(LVT_datastream_pluginMod)
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   defining routines that initialize various LVT-obss. 
!   The user defined functions are incorporated into 
!   the appropriate registry to be later invoked through generic calls. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  17 Feb 2004;   Sujay Kumar  Initial Specification
! 
!EOP
module LVT_datastream_pluginMod

  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LVT_datastream_plugin  
contains
!BOP
! 
! !ROUTINE: LVT_datastream_plugin
!  \label{LVT_datastream_plugin}
!
! !INTERFACE:
  subroutine LVT_datastream_plugin
    use LVT_pluginIndices
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing a new LVT-obs. 
!  The interface mandates that the following interfaces be implemented
!  and registered for each LVT-obs. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    use template_obsMod,        only : template_obsinit
    use LISoutputMod,           only : LISoutputInit
    use SCAN_obsMod,            only : SCAN_obsInit
    use LISda_obsMod,           only : LISda_obsInit
    use CEOP_obsMod,            only : ceop_obsInit
    use ISCCP_TskinobsMod,      only : ISCCP_TskinObsInit
    use MODIS_LSTobsMod,        only : MODIS_LSTobsInit
    use NASMD_obsMod,           only : NASMD_obsInit
    use SURFRAD_obsMod,         only : SURFRAD_obsInit
    use WGPBMRobsMod,           only : WGPBMRobsInit  
    use SNOTEL_obsMod,          only : SNOTEL_obsInit
    use FMISWE_obsMod,          only : FMISWE_obsInit
    use CMCSNWD_obsMod,         only : CMCSNWD_obsInit
    use SNODAS_obsMod,          only : SNODAS_obsInit
    use NASA_AMSREsm_obsMod,    only : NASA_AMSREsm_obsInit
    use LPRM_AMSREsm_obsMod,    only : LPRM_AMSREsm_obsInit
    use Ameriflux_obsMod,       only : Ameriflux_obsInit
    use ARM_obsMod,             only : ARM_obsInit
    use SMOSREX_obsMod,         only : SMOSREX_obsInit
    use AGRMET_dataMod,         only : AGRMET_dataInit
    use GlobSnow_obsMod,        only : GlobSnow_obsInit
    use SNODEP_metobsMod,       only : SNODEP_metobsInit
    use MOD10A1_obsMod,         only : MOD10A1_obsInit
    use ANSASNWD_obsMod,        only : ANSASNWD_obsInit
    use ANSASWE_obsMod,         only : ANSASWE_obsInit
    use CPCPRCP_obsMod,         only : CPCPRCP_obsInit
    use USGSSF_obsMod,          only : USGSSF_obsInit
    use NatSF_obsMod,           only : NatSF_obsInit
    use ISMN_obsMod,            only : ISMN_obsInit
    use FLUXNET_obsMod,         only : FLUXNET_obsInit
    use MOD16A2_obsMod,         only : MOD16A2_obsInit
    use UWET_obsMod,            only : UWET_obsInit
    use ARSsm_obsMod,           only : ARSsm_obsInit
    use NLDAS2_dataMod,         only : NLDAS2_dataInit
    use GHCN_obsMod,            only : GHCN_obsInit
    use ALEXI_obsMod,           only : ALEXI_obsInit
    use GRACE_obsMod,           only : GRACE_obsInit
    use USGSGWwell_obsMod,      only : USGSGWwell_obsInit
    use PBOH2O_obsMod,          only : PBOH2O_obsInit
    use SMOSL2sm_obsMod,        only : SMOSL2sm_obsInit
    use SMOSL1TB_obsMod,        only : SMOSL1TB_obsInit
    use GCOMW_AMSR2L3sm_obsMod, only : GCOMW_AMSR2L3sm_obsinit
    use SMOPSsm_obsMod,         only : SMOPSsm_obsInit
    use ESACCIsm_obsMod,        only : ESACCIsm_obsInit
    use GIMMS_NDVIobsMod,       only : GIMMS_NDVIobsinit
    use GLDAS2obsMod,           only : GLDAS2obsinit
    use MERRA2obsMod,           only : MERRA2obsinit
    use ERAinterimLandobsMod,   only : ERAinterimLandobsinit
    use SSEBop_obsMod,          only : SSEBop_obsinit
    use GRDC_obsMod,            only : GRDC_obsinit
    use GLERL_dataMod,          only : GLERL_obsinit
    use GL6JULES_obsMod,        only : GL6JULES_obsinit
    use LVTbenchmarkOUT_obsMod, only : LVTbenchmarkOUT_obsInit

    external readtemplateObs
    external readLISoutput
    external readSCANObs
    external readLISdaAsObs  
    external readCEOP
    external readISCCP_TskinObs
    external readMODIS_LSTObs
    external readNASMDObs
    external readSURFRADObs
    external readWGPBMRObs
    external readSNOTELObs
    external readFMISWEobs
    external readCMC_SNWDobs
    external readSNODASobs
    external readNASA_AMSREsmObs
    external readLPRM_AMSREsmObs
    external readAmerifluxObs
    external readARMobs
    external readSMOSREXobs
    external readAGRMETdata
    external readGlobSnowObs
    external readSNODEPmetobs
    external readMOD10A1obs
    external readANSASNWDobs
    external readANSASWEobs
    external readCPCPRCPobs
    external readUSGSSFobs
    external readNatSFobs
    external readISMNobs
    external readFLUXNETobs
    external readMOD16A2Obs
    external readUWETObs
    external readARSsmobs
    external readNLDAS2data
    external readGHCNObs
    external readALEXIobs
    external readGRACEObs
    external readUSGSGWwellobs
    external readPBOH2Oobs
    external readSMOSL2smobs
    external readSMOSL1TBobs
    external readGCOMW_AMSR2L3smobs
    external readSMOPSsmobs
    external readESACCIsmobs
    external readGLERLobs
    external readGLDAS2obs
    external readMERRA2obs
    external readERAinterimLandobs
    external readSSEBopObs
    external readGRDCobs
    external readGL6JULESobs
    external readGIMMS_NDVIobs
    external readLVTbenchmarkOUTobs

    call registerobsread(trim(LVT_LVTbenchmarkobsId)//char(0),&
         readLVTbenchmarkOUTobs)

    call registerobssetup(trim(LVT_templateobsId)//char(0), template_obsInit)
    call registerobsread(trim(LVT_templateobsId)//char(0),readtemplateObs)

    call registerobssetup(trim(LVT_LISoutputId)//char(0), LISoutputInit)
    call registerobsread(trim(LVT_LISoutputId)//char(0),readLISoutput)

    call registerobssetup(trim(LVT_SCANobsId)//char(0), SCAN_obsinit)
    call registerobsread(trim(LVT_SCANobsId)//char(0),readSCANObs)

    call registerobssetup(trim(LVT_LISdaobsId)//char(0), LISda_obsInit)
    call registerobsread(trim(LVT_LISdaobsId)//char(0),readLISdaAsObs)

    call registerobssetup(trim(LVT_ceopobsId)//char(0),CEOP_obsInit)
    call registerobsread(trim(LVT_ceopobsId)//char(0),readCEOP)

    call registerobssetup(trim(LVT_ISCCP_TskinobsId)//char(0),ISCCP_TskinobsInit)
    call registerobsread(trim(LVT_ISCCP_TskinobsId)//char(0),readISCCP_TskinObs)


    call registerobssetup(trim(LVT_NASMDobsId)//char(0), NASMD_obsinit)
    call registerobsread(trim(LVT_NASMDobsId)//char(0),readNASMDObs)

    call registerobssetup(trim(LVT_SURFRADobsId)//char(0), SURFRAD_obsinit)
    call registerobsread(trim(LVT_SURFRADobsId)//char(0),readSURFRADObs)

    call registerobssetup(trim(LVT_wgPBMRobsId)//char(0), WGPBMRobsinit)
    call registerobsread(trim(LVT_wgPBMRobsId)//char(0),readWGPBMRObs)

    call registerobssetup(trim(LVT_SNOTELobsId)//char(0), SNOTEL_obsinit)
    call registerobsread(trim(LVT_SNOTELobsId)//char(0),readSNOTELObs)

    call registerobssetup(trim(LVT_FMISWEobsId)//char(0), FMISWE_obsinit)
    call registerobsread(trim(LVT_FMISWEobsId)//char(0),readFMISWEobs)

    call registerobssetup(trim(LVT_CMCSNWDobsId)//char(0), CMCSNWD_obsinit)
    call registerobsread(trim(LVT_CMCSNWDobsId)//char(0),readCMC_SNWDobs)

    call registerobssetup(trim(LVT_SNODASobsId)//char(0), SNODAS_obsinit)
    call registerobsread(trim(LVT_SNODASobsId)//char(0),readSNODASobs)

    call registerobssetup(trim(LVT_NASAAMSREsmobsId)//char(0), NASA_AMSREsm_obsinit)
    call registerobsread(trim(LVT_NASAAMSREsmobsId)//char(0),readNASA_AMSREsmObs)

    call registerobssetup(trim(LVT_LPRMAMSREsmobsId)//char(0), LPRM_AMSREsm_obsinit)
    call registerobsread(trim(LVT_LPRMAMSREsmobsId)//char(0),readLPRM_AMSREsmObs)

    call registerobssetup(trim(LVT_AmerifluxobsId)//char(0), Ameriflux_obsinit)
    call registerobsread(trim(LVT_AmerifluxobsId)//char(0),readAmerifluxObs)

    call registerobssetup(trim(LVT_ARMobsId)//char(0), ARM_obsinit)
    call registerobsread(trim(LVT_ARMobsId)//char(0),readARMObs)

    call registerobssetup(trim(LVT_SMOSREXobsId)//char(0), SMOSREX_obsinit)
    call registerobsread(trim(LVT_SMOSREXobsId)//char(0),readSMOSREXObs)

    call registerobssetup(trim(LVT_AGRMETdataId)//char(0), AGRMET_datainit)
    call registerobsread(trim(LVT_AGRMETdataId)//char(0),readAGRMETdata)

    call registerobssetup(trim(LVT_GlobSnowObsId)//char(0), GlobSnow_obsinit)
    call registerobsread(trim(LVT_GlobSnowObsId)//char(0),readGlobSnowObs)

    call registerobssetup(trim(LVT_SNODEPmetObsId)//char(0), SNODEP_metobsinit)
    call registerobsread(trim(LVT_SNODEPmetObsId)//char(0),readSNODEPmetObs)

    call registerobssetup(trim(LVT_MOD10A1obsId)//char(0), MOD10A1_obsinit)
    call registerobsread(trim(LVT_MOD10A1obsId)//char(0),readMOD10A1obs)

    call registerobssetup(trim(LVT_ANSASNWDobsId)//char(0), ANSASNWD_obsinit)
    call registerobsread(trim(LVT_ANSASNWDobsId)//char(0),readANSASNWDobs)

    call registerobssetup(trim(LVT_ANSASWEobsId)//char(0), ANSASWE_obsinit)
    call registerobsread(trim(LVT_ANSASWEobsId)//char(0),readANSASWEobs)

    call registerobssetup(trim(LVT_CPCPRCPobsId)//char(0), CPCPRCP_obsinit)
    call registerobsread(trim(LVT_CPCPRCPobsId)//char(0),readCPCPRCPobs)

    call registerobssetup(trim(LVT_USGSSFobsId)//char(0), USGSSF_obsinit)
    call registerobsread(trim(LVT_USGSSFobsId)//char(0),readUSGSSFobs)

    call registerobssetup(trim(LVT_NatSFobsId)//char(0), NatSF_obsinit)
    call registerobsread(trim(LVT_NatSFobsId)//char(0),readNatSFobs)

    call registerobssetup(trim(LVT_ISMNobsId)//char(0), ISMN_obsinit)
    call registerobsread(trim(LVT_ISMNobsId)//char(0),readISMNobs)

    call registerobssetup(trim(LVT_FLUXNETobsId)//char(0), FLUXNET_obsinit)
    call registerobsread(trim(LVT_FLUXNETobsId)//char(0),readFLUXNETobs)

    call registerobssetup(trim(LVT_MOD16A2obsId)//char(0), MOD16A2_obsinit)
    call registerobsread(trim(LVT_MOD16A2obsId)//char(0),readMOD16A2obs)

    call registerobssetup(trim(LVT_UWETobsId)//char(0), UWET_obsinit)
    call registerobsread(trim(LVT_UWETobsId)//char(0),readUWETobs)

    call registerobssetup(trim(LVT_ARSsmobsId)//char(0), ARSsm_obsinit)
    call registerobsread(trim(LVT_ARSsmobsId)//char(0),readARSsmobs)

    call registerobssetup(trim(LVT_NLDAS2obsId)//char(0), NLDAS2_datainit)
    call registerobsread(trim(LVT_NLDAS2obsId)//char(0),readNLDAS2data)

    call registerobssetup(trim(LVT_GHCNobsId)//char(0), GHCN_obsinit)
    call registerobsread(trim(LVT_GHCNobsId)//char(0),readGHCNobs)

    call registerobssetup(trim(LVT_ALEXIobsId)//char(0), ALEXI_obsinit)
    call registerobsread(trim(LVT_ALEXIobsId)//char(0),readALEXIobs)

    call registerobssetup(trim(LVT_GRACEobsId)//char(0), GRACE_obsinit)
    call registerobsread(trim(LVT_GRACEobsId)//char(0),readGRACEObs)

    call registerobssetup(trim(LVT_USGSGWwellobsId)//char(0), USGSGWwell_obsinit)
    call registerobsread(trim(LVT_USGSGWwellobsId)//char(0),readUSGSGWwellobs)

    call registerobssetup(trim(LVT_PBOH2OobsId)//char(0), PBOH2O_obsinit)
    call registerobsread(trim(LVT_PBOH2OobsId)//char(0),readPBOH2Oobs)

    call registerobssetup(trim(LVT_SMOSL2smobsId)//char(0), SMOSL2sm_obsinit)
    call registerobsread(trim(LVT_SMOSL2smobsId)//char(0),readSMOSL2smobs)

    call registerobssetup(trim(LVT_SMOSL1TBobsId)//char(0), SMOSL1TB_obsinit)
    call registerobsread(trim(LVT_SMOSL1TBobsId)//char(0),readSMOSL1TBobs)

    call registerobssetup(trim(LVT_GCOMW_AMSR2L3smobsId)//char(0), &
         GCOMW_AMSR2L3sm_obsinit)
    call registerobsread(trim(LVT_GCOMW_AMSR2L3smobsId)//char(0),&
         readGCOMW_AMSR2L3smobs)

    call registerobssetup(trim(LVT_SMOPSsmobsId)//char(0), &
         SMOPSsm_obsinit)
    call registerobsread(trim(LVT_SMOPSsmobsId)//char(0),&
         readSMOPSsmobs)

    call registerobssetup(trim(LVT_ESACCIsmobsId)//char(0), &
         ESACCIsm_obsinit)
    call registerobsread(trim(LVT_ESACCIsmobsId)//char(0),&
         readESACCIsmobs)

    call registerobssetup(trim(LVT_GIMMS_NDVIobsId)//char(0), &
         GIMMS_NDVIobsinit)
    call registerobsread(trim(LVT_GIMMS_NDVIobsId)//char(0),&
         readGIMMS_NDVIobs)

    call registerobssetup(trim(LVT_MODIS_LSTobsId)//char(0), &
         MODIS_LSTobsinit)
    call registerobsread(trim(LVT_MODIS_LSTobsId)//char(0),&
         readMODIS_LSTobs)

    call registerobssetup(trim(LVT_GLDAS2obsId)//char(0), &
         GLDAS2obsinit)
    call registerobsread(trim(LVT_GLDAS2obsId)//char(0),&
         readGLDAS2obs)

    call registerobssetup(trim(LVT_MERRA2obsId)//char(0), &
         MERRA2obsinit)
    call registerobsread(trim(LVT_MERRA2obsId)//char(0),&
         readMERRA2obs)

    call registerobssetup(trim(LVT_ERAIlandobsId)//char(0), &
         ERAinterimLandobsinit)
    call registerobsread(trim(LVT_ERAIlandobsId)//char(0),&
         readERAinterimLandobs)

    call registerobssetup(trim(LVT_SSEBopobsId)//char(0), &
         SSEBop_obsinit)
    call registerobsread(trim(LVT_SSEBopobsId)//char(0),&
         readSSEBopObs)

    call registerobssetup(trim(LVT_GRDCobsId)//char(0), &
         GRDC_obsinit)
    call registerobsread(trim(LVT_GRDCobsId)//char(0),&
         readGRDCObs)

    call registerobssetup(trim(LVT_GLERLobsId)//char(0), &
         GLERL_obsinit)
    call registerobsread(trim(LVT_GLERLobsId)//char(0),&
         readGLERLobs)

    call registerobssetup(trim(LVT_GL6JULESobsId)//char(0), &
         GL6JULES_obsinit)
    call registerobsread(trim(LVT_GL6JULESobsId)//char(0),&
         readGL6JULESObs)

    call registerobssetup(trim(LVT_LVTbenchmarkobsId)//char(0), &
         LVTbenchmarkOUT_obsInit)
    call registerobsread(trim(LVT_LVTbenchmarkobsId)//char(0),&
         readLVTbenchmarkOUTobs)

  end subroutine LVT_datastream_plugin
end module LVT_datastream_pluginMod
