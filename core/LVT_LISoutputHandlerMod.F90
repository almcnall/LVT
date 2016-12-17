!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"

!BOP
! 
! !MODULE: LVT_LISoutputHandlerMod
! \label(LVT_LISoutputHandlerMod)
!
! !INTERFACE:
module LVT_LISoutputHandlerMod
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_constantsMod
  use grib_api
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  public :: LVT_LISoutputInit
  
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
! MOC - Model Output Convention
!-----------------------------------------------------------------------------
  public :: LVT_LISoutput
  public :: LVT_resetLISoutputContainers

  integer :: LVT_LIS_MOC_SWNET(2)      = -9999
  integer :: LVT_LIS_MOC_LWNET(2)      = -9999
  integer :: LVT_LIS_MOC_QLE(2)        = -9999
  integer :: LVT_LIS_MOC_QH(2)         = -9999
  integer :: LVT_LIS_MOC_QG(2)         = -9999
  integer :: LVT_LIS_MOC_QF(2)         = -9999
  integer :: LVT_LIS_MOC_QV(2)         = -9999
  integer :: LVT_LIS_MOC_QTAU(2)       = -9999
  integer :: LVT_LIS_MOC_QA(2)         = -9999
  integer :: LVT_LIS_MOC_DELSURFHEAT(2) = -9999
  integer :: LVT_LIS_MOC_DELCOLDCONT(2) = -9999

   ! ALMA WATER BALANCE COMPONENTS
  integer :: LVT_LIS_MOC_SNOWF(2)      = -9999
  integer :: LVT_LIS_MOC_RAINF(2)      = -9999
  integer :: LVT_LIS_MOC_EVAP(2)       = -9999
  integer :: LVT_LIS_MOC_QS(2)         = -9999
  integer :: LVT_LIS_MOC_QREC(2)       = -9999
  integer :: LVT_LIS_MOC_QSB(2)        = -9999
  integer :: LVT_LIS_MOC_QSM(2)        = -9999
  integer :: LVT_LIS_MOC_QFZ(2)        = -9999
  integer :: LVT_LIS_MOC_QST(2)        = -9999
  integer :: LVT_LIS_MOC_DELSOILMOIST(2) = -9999
  integer :: LVT_LIS_MOC_DELSWE(2)     = -9999
  integer :: LVT_LIS_MOC_DELSURFSTOR(2)  = -9999
  integer :: LVT_LIS_MOC_DELINTERCEPT(2) = -9999

   ! ALMA SURFACE STATE VARIABLES
  integer :: LVT_LIS_MOC_SNOWT(2)      = -9999
  integer :: LVT_LIS_MOC_VEGT(2)       = -9999
  integer :: LVT_LIS_MOC_BARESOILT(2)  = -9999
  integer :: LVT_LIS_MOC_AVGSURFT(2)   = -9999
  integer :: LVT_LIS_MOC_RADT(2)       = -9999
  integer :: LVT_LIS_MOC_ALBEDO(2)     = -9999
  integer :: LVT_LIS_MOC_SWE(2)        = -9999
  integer :: LVT_LIS_MOC_SWEVEG(2)     = -9999 
  integer :: LVT_LIS_MOC_SURFSTOR(2)   = -9999

   ! ALMA SUBSURFACE STATE VARIABLES
   integer :: LVT_LIS_MOC_SOILMOIST(2)  = -9999
   integer :: LVT_LIS_MOC_SOILTEMP(2)   = -9999
   integer :: LVT_LIS_MOC_SMLIQFRAC(2)  = -9999
   integer :: LVT_LIS_MOC_SMFROZFRAC(2) = -9999
   integer :: LVT_LIS_MOC_SOILWET(2)    = -9999

   ! ALMA EVAPORATION COMPONENTS
   integer :: LVT_LIS_MOC_POTEVAP(2)    = -9999
   integer :: LVT_LIS_MOC_ECANOP(2)     = -9999
   integer :: LVT_LIS_MOC_TVEG(2)       = -9999
   integer :: LVT_LIS_MOC_ESOIL(2)      = -9999
   integer :: LVT_LIS_MOC_EWATER(2)     = -9999
   integer :: LVT_LIS_MOC_ROOTMOIST(2)  = -9999
   integer :: LVT_LIS_MOC_CANOPINT(2)   = -9999
   integer :: LVT_LIS_MOC_EVAPSNOW(2)   = -9999
   integer :: LVT_LIS_MOC_SUBSNOW(2)    = -9999
   integer :: LVT_LIS_MOC_SUBSURF(2)    = -9999
   integer :: LVT_LIS_MOC_ACOND(2)      = -9999

   ! ALMA OTHER HYDROLOGIC VARIABLES
  integer :: LVT_LIS_MOC_WATERTABLED(2)= -9999
  integer :: LVT_LIS_MOC_TWS(2)        = -9999

   ! ALMA COLD SEASON PROCESSES
  integer :: LVT_LIS_MOC_SNOWCOVER(2)  = -9999
  integer :: LVT_LIS_MOC_SALBEDO(2)    = -9999
  integer :: LVT_LIS_MOC_SNOWTPROF(2)  = -9999
  integer :: LVT_LIS_MOC_SNOWDEPTH(2)  = -9999
  integer :: LVT_LIS_MOC_SLIQFRAC(2)   = -9999

   ! ALMA VARIABLES TO BE COMPARED WITH REMOTE SENSED DATA
   integer :: LVT_LIS_MOC_LWUP(2)       = -9999

   ! ALMA CARBON VARIABLES
   integer :: LVT_LIS_MOC_GPP(2)        = -9999
   integer :: LVT_LIS_MOC_NPP(2)        = -9999
   integer :: LVT_LIS_MOC_NEE(2)        = -9999
   integer :: LVT_LIS_MOC_AUTORESP(2)   = -9999
   integer :: LVT_LIS_MOC_HETERORESP(2) = -9999
   integer :: LVT_LIS_MOC_LEAFRESP(2)   = -9999
   integer :: LVT_LIS_MOC_TOTSOILCARB(2)= -9999
   integer :: LVT_LIS_MOC_TOTLIVBIOM(2) = -9999

   ! ALMA FORCING VARIABLES
   integer :: LVT_LIS_MOC_WINDFORC(2)   = -9999
   integer :: LVT_LIS_MOC_RAINFFORC(2)  = -9999
   integer :: LVT_LIS_MOC_SNOWFFORC(2)  = -9999
   integer :: LVT_LIS_MOC_CRAINFFORC(2) = -9999
   integer :: LVT_LIS_MOC_TAIRFORC(2)   = -9999
   integer :: LVT_LIS_MOC_QAIRFORC(2)   = -9999
   integer :: LVT_LIS_MOC_PSURFFORC(2)  = -9999
   integer :: LVT_LIS_MOC_SWDOWNFORC(2) = -9999
   integer :: LVT_LIS_MOC_LWDOWNFORC(2) = -9999

   ! CLSM FORCING VARIABLES
   integer :: LVT_LIS_MOC_PARDRFORC(2)  = -9999
   integer :: LVT_LIS_MOC_PARDFFORC(2)  = -9999

   ! PARAMETER OUTPUT - EXPERIMENTAL (USE W/WRF-WPS)
   integer :: LVT_LIS_MOC_LANDMASK(2)   = -9999
   integer :: LVT_LIS_MOC_LANDCOVER(2)  = -9999
   integer :: LVT_LIS_MOC_SOILTYPE(2)   = -9999
   integer :: LVT_LIS_MOC_SANDFRAC(2)   = -9999
   integer :: LVT_LIS_MOC_CLAYFRAC(2)   = -9999
   integer :: LVT_LIS_MOC_SILTFRAC(2)   = -9999
   integer :: LVT_LIS_MOC_POROSITY(2)   = -9999
   integer :: LVT_LIS_MOC_SOILCOLOR(2)  = -9999
   integer :: LVT_LIS_MOC_ELEVATION(2)  = -9999
   integer :: LVT_LIS_MOC_SLOPE(2)      = -9999
   integer :: LVT_LIS_MOC_LAI(2)        = -9999
   integer :: LVT_LIS_MOC_SAI(2)        = -9999
   integer :: LVT_LIS_MOC_SNFRALBEDO(2) = -9999
   integer :: LVT_LIS_MOC_MXSNALBEDO(2) = -9999
   integer :: LVT_LIS_MOC_GREENNESS(2)  = -9999
   integer :: LVT_LIS_MOC_TEMPBOT(2)   = -9999

   ! NLDAS OUTPUT
   integer :: LVT_LIS_MOC_CCOND(2)    = -9999

   ! ADDITIONAL AFWA VARIABLES
   integer :: LVT_LIS_MOC_RELSMC(2)       = -9999
   integer :: LVT_LIS_MOC_RHMIN(2)        = -9999
   integer :: LVT_LIS_MOC_ROOTTEMP(2)  = -9999
   integer :: LVT_LIS_MOC_TOTALPRECIP(2) = -9999
   integer :: LVT_LIS_MOC_RAINFCONV(2) = -9999

   ! multivariate diagnostics (Bowen Ratio, Evaporative fraction)
   integer :: LVT_LIS_MOC_BR(2) = -9999
   integer :: LVT_LIS_MOC_EF(2) = -9999

   ! ADDITIONAL COUPLING FORCING VARIABLES
   integer :: LVT_LIS_MOC_DIRECTSWFORC(2)  = -9999
   integer :: LVT_LIS_MOC_DIFFUSESWFORC(2) = -9999
   integer :: LVT_LIS_MOC_NWINDFORC(2)     = -9999
   integer :: LVT_LIS_MOC_EWINDFORC(2)     = -9999
   integer :: LVT_LIS_MOC_FHEIGHTFORC(2)   = -9999
   integer :: LVT_LIS_MOC_CHFORC(2)        = -9999
   integer :: LVT_LIS_MOC_CMFORC(2)        = -9999
   integer :: LVT_LIS_MOC_EMISSFORC(2)     = -9999
   integer :: LVT_LIS_MOC_MIXRATIOFORC(2)  = -9999
   integer :: LVT_LIS_MOC_COSZENFORC(2)    = -9999
   integer :: LVT_LIS_MOC_ALBEDOFORC(2)    = -9999

   ! ADDITIONAL Noah3.x variables
   integer :: LVT_LIS_MOC_SOILET(2)  = -9999
   integer :: LVT_LIS_MOC_Z0BRD(2)   = -9999
   integer :: LVT_LIS_MOC_ROUGHNESS(2)   = -9999

   !t2,q2 diagnostics
   integer :: LVT_LIS_MOC_T2DIAG(2) = -9999
   integer :: LVT_LIS_MOC_Q2DIAG(2) = -9999
   integer :: LVT_LIS_MOC_RNET(2) = -9999
   integer :: LVT_LIS_MOC_CH(2)     = -9999
   integer :: LVT_LIS_MOC_CM(2)     = -9999
   integer :: LVT_LIS_MOC_MIXRATIO(2) = -9999

!<for vic>
   !Additional VIC forcing variables
   integer :: LVT_LIS_MOC_SNOWFLAGFORC(2)          = -9999
   integer :: LVT_LIS_MOC_DENSITYFORC(2)           = -9999
   integer :: LVT_LIS_MOC_VAPORPRESSFORC(2)        = -9999
   integer :: LVT_LIS_MOC_VAPORPRESSDEFICITFORC(2) = -9999
   integer :: LVT_LIS_MOC_ARESIST(2)               = -9999
!</for vic>

   ! Sacrament/Snow17 variables
   integer :: LVT_LIS_MOC_SACQS(2)             = -9999
   integer :: LVT_LIS_MOC_SACQSB(2)            = -9999
   integer :: LVT_LIS_MOC_SACSWE(2)            = -9999
   integer :: LVT_LIS_MOC_SACPOTEVAP(2)        = -9999
   integer :: LVT_LIS_MOC_SACEVAP(2)           = -9999
   integer :: LVT_LIS_MOC_SACUZTWC(2)          = -9999
   integer :: LVT_LIS_MOC_SACUZFWC(2)          = -9999
   integer :: LVT_LIS_MOC_SACLZTWC(2)          = -9999
   integer :: LVT_LIS_MOC_SACLZFSC(2)          = -9999
   integer :: LVT_LIS_MOC_SACLZFPC(2)          = -9999
   integer :: LVT_LIS_MOC_SACADIMPC(2)         = -9999
   !integer :: LVT_LIS_MOC_SACSOILMOIST1     = -9999
   integer :: LVT_LIS_MOC_SACSOILMOIST2(2)     = -9999
   !integer :: LVT_LIS_MOC_SACSOILMOIST3     = -9999
   integer :: LVT_LIS_MOC_SACSOILMOIST4(2)     = -9999
   integer :: LVT_LIS_MOC_SACSOILMOIST5(2)     = -9999
   !integer :: LVT_LIS_MOC_SACSOILMOIST6     = -9999
   integer :: LVT_LIS_MOC_SACSOILTEMP1(2)      = -9999
   integer :: LVT_LIS_MOC_SACSOILTEMP2(2)      = -9999
   integer :: LVT_LIS_MOC_SACSFRAC(2)          = -9999
   integer :: LVT_LIS_MOC_SACSDEPTH(2)         = -9999
   integer :: LVT_LIS_MOC_SACUZTWM(2)          = -9999
   integer :: LVT_LIS_MOC_SACUZFWM(2)          = -9999
   integer :: LVT_LIS_MOC_SACLZTWM(2)          = -9999
   integer :: LVT_LIS_MOC_SACLZFSM(2)          = -9999
   integer :: LVT_LIS_MOC_SACLZFPM(2)          = -9999
   integer :: LVT_LIS_MOC_SACADIMP(2)          = -9999
   integer :: LVT_LIS_MOC_SACRAINF(2)          = -9999
   integer :: LVT_LIS_MOC_SACTAIR(2)           = -9999
   integer :: LVT_LIS_MOC_SNOW17SWE(2)         = -9999
   integer :: LVT_LIS_MOC_SNOW17TWE(2)         = -9999
   integer :: LVT_LIS_MOC_SNOW17AEADJ(2)       = -9999
   integer :: LVT_LIS_MOC_SNOW17SWEINCR(2)     = -9999
   integer :: LVT_LIS_MOC_SNOW17WEINCR(2)      = -9999
   integer :: LVT_LIS_MOC_SNOW17TWEINCR(2)     = -9999
   integer :: LVT_LIS_MOC_SNOW17SFRAC1(2)      = -9999
   integer :: LVT_LIS_MOC_SNOW17NEGHS(2)       = -9999
   integer :: LVT_LIS_MOC_SNOW17LIQW(2)        = -9999
   integer :: LVT_LIS_MOC_SNOW17DS(2)          = -9999
   integer :: LVT_LIS_MOC_SNOW17AESCINCR(2)    = -9999
   integer :: LVT_LIS_MOC_SNOW17CNHS(2)        = -9999
   integer :: LVT_LIS_MOC_SNOW17SFRAC4(2)      = -9999
   integer :: LVT_LIS_MOC_SNOW17ACCMAX(2)      = -9999
   integer :: LVT_LIS_MOC_SNOW17SFRAC6(2)      = -9999
   integer :: LVT_LIS_MOC_SNOW17SFRAC7(2)      = -9999
   integer :: LVT_LIS_MOC_SNOW17SFRAC8(2)      = -9999
   integer :: LVT_LIS_MOC_SNOW17SFRAC9(2)      = -9999
   integer :: LVT_LIS_MOC_SNOW17SFRAC10(2)     = -9999
   integer :: LVT_LIS_MOC_SNOW17SFRAC11(2)     = -9999
   integer :: LVT_LIS_MOC_SNOW17SFRAC12(2)     = -9999
   integer :: LVT_LIS_MOC_SNOW17SFRAC13(2)     = -9999
   integer :: LVT_LIS_MOC_SNOW17SFRAC14(2)     = -9999
   integer :: LVT_LIS_MOC_SNOW17SFRAC15(2)     = -9999
   integer :: LVT_LIS_MOC_SNOW17RAINF(2)       = -9999
   integer :: LVT_LIS_MOC_SNOW17SNOWF(2)       = -9999
   integer :: LVT_LIS_MOC_SNOW17TAIR(2)        = -9999
   integer :: LVT_LIS_MOC_SNOW17RMLT(2)        = -9999
   integer :: LVT_LIS_MOC_SNOW17RMINCR(2)      = -9999
   
   integer ::   LVT_LIS_MOC_SACUZTWH(2)  = -9999
   integer ::   LVT_LIS_MOC_SACUZFWH(2)  = -9999
   integer ::   LVT_LIS_MOC_SACLZTWH(2)  = -9999
   integer ::   LVT_LIS_MOC_SACLZFSH(2)  = -9999
   integer ::   LVT_LIS_MOC_SACLZFPH(2)  = -9999
   
   integer ::   LVT_LIS_MOC_SACSWINT(2)  = -9999
   integer ::   LVT_LIS_MOC_SACTSINT(2)  = -9999
   integer ::   LVT_LIS_MOC_SACSWHINT(2)  = -9999
   integer ::   LVT_LIS_MOC_SACFROST(2)  = -9999

!<for vic>
   integer :: LVT_LIS_MOC_VIC_PET_SATSOIL(2)   = -9999
   integer :: LVT_LIS_MOC_VIC_PET_H2OSURF(2)   = -9999
   integer :: LVT_LIS_MOC_VIC_PET_SHORT(2)     = -9999
   integer :: LVT_LIS_MOC_VIC_PET_TALL(2)      = -9999
   integer :: LVT_LIS_MOC_VIC_PET_NATVEG(2)    = -9999
   integer :: LVT_LIS_MOC_VIC_PET_VEGNOCR(2)   = -9999
!</for vic>

   !FLDAS
   integer :: LVT_LIS_MOC_PETFORC(2)         = -9999
   integer :: LVT_LIS_MOC_REFETFORC(2)       = -9999

   ! FLDAS-WRSI OUTPUTS LIST
   integer :: LVT_LIS_MOC_SOS(2) = -9999
   integer :: LVT_LIS_MOC_WRSI(2) = -9999
   integer :: LVT_LIS_MOC_KF2(2) = -9999
   integer :: LVT_LIS_MOC_SumWR(2) = -9999
   integer :: LVT_LIS_MOC_SumET(2) = -9999
   integer :: LVT_LIS_MOC_SWI(2) = -9999
   integer :: LVT_LIS_MOC_SOSa(2) = -9999
   integer :: LVT_LIS_MOC_TotalSurplusWater(2) = -9999
   integer :: LVT_LIS_MOC_MaxSurplusWater(2) = -9999
   integer :: LVT_LIS_MOC_TotalWaterDeficit(2) = -9999
   integer :: LVT_LIS_MOC_MaxWaterDeficit(2) = -9999
   integer :: LVT_LIS_MOC_TotalAETInitial(2) = -9999
   integer :: LVT_LIS_MOC_TotalWRInitial(2) = -9999
   integer :: LVT_LIS_MOC_TotalSurplusWaterInitial(2) = -9999
   integer :: LVT_LIS_MOC_TotalWaterDeficitInitial(2) = -9999
   integer :: LVT_LIS_MOC_TotalAETVeg(2) = -9999
   integer :: LVT_LIS_MOC_TotalWRVeg(2) = -9999
   integer :: LVT_LIS_MOC_TotalSurplusWaterVeg(2) = -9999
   integer :: LVT_LIS_MOC_TotalWaterDeficitVeg(2) = -9999
   integer :: LVT_LIS_MOC_TotalAETFlower(2) = -9999
   integer :: LVT_LIS_MOC_TotalWRFlower(2) = -9999
   integer :: LVT_LIS_MOC_TotalSurplusWaterFlower(2) = -9999
   integer :: LVT_LIS_MOC_TotalWaterDeficitFlower(2) = -9999
   integer :: LVT_LIS_MOC_TotalAETRipe(2) = -9999
   integer :: LVT_LIS_MOC_TotalWRRipe(2) = -9999
   integer :: LVT_LIS_MOC_TotalSurplusWaterRipe(2) = -9999
   integer :: LVT_LIS_MOC_TotalWaterDeficitRipe(2) = -9999
   integer :: LVT_LIS_MOC_PermWiltDate(2) = -9999
   integer :: LVT_LIS_MOC_Wilting1(2) = -9999
   integer :: LVT_LIS_MOC_Wilting2(2) = -9999
   integer :: LVT_LIS_MOC_WRSIa(2) = -9999
   integer :: LVT_LIS_MOC_growing_season(2) = -9999
   integer :: LVT_LIS_MOC_WHC(2)                      = -9999
   integer :: LVT_LIS_MOC_LGP(2)                      = -9999
   integer :: LVT_LIS_MOC_WR_TimeStep(2)              = -9999 ! SY
   integer :: LVT_LIS_MOC_AET_TimeStep(2)             = -9999 ! SY
   integer :: LVT_LIS_MOC_WRSI_TimeStep(2)            = -9999 ! SY
   integer :: LVT_LIS_MOC_SurplusWater_TimeStep(2)    = -9999 ! SY

   !NLDAS
   integer :: LVT_LIS_MOC_CAPEFORC(2)         = -9999

   integer :: LVT_LIS_MOC_LAKE_T_SNOW(2)    =   -9999
   integer :: LVT_LIS_MOC_LAKE_T_ICE(2) =   -9999
   integer :: LVT_LIS_MOC_LAKE_T_MNW(2) =   -9999
   integer :: LVT_LIS_MOC_LAKE_T_WML(2) =   -9999
   integer :: LVT_LIS_MOC_LAKE_T_BOT(2) =   -9999
   integer :: LVT_LIS_MOC_LAKE_T_B1(2)  =   -9999
   integer :: LVT_LIS_MOC_LAKE_C_T(2)   =   -9999
   integer :: LVT_LIS_MOC_LAKE_H_SNOW(2)    =   -9999
   integer :: LVT_LIS_MOC_LAKE_H_ICE(2) =   -9999
   integer :: LVT_LIS_MOC_LAKE_H_ML(2)  =   -9999
   integer :: LVT_LIS_MOC_LAKE_H_B1(2)  =   -9999
   integer :: LVT_LIS_MOC_LAKE_T_SFC(2) =   -9999
   integer :: LVT_LIS_MOC_LAKE_ALBEDO_WATER(2)  =   -9999
   integer :: LVT_LIS_MOC_LAKE_ALBEDO_ICE(2)    =   -9999
   integer :: LVT_LIS_MOC_LAKE_ALBEDO_SNOW(2)   =   -9999
   integer :: LVT_LIS_MOC_LAKE_UFR_A(2) =   -9999
   integer :: LVT_LIS_MOC_LAKE_UFR_W(2) =   -9999
   integer :: LVT_LIS_MOC_LAKE_WCONV(2) =   -9999
   integer :: LVT_LIS_MOC_LAKE_Q_SE(2)  =   -9999
   integer :: LVT_LIS_MOC_LAKE_Q_LA(2)  =   -9999
   integer :: LVT_LIS_MOC_LAKE_I_W(2)   =   -9999
   integer :: LVT_LIS_MOC_LAKE_Q_LWA(2) =   -9999
   integer :: LVT_LIS_MOC_LAKE_Q_LWW(2) =   -9999
   integer :: LVT_LIS_MOC_LAKE_Q_BOT(2) =   -9999


   integer :: LVT_LIS_MOC_EBAL(2)                     = -9999 
   integer :: LVT_LIS_MOC_WBAL(2)                     = -9999 
   integer :: LVT_LIS_MOC_EVAPBAL(2)                  = -9999 
   integer :: LVT_LIS_MOC_SWEOVERP(2)                 = -9999
   integer :: LVT_LIS_MOC_ETOVERP(2)                  = -9999
   integer :: LVT_LIS_MOC_QSOVERP(2)                  = -9999
   integer :: LVT_LIS_MOC_QSBOVERP(2)                 = -9999
   integer :: LVT_LIS_MOC_RUNOFF(2)                   = -9999

   integer :: LVT_LIS_MOC_STREAMFLOW(2) = -9999
   integer :: LVT_LIS_MOC_RIVSTO(2) = -9999
   integer :: LVT_LIS_MOC_RIVDPH(2) = -9999
   integer :: LVT_LIS_MOC_RIVVEL(2) = -9999
   integer :: LVT_LIS_MOC_FLDOUT(2) = -9999
   integer :: LVT_LIS_MOC_FLDEVAP(2) = -9999
   integer :: LVT_LIS_MOC_FLDSTO(2) = -9999
   integer :: LVT_LIS_MOC_FLDDPH(2) = -9999
   integer :: LVT_LIS_MOC_FLDVEL(2) = -9999
   integer :: LVT_LIS_MOC_FLDFRC(2) = -9999
   integer :: LVT_LIS_MOC_FLDARE(2) = -9999
   integer :: LVT_LIS_MOC_SFCELV(2) = -9999
   integer :: LVT_LIS_MOC_RNFSTO(2) = -9999
   integer :: LVT_LIS_MOC_BSFSTO(2) = -9999

   !variables related to RTMs
   integer :: LVT_LIS_MOC_RTM_EMISSIVITY(2)           = -9999 
   integer :: LVT_LIS_MOC_RTM_TB(2)                   = -9999 
   integer :: LVT_LIS_MOC_RTM_SM(2)                   = -9999 

   integer :: LVT_LIS_MOC_IRRIGATEDWATER(2)           = -9999 

   integer :: LVT_LIS_MOC_LSM_COUNT(2)
   integer :: LVT_LIS_MOC_ROUTING_COUNT(2)
   integer :: LVT_LIS_MOC_RTM_COUNT(2)
   integer :: LVT_LIS_MOC_IRRIG_COUNT(2)

#if 0 
   ! SPECIAL CASE INDICES
   ! These are required because Min/Max support cannot be generically
   ! handled for GRIB output.  The routine writeSingleGrib1Var maps
   ! these two entries to LVT_LIS_MOC_TAIRFORC.
   ! They should not be counted in the LVT_LIS_MOC_COUNT total count.
   integer :: LVT_LIS_MOC_TAIRFORC_MIN = 142
   integer :: LVT_LIS_MOC_TAIRFORC_MAX = 143

   ! READ ABOVE NOTE ABOUT SPECIAL CASE INDICES
   integer :: LVT_LIS_MOC_COUNT      = 148
   ! Add the special cases.  LVT_LIS_MOC_GRIB_COUNT should be used only in
   ! LVT_gribMod.F90.
   integer :: LVT_LIS_MOC_GRIB_COUNT = 148 
#endif
   real :: LVT_LIS_MOC_MAX_NUM =  999999.0
   real :: LVT_LIS_MOC_MIN_NUM = -999999.0

!EOP

  type, public :: LVT_LISmetadataEntry
     integer                   :: index
     character*100             :: long_name
     character*100             :: standard_name
     character*100             :: short_name
     character*20              :: units  !unit of the the variable
     integer                   :: nunits !number of units supported
     character*20, allocatable :: unittypes(:) !supported unit types
     integer                   :: ndirs
     character*20              :: dir
     character*20, allocatable :: dirtypes(:)
     character*1               :: format          ! (scientific - E, else - F)
     real,         allocatable :: valid_min(:)
     real,         allocatable :: valid_max(:)
     integer                   :: vlevels
     integer                   :: varid_def
     integer                   :: gribSF          ! grib scale factor
     integer                   :: gribCat
     integer                   :: timeAvgOpt
     integer                   :: selectOpt
     integer                   :: selectStats      !whether to output stats
     integer                   :: computeVar  !variable not present in the LIS 
                                          !output, needs to be computed
     logical                   :: stdev_flag !whether stdev of the measurement is specified
     real,        allocatable  :: stdev(:,:)  !standard deviation of the measurement
     integer,     allocatable  :: count_stdev(:,:)

     integer                   :: minMaxOpt
     integer                   :: stdevOpt
     integer, allocatable      :: count(:,:)
     real, allocatable         :: value(:,:) 

     type(LVT_LISmetadataEntry), pointer :: next
  end type LVT_LISmetadataEntry

  ! To create an array of allocatables, you must create a derived
  ! datatype containing the allocatable.  Then you create
  ! and array of this datatype.
  type, public :: lisdep
     type(LVT_LISmetadataEntry), pointer :: dataEntryPtr
  end type lisdep

  type, public :: lis_output_meta
     integer                 :: xtimeID
     type(LVT_LISmetadataEntry) :: xlat
     type(LVT_LISmetadataEntry) :: xlon

     type(LVT_LISmetadataEntry), pointer :: head_lsm_list
     type(LVT_LISmetadataEntry), pointer :: head_routing_list
     type(LVT_LISmetadataEntry), pointer :: head_rtm_list
     type(LVT_LISmetadataEntry), pointer :: head_irrig_list

     type(lisdep), allocatable, dimension(:) :: ptr_into_lsm_list
     type(lisdep), allocatable, dimension(:) :: ptr_into_routing_list
     type(lisdep), allocatable, dimension(:) :: ptr_into_rtm_list
     type(lisdep), allocatable, dimension(:) :: ptr_into_irrig_list

  end type lis_output_meta

  type(lis_output_meta)     :: LVT_LISoutput(2)

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the processing of the LIS model output
!  as a special LVT datastream
!  
! !FILES USED:
!
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP

!BOP
! 
! !ROUTINE: LVT_readvar_gridded
! \label{LVT_readvar_gridded}
!
! !INTERFACE:
  interface LVT_readvar_gridded
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! This interface provides routines for reading variables (real)
! from a gridded (output) file into a 1d or 2d local gridded space. 
! The gridded files are read by the master processor. 
! The domain decomposition of the "global" grid space on the master processor
! to individual processors are also performed by this routine. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP 
! !PRIVATE MEMBER FUNCTIONS:

     module procedure readvar_2dgridded_real
     module procedure readvar_1dgridded_fromvector_real
!EOP
  end interface

  contains

!BOP
! 
! !ROUTINE: LVT_LISoutputInit
! \label{LVT_LISoutputInit}
!
! !INTERFACE: 
  subroutine LVT_LISoutputInit()
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 

! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine initializes the required arrays to hold the selected list of 
!   LSM variables
!
!   The routines invoked are: 
!   \begin{description}
!    \item[get\_moc\_attributes] (\ref{get_moc_attributes}) \newline
!      reads the LIS model output attributes file entries
!    \item[register\_dataEntry] (\ref{register_dataEntry}) \newline
!      allocates memory and initializes the data structures for
!      specified variables in the LIS model output attributes file
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    type(ESMF_Config) :: modelSpecConfig
    integer           :: k,kk,source
    integer           :: ftn
    integer           :: c,r
    integer           :: nsize
    integer           :: rc
    type(LVT_metadataEntry),    pointer :: dataEntry
    type(LVT_LISmetadataEntry), pointer :: lisdataEntry
    logical                             :: strat_var_found
    logical                             :: var_found
    
    if(LVT_rc%lis_output_obs) then 
       source = 1
       if((LVT_rc%obssource(1).eq."LIS output").and.&
            (LVT_rc%obssource(2).eq."LIS output")) then 
          source = 2
       endif
       
       do kk=1,source
    
          if(LVT_LIS_rc(kk)%wopt.eq."1d tilespace") then 
             nsize = LVT_LIS_rc(kk)%ntiles
          else
             nsize = LVT_rc%ngrid
          endif
          
          LVT_LISoutput(kk)%head_lsm_list     => null()
          LVT_LISoutput(kk)%head_routing_list => null()
          LVT_LISoutput(kk)%head_rtm_list     => null()
          LVT_LISoutput(kk)%head_irrig_list     => null()
          
          LVT_LIS_MOC_LSM_COUNT(kk)     = 0
          LVT_LIS_MOC_ROUTING_COUNT(kk) = 0
          LVT_LIS_MOC_IRRIG_COUNT(kk)   = 0 
          LVT_LIS_MOC_RTM_COUNT(kk)     = 0
          
          modelSpecConfig = ESMF_ConfigCreate(rc=rc)
          call LVT_verify(rc,'config create in readconfig ')
          call ESMF_ConfigLoadFile(modelSpecConfig,&
               trim(LVT_LIS_rc(kk)%outputSpecFile),rc=rc)     
          if(rc.ne.0) then
             write(LVT_logunit,*) '[ERR] '
             write(LVT_logunit,*) '[ERR] Loading LIS output attribute file failed'
             write(LVT_logunit,*) '[ERR] '
             call LVT_endrun()
          endif
    !-------------------------------------------------------------------------
    ! read the meta data attributes for each variable
    !-------------------------------------------------------------------------

          
          call ESMF_ConfigFindLabel(modelSpecConfig,"Swnet:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list, &
               "Swnet",&
               "surface_net_downward_shortwave_flux",&
               "net downward shortwave radiation","F",rc)
          
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk), LVT_LIS_MOC_SWNET(kk), &
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,&
                  nsize,(/"W/m2"/),2,(/"UP","DN"/),valid_min=(/0.0/), &
                  valid_max=(/1200.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lwnet:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Lwnet",&
               "surface_net_downward_longwave_flux",&
               "net downward longwave radiation","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LWNET(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,nsize,&
                  (/"W/m2"/),2,(/"UP","DN"/),valid_min =(/-500.0/), &
                  valid_max=(/510.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Rnet:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
                "Rnet",&
               "net_radiation_flux",&
               "total net radiation","F",rc)

          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_RNET(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,nsize,&
                  (/"W/m2"/),2,(/"UP","DN"/),valid_min =(/-9999.0/), &
                  valid_max=(/-9999.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Qle:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
                "Qle",&
               "surface_upward_latent_heat_flux",&
               "latent heat flux","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_QLE(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"W/m2"/),2,(/"UP","DN"/),valid_min =(/-700.0/),&
                  valid_max=(/700.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Qh:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
                "Qh",&
               "surface_upward_sensible_heat_flux",&
               "sensible heat flux","F",rc)

          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_QH(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,nsize,&
                  (/"W/m2"/),2,(/"UP","DN"/),valid_min=(/-600.0/),&
                  valid_max=(/600.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Qg:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
                "Qg",&
               "downward_heat_flux_in_soil",&
               "soil heat flux","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk), LVT_LIS_MOC_QG(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"W/m2"/), &
                  2,(/"UP","DN"/),valid_min=(/-500.0/),&
                  valid_max=(/500.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Qf:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
                "Qf",&
               "energy_of_fusion",&
               "energy of fusion","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_QF(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,nsize,&
                  (/"W/m2"/),2,(/"S2L","L2S"/),valid_min=(/-1200.0/),&
                  valid_max=(/1200.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Qv:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
                "Qv",&
               "surface_snow_sublimation_heat_flux",&
               "energy of sublimation","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_QV(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,nsize,&
                  (/"W/m2"/),2,(/"S2V","V2S"/),valid_min=(/-600.0/),&
                  valid_max=(/600.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Qtau:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
                "Qtau",&
               "momentum_flux",&
               "momentum flux","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_QTAU(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,nsize,&
                  (/"W/m2"/),2,(/"UP","DN"/),valid_min=(/-100.0/),&
                  valid_max=(/100.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Qa:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
                "Qa",&
               "advective_energy",&
               "advective energy","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_QA(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,nsize,&
                  (/"W/m2"/),2,(/"UP","DN"/),valid_min=(/-50.0/),&
                  valid_max=(/50.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"DelSurfHeat:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "DelSurfHeat",&
               "change_in_heat_storage",&
               "change in heat storage","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_DELSURFHEAT(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"J/m2"/),2,(/"INC","DEC"/),valid_min=(/-500.0/),&
                  valid_max=(/500.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"DelColdCont:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "DelColdCont",&
               "change_in_cold_content",&
               "change in cold content","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_DELCOLDCONT(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"J/m2"/),2,(/"INC","DEC"/),valid_min=(/-600.0/),&
                  valid_max=(/1000.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"BR:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
                "BR",&
               "bowen_ratio",&
               "bowen ratio","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_BR(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,nsize,(/"-"/),&
                  1,("-"),&
                  valid_min=(/0.0/),valid_max=(/1.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"EF:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
                "EF",&
               "evaporative_fraction",&
               "evaporative fraction","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_EF(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),&
                  1,("-"),&
                  valid_min=(/0.0/),valid_max=(/1.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Snowf:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Snowf","snowfall_rate","snowfall rate","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SNOWF(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,2,&
                  nsize,(/"kg/m2s","kg/m2 "/),2,(/"UP","DN"/),&
                  valid_min=(/0.0,0.0/),valid_max=(/0.0085,750.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Rainf:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Rainf","rainfall_rate", "rainfall rate","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_RAINF(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,2,nsize,&
                  (/"kg/m2s","kg/m2 "/),2,(/"UP","DN"/),&
                  valid_min=(/0.0,0.0/),valid_max=(/0.02,2000.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"RainfConv:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "RainfConv","convective_rainfall_rate","convective rainfall","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_RAINFCONV(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2s", "kg/m2 "/),2,(/"UP","DN"/),&
                  valid_min=(/0.0,-9999.0/),valid_max=(/0.02,-9999.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Evap:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Evap","total_evapotranspiration",&
               "total evapotranspiration","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_EVAP(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,3,nsize,&
                  (/"kg/m2s","kg/m2 ","mm/hr "/),2,(/"UP","DN"/),&
                  valid_min=(/-0.0003/),valid_max=(/0.0003/))     
          endif


          call ESMF_ConfigFindLabel(modelSpecConfig,"Qs:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Qs","surface_runoff_amount",&
               "surface runoff","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_QS(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,2,nsize,&
                  (/"kg/m2s", "kg/m2 "/),2,(/"IN ","OUT"/),&
                  valid_min=(/0.0,0.0/),valid_max=(/5.0,43200.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Qsb:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Qsb","subsurface_runoff_amount",&
               "subsurface runoff amount","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_QSB(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,2,nsize,&
                  (/"kg/m2s","kg/m2 "/),2,(/"IN ","OUT"/),&
                  valid_min=(/0.0,0.0/),valid_max=(/5.0,43200.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Qrec:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Qrec","recharge","recharge","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_QREC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,nsize,&
                  (/"kg/m2s"/),2,(/"IN ","OUT"/),&
                  valid_min=(/0.0/),valid_max=(/5.0/))     
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Qsm:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Qsm","snowmelt","snowmelt","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_QSM(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,2,nsize,&
                  (/"kg/m2s", "kg/m2 "/),2,(/"S2L","L2S"/),&
                  valid_min=(/0.0,0.0/),valid_max=(/0.005,450.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Qfz:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Qfz","refreezing_of_water",&
               "refreezing water","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_QFZ(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,nsize,&
                  (/"kg/m2s"/),2,(/"L2S","S2L"/),&
                  valid_min=(/0.0/),valid_max=(/0.005/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Qst:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Qst",&
               "snow_throughfall","snow throughfall","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_QST(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,nsize,&
                  (/"kg/m2s"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/0.005/))     
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"DelSoilMoist:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "DelSoilMoist",&
               "change_in_soil_moisture",&
               "change in soil moisture","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_DELSOILMOIST(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2"/),2,(/"INC","DEC"/),&
                  valid_min=(/-2000.0/),valid_max=(/2000.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"DelSWE:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "DelSWE",&
               "change_in_swe",&
               "change in snow water equivalent","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_DELSWE(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2"/),2,(/"INC","DEC"/),&
                  valid_min=(/-2000.0/),valid_max=(/2000.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"DelSurfStor:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "DelSurfStor",&
               "DelSurfStor",&
               "change in surface water storage","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_DELSURFSTOR(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2"/),2,(/"INC","DEC"/),&
                  valid_min=(/-2000.0/),valid_max=(/2000.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"DelIntercept:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "DelIntercept",&
               "change_in_interception_storage",&
               "change in interception storage","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_DELINTERCEPT(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2"/),2,(/"INC","DEC"/),&
                  valid_min=(/-100.0/),valid_max=(/100.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SnowT:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SnowT",&
               "temperature_in_surface_snow","surface snow temperature","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SNOWT(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,nsize,&
                  (/"K"/),1,(/"-"/),&
                  valid_min=(/213.0/),valid_max=(/280.0/))     
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"VegT:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "VegT","canopy_temperature","canopy temperature","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_VEGT(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,nsize,&
                  (/"K"/),1,(/"-"/),&
                  valid_min=(/213.0/),valid_max=(/333.0/))     
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"BareSoilT:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "BareSoilT","bare_soil_temperature","bare soil temperature","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_BARESOILT(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"K"/),1,(/"-"/),&
                  valid_min=(/213.0/),valid_max=(/343.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"AvgSurfT:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "AvgSurfT","surface_temperature","surface temperature","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_AVGSURFT(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"K"/),1,(/"-"/),&
                  valid_min=(/213.0/),valid_max=(/333.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"RadT:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "RadT",&
               "surface_radiative_temperature",&
               "surface radiative temperature","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_RADT(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,nsize,&
                  (/"K"/),1,(/"-"/),&
                  valid_min=(/213.0/),valid_max=(/353.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Albedo:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Albedo",&
               "surface_albedo","surface albedo","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_ALBEDO(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,2,&
                  nsize,(/"-","%"/),1,(/"-"/),&
                  valid_min=(/0.0,0.0/),valid_max=(/1.0,100.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SWE:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
                "SWE",&
               "liquid_water_content_of_surface_snow",&
               "snow water equivalent","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SWE(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,2,nsize,&
                  (/"kg/m2","m    "/),1,(/"-"/),&
                  valid_min=(/0.0,0.0/),valid_max=(/2000.0,2.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SnowDepth:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SnowDepth","snow_depth","snow depth","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SNOWDEPTH(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  3,nsize,(/"m ","cm","mm"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/10.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SWEVeg:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SWEVeg",&
               "swe_intercepted_by_vegetation",&
               "swe intercepted by vegetation","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SWEVEG(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"kg/m2"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/500.0/))
          endif
          call ESMF_ConfigFindLabel(modelSpecConfig,"SurfStor:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SurfStor",&
               "surface_water_storage",&
               "surface water storage","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SURFSTOR(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/2000.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SoilMoist:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SoilMoist",&
               "soil_moisture_content","soil moisture content","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SOILMOIST(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2", "m3/m3"/),1,(/"-"/),&
                  valid_min=(/0.0,0.0/),valid_max=(/2000.0,0.5/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SoilTemp:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SoilTemp",&
               "soil_temperature","soil temperature","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SOILTEMP(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"K"/),1,(/"-"/),&
                  valid_min=(/213.0/),valid_max=(/333.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SmLiqFrac:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SmLiqFrac",&
               "liquid_fraction_of_soil_moisture",&
               "average layer fraction of liquid moisture","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SMLIQFRAC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"-    ","m3/m3"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))     
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SmFrozFrac:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SmFrozFrac",&
               "frozen_fraction_of_soil_moisture",&
               "average layer fraction of frozen moisture","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SMFROZFRAC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SoilWet:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SoilWet",&
               "total_soil_wetness",&
               "total soil wetness","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SOILWET(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,2,&
                  nsize,(/"- ","% ","mm"/),1,(/"-"/),&
                  valid_min=(/-0.2/),valid_max=(/1.2/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"PotEvap:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "PotEvap","potential_evapotranspiration",&
               "potential evapotranspiration","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_POTEVAP(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,3,&
                  nsize,(/"kg/m2s","mm/hr ","W/m2  "/),2,(/"UP","DN"/),&
                  valid_min=(/-0.0006,-9999.0/),valid_max=(/0.0006,-9999.0/))     
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"ECanop:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "ECanop",&
               "interception_evaporation",&
               "interception evaporation","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_ECANOP(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,2,&
                  nsize,(/"kg/m2s","W/m2  "/),2,(/"UP","DN"/),&
                  valid_min=(/-0.0003/),valid_max=(/0.0003/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TVeg:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "TVeg",&
               "vegetation_transpiration",&
               "vegetation transpiration","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TVEG(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,3,&
                  nsize,(/"kg/m2s","mm/hr ","W/m2  "/),2,(/"UP","DN"/),&
                  valid_min=(/-0.0006,-9999.0/),valid_max=(/0.0006,-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"ESoil:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "ESoil",&
               "bare_soil_evaporation",&
               "bare soil evaporation","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_ESOIL(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,2,&
                  nsize,(/"kg/m2s","W/m2  "/),2,(/"UP","DN"/),&
                  valid_min=(/-0.0003/),valid_max=(/0.0003/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"EWater:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "EWater",&
               "open_water_evaporation",&
               "open water evaporation","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_EWATER(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"kg/m2s"/),2,(/"UP","DN"/),&
                  valid_min=(/-0.0003/),valid_max=(/0.0003/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"RootMoist:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "RootMoist",&
               "root_zone_soil_moisture",&
               "root zone soil moisture","F",rc)

          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_ROOTMOIST(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"m3/m3", "kg/m2"/),1,(/"-"/),&
                  valid_min=(/0.0,0.0/),valid_max=(/0.5,2000.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"CanopInt:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "CanopInt",&
               "total_canopy_water_storage",&
               "total canopy water storage","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_CANOPINT(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/100.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"EvapSnow:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "EvapSnow",&
               "snow_evaporation",&
               "snow evaporation","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_EVAPSNOW(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2s"/),1,(/"-"/),&
                  valid_min=(/-0.0003/),valid_max=(/0.0003/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SubSnow:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SubSnow",&
               "snow_sublimation",&
               "snow sublimation","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SUBSNOW(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,2,&
                  nsize,(/"kg/m2s","W/m2  "/),1,(/"-"/),&
                  valid_min=(/-0.0003/),valid_max=(/0.0003/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SubSurf:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SubSurf",&
               "sublimation_of_the_snow_free_area",&
               "sublimation of the snow free area","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SUBSURF(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"kg/m2s"/),1,(/"-"/),&
                  valid_min=(/-0.0003/),valid_max=(/0.0003/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"ACond:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "ACond","aerodynamic_conductance",&
               "aerodynamic conductance","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_ACOND(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"m/s"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"CCond:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "CCond","canopy_conductance","canopy conductance","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_CCOND(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"m/s"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"WaterTableD:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "WaterTableD",&
               "water_table_depth",&
               "water table depth","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_WaterTableD(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"m"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TWS:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "TWS",&
               "terrestrial_water_storage",&
               "terrestrial water storage","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TWS(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"mm"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Snowcover:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SnowCover","surface_snow_area_fraction","snow cover","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SNOWCOVER(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SAlbedo:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SAlbedo","snow_albedo","snow albedo","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SALBEDO(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SnowTProf:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SnowTProf","snow_temperature_profile","snow temperature profile","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SNOWTPROF(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"K"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SLiqFrac:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SLiqFrac",&
               "snow_liquid_fraction_on_ground",&
               "snow liquid fraction on ground","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SLIQFRAC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"LWup:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LWup",&
               "longwave_radiation_up",&
               "longwave radiation up","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LWUP(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"W/m2"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"GPP:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "GPP",&
               "gross_primary_production",&
               "gross primary production","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_GPP(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2s2 ","umol/m2s"/),2,(/"UP","DN"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"NPP:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "NPP",&
               "net_primary_productivity",&
               "net primary productivity","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_NPP(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2s2 ","umol/m2s"/),2,(/"UP","DN"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif


          call ESMF_ConfigFindLabel(modelSpecConfig,"NEE:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "NEE",&
               "net_ecosystem_exchange",&
               "net ecosystem exchange","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_NEE(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2s2 ","umol/m2s"/),2,(/"UP","DN"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"AutoResp:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "AutoResp",&
               "autoresp",&
               "","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_AUTORESP(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2s2 ","umol/m2s"/),2,(/"UP","DN"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"HeteroResp:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "HeteroResp",&
               "heteroresp",&
               "","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_HETERORESP(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2s2 ","umol/m2s"/),2,(/"UP","DN"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"LeafResp:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "LeafResp",&
               "leafresp",&
               "","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LEAFRESP(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2s2 ","umol/m2s"/),2,(/"UP","DN"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotSoilCarb:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "total_soil_and_litter_carbon_content",&
               "total soil and litter carbon content",&
               "","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TOTSOILCARB(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotLivBiom:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "total_living_biomass_carbon_content",&
               "total living biomass carbon content",&
               "","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TOTLIVBIOM(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SoilET:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "SoilET",&
               "soil_evaporation",&
               "soil evaporation","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SOILET(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"W/m2"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))     
          endif


          call ESMF_ConfigFindLabel(modelSpecConfig,"Z0brd:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "z0brd",&
               "z0brd",&
               "z0brd","F",rc) 
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_Z0BRD(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"m"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))     
          endif


          call ESMF_ConfigFindLabel(modelSpecConfig,"Ch:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "Ch",&
               "heat_exchange_coefficient",&
               "heat exchange coefficient","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_CH(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))     
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Cm:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "Cm",&
               "momentum_exchange_coefficient",&
               "momentum exchange coefficient","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_CM(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))     
          endif


          call ESMF_ConfigFindLabel(modelSpecConfig,"T2diag:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "T2diag","diagnostic_t2","diagnostic t2","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_T2DIAG(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"K"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))       
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Q2diag:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Q2diag","diagnostic_q2","diagnostic q2","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_Q2DIAG(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/kg"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"RootTemp:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "RootTemp","root_zone_temperature","root zone temperature","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_ROOTTEMP(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"K"/),1,(/"-"/),&
                  valid_min=(/213.0/),valid_max=(/350.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Wind_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Wind_f","wind_speed","wind speed","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_WINDFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"m/s   ","km/day"/),1,(/"-"/),&
                  valid_min=(/-75.0,-6500.0/),valid_max=(/75.0,6500.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Rainf_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Rainf_f","rainfall_flux","rainfall flux","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_RAINFFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2s", "kg/m2 "/),2,(/"UP","DN"/),&
                  valid_min=(/0.0/),valid_max=(/0.02/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Snowf_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Snowf_f","snowfall_flux","snowfall flux","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SNOWFFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2s"/),2,(/"UP","DN"/),&
                  valid_min=(/0.0/),valid_max=(/0.02/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"CRainf_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "CRainf_f","convective_rainfall_flux","convective rainfall flux","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_CRAINFFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2s", "kg/m2 "/),2,(/"UP","DN"/),&
                  valid_min=(/0.0/),valid_max=(/0.02/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Tair_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Tair_f","air_temperature","air temperature","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SNOWFFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"K"/),2,(/"UP","DN"/),&
                  valid_min=(/0.0/),valid_max=(/0.02/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Qair_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Qair_f","specific_humidity","specific humidity","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_QAIRFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/kg"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/0.03/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Psurf_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Psurf_f","surface_air_pressure","surface pressure","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_PSURFFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"Pa"/),1,("-"),&
                  valid_min=(/5000.0/),valid_max=(/110000.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SWdown_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SWdown_f",&
               "surface_downwelling_shortwave_flux_in_air",&
               "downward shortwave radiation","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SWDOWNFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"W/m2"/),2,(/"UP","DN"/),&
                  valid_min=(/0.0/),valid_max=(/1360.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"LWdown_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LWdown_f",&
               "surface_downwelling_longwave_flux_in_air",&
               "downward longwave radiation","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LWDOWNFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"W/m2"/),2,(/"UP","DN"/),&
                  valid_min=(/0.0/),valid_max=(/750.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"DirectSW_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "DirectSW_f",&
               "surface_direct_downwelling_shortwave_flux_in_air",&
               "direct shortwave flux","F", rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_DIRECTSWFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"W/m2"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/1360.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"DiffuseSW_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "DiffuseSW_f",&
               "surface_diffuse_downwelling_shortwave_flux_in_air",&
               "diffuse shortwave flux","F", rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_DIFFUSESWFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"W/m2"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/1360.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"NWind_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "NWind_f",&
               "northward_wind",&
               "northward wind","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_NWINDFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"m/s"/),2,(/"N","E"/),&
                  valid_min=(/-75.0/),valid_max=(/75.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"EWind_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "EWind_f",&
               "eastward_wind",&
               "eastward wind","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_EWINDFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"m/s"/),2,(/"N","E"/),&
                  valid_min=(/-75.0/),valid_max=(/75.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"FHeight_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "FHeight_f","forcing_height","forcing height","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_FHEIGHTFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"m"/),1,("-"),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Ch_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "CH_f","heat_exchange_coefficient",&
               "heat exchange coefficient","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_CHFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"-"/),1,("-"),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Cm_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "CM_f",&
               "momentum_exchange_coefficient",&
               "momentum exchange coefficient","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_CMFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"-"/),1,("-"),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"MixRatio_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "MixRatio_f","mixing_ratio","mixing ratio","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_MIXRATIOFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,("-"),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"CosZenith_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "CosZenith_f",&
               "cosine_solar_zenith_angle","cosine of solar zenith angle",&
               "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_COSZENFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,("-"),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Albedo_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Albedo_f","albedo","albedo","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_ALBEDOFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/1.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"PARDR_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "PARDR_f",&
               "surface_downward_PAR_direct",&
               "surface downward PAR direct","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_PARDRFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"W/m2"/),1,("DN"),&
                  valid_min=(/0.0/),valid_max=(/1.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"PARDF_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "PARDF_f",&
               "surface_downward_PAR_diffuse",&
               "surface downward PAR diffuse","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_PARDFFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"W/m2"/),1,("DN"),&
                  valid_min=(/0.0/),valid_max=(/1.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Landmask:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Landmask","landmask","landmask","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LANDMASK(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Landcover:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Landcover","landcover","landcover","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LANDCOVER(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/24.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Soiltype:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Soiltype","soiltype","soiltype","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SOILTYPE(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/20.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SandFrac:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SandFrac","sandfrac","fraction of sand","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SANDFRAC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"ClayFrac:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "ClayFrac","clayfrac","fraction of clay","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_CLAYFRAC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SiltFrac:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SiltFrac","siltfrac","fraction of silt","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SILTFRAC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Porosity:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Porosity","porosity","porosity","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_POROSITY(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Soilcolor:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SoilColor","soilcolor","soil color","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SOILCOLOR(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Elevation:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Elevation","elevation","elevation","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_ELEVATION(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"m"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Slope:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Slope","slope","slope","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SLOPE(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"LAI:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LAI","leaf_area_index","leaf area index","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAI(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif


          call ESMF_ConfigFindLabel(modelSpecConfig,"SAI:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SAI","stem_area_index","stem area index","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SAI(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,1,&
                  nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Snfralbedo:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Snfralbedo","snow_free_albedo","snow free albedo","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SNFRALBEDO(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Mxsnalbedo:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Mxsnalbedo","maximum_snow_free_albedo","maximum snow free albedo","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_MXSNALBEDO(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Greenness:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Greenness",&
               "green_vegetation_fraction","green vegetation fraction","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_GREENNESS(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"-","%"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Tempbot:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Tempbot","bottom_temperature","bottom temperature","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TEMPBOT(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"K"/),1,(/"-"/),&
                  valid_min=(/213.0/),valid_max=(/333.0/))

          endif

          !fldas 
          call ESMF_ConfigFindLabel(modelSpecConfig,"PET_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "PET_f","potential_evaporation","potential evaporation","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_PETFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2s"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"RefET_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "RefET_f","reference_ET","reference ET","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_REFETFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2s"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"CAPE_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "CAPE_f",&
               "convective_available_potential_energy",&
               "convective available potential energy","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_CAPEFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"J/kg"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif
          call ESMF_ConfigFindLabel(modelSpecConfig,"SOS:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SOS","start_of_season","start of season","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SOS(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"WRSI:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "WRSI","water_requirements_satisfaction_index",&
               "water requirements satisfaction index","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_WRSI(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"KF2:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
                "KF2",&
               "KF2","KF2","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_KF2(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SumWR:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SumWR","SumWR","SumWR","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SumWR(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SumET:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SumET","SumET", "SumET","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SumET(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SWI:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SWI","soil_water_index","soil water index","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SWI(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SOSa:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "SOSa","SOSa", "SOSa","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SOSa(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalSurplusWater:",rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list, &
               "TotalSurplusWater",&
               "total_surplus_water","total surplus water", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TotalSurplusWater(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"MaxSurplusWater:",rc=rc)
          call get_moc_attributes(modelSpecConfig,&
               LVT_LISoutput(kk)%head_lsm_list, &
               "MaxSurplusWater",&
               "max_surplus_water","max surplus water", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_MaxSurplusWater(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWaterDeficit:",rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list, &
               "TotalWaterDeficit",&
               "total_water_deficit", "total water deficit","F", rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TotalWaterDeficit(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif


          call ESMF_ConfigFindLabel(modelSpecConfig,"MaxWaterDeficit:",rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list, &
               "MaxWaterDeficit",&
               "max_water_deficit", "max water deficit","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_MaxWaterDeficit(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalAETInitial:",rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list, &
               "TotalAETInitial",&
               "total_AET_initial", "total AET initial","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TotalAETInitial(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWRInitial:",rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list, &
               "TotalWRInitial",&
               "total_WR_initial","total WR initial","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TotalWRInitial(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalSurplusWaterInitial:",&
               rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list, &
               "TotalSurplusWaterInitial","total_surplus_water_initial",&
               "total surplus water initial", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),&
                  LVT_LIS_MOC_TotalSurplusWaterInitial(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWaterDeficitInitial:",&
               rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list, &
               "TotalWaterDeficitInitial",&
               "total_water_deficit_initial","total water deficit initial","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),&
                  LVT_LIS_MOC_TotalWaterDeficitInitial(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalAETVeg:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "TotalAETVeg","total_AET_veg","total AET veg","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TotalAETVeg(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWRVeg:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "TotalWRVeg","total_WR_veg","total WR veg","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TotalWRVeg(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalSurplusWaterVeg:",rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list, &
               "TotalSurplusWaterVeg",&
               "total_surplus_water_veg","total surplus water veg","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TotalSurplusWaterVeg(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWaterDeficitVeg:",rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list, &
               "TotalWaterDeficitVeg",&
               "total_water_deficit_veg","total water deficit veg","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TotalWaterDeficitVeg(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalAETFlower:",rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list,&
               "TotalAETFlower",&
               "totalAETFlower","totalAETFlower","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TotalAETFlower(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWRFlower:",rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list, &
               "TotalWRFlower", "TotalWRFlower",&
               "TotalWRFlower","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TotalWRFlower(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalSurplusWaterFlower:",rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list, &
               "TotalSurplusWaterFlower","TotalSurplusWaterFlower",&
               "TotalSurplusWaterFlower","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),&
                  LVT_LIS_MOC_TotalSurplusWaterFlower(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWaterDeficitFlower:",rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list, &
               "TotalWaterDeficitFlower","TotalWaterDeficitFlower",&
               "TotalWaterDeficitFlower","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),&
                  LVT_LIS_MOC_TotalWaterDeficitFlower(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalAETRipe:",rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list, &
               "TotalAETRipe",&
               "TotalAETRipe","TotalAETRipe","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TotalAETRipe(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWRRipe:",rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list, &
               "TotalWRRipe","TotalWRRipe",&
               "TotalWRRipe","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TotalWRRipe(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalSurplusWaterRipe:",rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list,&
               "TotalSurplusWaterRipe",&
               "TotalSurplusWaterRipe","TotalSurplusWaterRipe","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TotalSurplusWaterRipe(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWaterDeficitRipe:",rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list, &
               "TotalWaterDeficitRipe","TotalWaterDeficitRipe",&
               "TotalWaterDeficitRipe","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),&
                  LVT_LIS_MOC_TotalWaterDeficitRipe(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"PermWiltDate:",rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list, &
               "PermWiltDate",&
               "PermWiltDate","PermWiltDate","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_PermWiltDate(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Wilting1:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Wilting1","Wilting1","Wilting1","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_Wilting1(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))


          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Wilting2:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Wilting2","Wilting2","Wilting2","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_Wilting2(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"WRSIa:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "WRSIa","WRSIa", "WRSIa","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_WRSIa(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"growing_season:",rc=rc)
          call get_moc_attributes(modelSpecConfig, &
               LVT_LISoutput(kk)%head_lsm_list,&
               "growing_season",&
               "growing_season","growing season","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_growing_season(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"WHC:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "WHC",&
               "WHC",&
               "WHC","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_WHC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))     
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"LGP:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "LGP",&
               "LGP",&
               "LGP","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LGP(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))     
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"WR_TimeStep:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "WR_TimeStep",&
               "WR_TimeStep",&
               "WR_TimeStep","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_WR_TimeStep(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))     
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"AET_TimeStep:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "AET_TimeStep",&
               "AET_TimeStep",&
               "AET_TimeStep","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_AET_TimeStep(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))     
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"WRSI_TimeStep:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "WRSI_TimeStep",&
               "WRSI_TimeStep",&
               "WRSI_TimeStep","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_WRSI_TimeStep(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))     
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"SurplusWater_TimeStep:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "SurplusWater_TimeStep",&
               "SurplusWater_TimeStep",&
               "SurplusWater_TimeStep","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SurplusWater_TimeStep(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))     
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Snowflag_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "Snowflag",&
               "Snowflag",&
               "Snowflag","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SNOWFLAGFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))     
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Density_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "Atmospheric_Density",&
               "Atmospheric_Density",&
               "Atmospheric_Density","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_DENSITYFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m3"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))     
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"VaporPress_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "Vapor_Pressure",&
               "Vapor_Pressure",&
               "Vapor_Pressure","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_VAPORPRESSFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))     
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"VaporPressDeficit_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "Vapor_Pressure_Deficit",&
               "Vapor_Pressure_Deficit",&
               "Vapor_Pressure_Deficit","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_VAPORPRESSDEFICITFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))     
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"AResist:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "Aerodynamic_Resistance",&
               "Aerodynamic_Resistance",&
               "Aerodynamic_Resistance","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_ARESIST(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"s/m"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))     
          endif


          call ESMF_ConfigFindLabel(modelSpecConfig,"sac_tsint:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "sac_tsint",&
               "sac_tsint",&
               "sac soil temperature of intended layer", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SACTSINT(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"K"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"sac_swint:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "sac_swint",&
               "sac_swint",&
               "sac total volumetric soil moisture content of intented layer", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SACSWINT(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"m3/m3"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif
          call ESMF_ConfigFindLabel(modelSpecConfig,"sac_swhint:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "sac_swhint",&
               "sac_swhint",&
               "sac liquid volumetric soil moisture content of intented layer", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SACSWHINT(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"m3/m3"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"sac_frost:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "sac_frost",&
               "sac_frost",&
               "sac_frost", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SACFROST(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"sac_uztwc:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "sac_uztwc",&
               "sac_uztwc",&
               "sac_uztwc", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SACUZTWC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"mm","- "/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif


          call ESMF_ConfigFindLabel(modelSpecConfig,"sac_uzfwc:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "sac_uzfwc",&
               "sac_uzfwc",&
               "sac_uzfwc", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SACUZFWC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"mm","- "/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"sac_lztwc:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "sac_lztwc",&
               "sac_lztwc",&
               "sac_lztwc", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SACLZTWC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"mm","- "/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"sac_lzfsc:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "sac_lzfsc",&
               "sac_lzfsc",&
               "sac_lzfsc", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SACLZFSC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"mm","- "/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"sac_lzfpc:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "sac_lzfpc",&
               "sac_lzfpc",&
               "sac_lzfpc", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SACLZFPC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"mm","- "/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"sac_adimpc:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "sac_adimpc",&
               "sac_adimpc",&
               "sac_adimpc", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SACADIMPC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"mm","- "/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"sac_uztwh:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "sac_uztwh",&
               "sac_uztwh",&
               "sac_uztwh", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SACUZTWH(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"mm"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"sac_uzfwh:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "sac_uzfwh",&
               "sac_uzfwh",&
               "sac_uzfwh", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SACUZFWH(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"mm"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"sac_lztwh:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "sac_lztwh",&
               "sac_lztwh",&
               "sac_lztwh", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SACLZTWH(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"mm"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"sac_lzfsh:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "sac_lzfsh",&
               "sac_lzfsh",&
               "sac_lzfsh", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SACLZFSH(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"mm"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"sac_lzfph:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "sac_lzfph",&
               "sac_lzfph",&
               "sac_lzfph", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SACLZFPH(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"mm"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"snow17_swe:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "snow17_swe",&
               "snow17_swe",&
               "snow17_swe", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SNOW17SWE(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"snow17_aeadj:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "snow17_aeadj",&
               "snow17_aeadj",&
               "snow17_aeadj", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SNOW17AEADJ(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"mm"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"snow17_neghs:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "snow17_neghs",&
               "snow17_neghs",&
               "snow17_neghs", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SNOW17NEGHS(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"mm"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"snow17_liqw:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "snow17_liqw",&
               "snow17_liqw",&
               "snow17_liqw", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SNOW17LIQW(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"kg/m2"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"snow17_accmax:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "snow17_accmax",&
               "snow17_accmax",&
               "snow17_accmax", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SNOW17ACCMAX(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"mm"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"snow17_rmlt:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "snow17_rmlt",&
               "snow17_rmlt",&
               "snow17_rmlt", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SNOW17RMLT(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2 ","kg/m2s"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif
          call ESMF_ConfigFindLabel(modelSpecConfig,"vic_pet_satsoil:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "vic_pet_satsoil",&
               "vic_pet_satsoil",&
               "vic_pet_satsoil", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_VIC_PET_SATSOIL(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2 ","kg/m2s"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"vic_pet_h2osurf:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "vic_pet_h2osurf",&
               "vic_pet_h2osurf",&
               "vic_pet_h2osurf", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_VIC_PET_H2OSURF(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2 ","kg/m2s"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"vic_pet_short:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "vic_pet_short",&
               "vic_pet_short",&
               "vic_pet_short", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_VIC_PET_SHORT(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2 ","kg/m2s"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"vic_pet_tall:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "vic_pet_tall",&
               "vic_pet_tall",&
               "vic_pet_tall", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_VIC_PET_TALL(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2 ","kg/m2s"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"vic_pet_natveg:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "vic_pet_natveg",&
               "vic_pet_natveg",&
               "vic_pet_natveg", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_VIC_PET_NATVEG(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2 ","kg/m2s"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"vic_pet_vegnocr:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "vic_pet_vegnocr",&
               "vic_pet_vegnocr",&
               "vic_pet_vegnocr", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_VIC_PET_VEGNOCR(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2 ","kg/m2s"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Tsnow:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeTsnow",&
               "Lake_temperature_at_air_snow_interface",&
               "Lake temperature at the air snow interface", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_T_SNOW(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"K"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Tice:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeTice",&
               "Lake_temperature_at_snow_ice_interface",&
               "Lake temperature at the snow ice interface", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_T_ICE(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"K"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Tmnw:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeTmnw",&
               "Mean_temperature_of_the_water_column",&
               "Mean temperature of the water column", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_T_MNW(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"K"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Twml:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeTwml",&
               "Lake_temperature_of_the_mixed_layer",&
               "Lake temperature of the mixed layer", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_T_WML(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"K"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Tbot:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeTbot",&
               "Lake_temperature_at_the_water_bottom",&
               "Lake temperature at the water bottom", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_T_BOT(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"K"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Tb1:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeTb1",&
               "Temperature_at_the_bottom_of_upper_layer_of_sediments",&
               "Temperature at the bottom of upper layer of sediments", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_T_B1(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"K"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_CT:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeCT",&
               "Thermocline_shape_factor_of_lake",&
               "Thermocline shape factor of lake", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_C_T(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"1"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Hice:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeHice",&
               "Ice_thickness_above_lake",&
               "Ice thickness above lake", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_H_ICE(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"m"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Hml:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeHml",&
               "Thickness_of_mixed_layer_of_lake",&
               "Thickness of mixed layer of lake", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_H_ML(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"m"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Hb1:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeHb1",&
               "Thickness_of_upper_layer_of_bottom_sediments",&
               "Thickness of upper layer of bottom sediments", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_H_B1(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"m"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Walbedo:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeWalbedo",&
               "Water_surface_albedo_over_lake",&
               "Water surface albedo over lake", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_ALBEDO_WATER(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif


          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_IceAlbedo:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeIceAlbedo",&
               "Ice_surface_albedo_over_lake",&
               "Ice_surface_albedo_over_lake", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_ALBEDO_ICE(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"m"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_SnowAlbedo:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeSnowAlbedo",&
               "Snow_surface_albedo_over_lake",&
               "Snow surface albedo over lake", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_ALBEDO_SNOW(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_UFRa:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeUFRa",&
               "Lake_friction_velocity_in_air",&
               "Lake friction velocity in air", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_UFR_A(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"m/s"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_UFRw:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeUFRw",&
               "Lake_friction_velocity_in_surface_water",&
               "Lake friction velocity in surface water", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_UFR_W(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"m/s"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_WConv:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeWConv",&
               "Lake_convective_velocity_scale",&
               "Lake convective velocity scale", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_WCONV(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"m/s"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_IW:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeIW",&
               "Lake_radiation_flux_at_the_interface",&
               "Lake radiation flux at the interface", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_I_W(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"W/m2"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Qbot:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "LakeQbot",&
               "Lake_heat_flux_across_water_sediment_boundary",&
               "Lake heat flux across water sediment boundary", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_LAKE_Q_BOT(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"W/m2"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"RelSMC:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "RelSMC",&
               "relative_soil_moisture",&
               "relative soil moisture", "F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_RELSMC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"%    ","m3/m3"/),1,(/"-"/),&
                  valid_min=(/0.0/),valid_max=(/1.0/))    
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"RHMin:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "RHMin","min_relative_humidity","minimum relative humidity","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_RHMIN(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,2,&
                  nsize,(/"-","%"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))
          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"TotalPrecip:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "TotalPrecip","total_precipitation_amount","precipitation amount","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_TOTALPRECIP(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  2,nsize,(/"kg/m2s", "kg/m2 "/),2,(/"UP","DN"/),&
                  valid_min=(/0.0,-9999.0/),valid_max=(/0.02,-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Emiss_f:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Emiss_f","emissivity","emissivity","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_EMISSFORC(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,("-"),&
                  valid_min=(/0.7/),valid_max=(/1.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig,"Roughness:",rc=rc)
          call get_moc_attributes(modelSpecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "Roughness",&
               "surface_roughness_length",&
               "surface roughness length","F",rc)  
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_ROUGHNESS(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"m"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))     
          endif


          call ESMF_ConfigFindLabel(modelSpecConfig,"Runoff:",rc=rc)
          call get_moc_attributes(modelSpecConfig, LVT_LISoutput(kk)%head_lsm_list,&
               "Runoff","Total_runoff_amount",&
               "Total runoff amount","E",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_RUNOFF(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,2,nsize,&
                  (/"kg/m2s","kg/m2 "/),2,(/"IN ","OUT"/),&
                  valid_min=(/0.0,0.0/),valid_max=(/5.0,43200.0/))

          endif


          call ESMF_ConfigFindLabel(modelSpecConfig, "EBAL:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "EBAL","energy_balance","energy balane","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_EBAL(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/5.0/),valid_max=(/5.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "WBAL:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "WBAL","water_balance","water balance","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_WBAL(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "EVAPBAL:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "EVAPBAL","evaporation_balance","evaporation balance","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_EVAPBAL(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "SWE/P:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "SWE/P","SWE_over_P","SWE normalized by precipitation","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_SWEOVERP(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "ET/P:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "ET/P","ET_over_P","Evapotranspiration efficiency","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_ETOVERP(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "Qs/P:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "Qs/P","surface_runoff_ratio","Surface runoff ratio","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_QSOVERP(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "Qsb/P:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_lsm_list,&
               "Qsb/P","subsurface_runoff_ratio","Subsurface_runoff_ratio","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_LSM_COUNT(kk),LVT_LIS_MOC_QSBOVERP(kk),&
                  LVT_LISoutput(kk)%head_lsm_list,&
                  1,nsize,(/"-"/),1,(/"-"/),&
                  valid_min=(/-9999.0/),valid_max=(/-9999.0/))

          endif

          !routing outputs
          call ESMF_ConfigFindLabel(modelSpecConfig, "Streamflow:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_routing_list,&
               "Streamflow","streamflow","streamflow","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_ROUTING_COUNT(kk), LVT_LIS_MOC_STREAMFLOW(kk), &
                  LVT_LISoutput(kk)%head_routing_list, &
                  1,nsize,(/"m3/s"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/500000.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "RiverStor:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_routing_list,&
               "RiverStor",&
               "River_Water_Storage",&
               "River Water Storage","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_ROUTING_COUNT(kk), LVT_LIS_MOC_RIVSTO(kk), &
                  LVT_LISoutput(kk)%head_routing_list, &
                  1,nsize,(/"m3"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/500000.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "RiverDepth:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_routing_list,&
               "RiverDepth",&
               "River_Depth",&
               "River Depth","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_ROUTING_COUNT(kk), LVT_LIS_MOC_RIVDPH(kk), &
                  LVT_LISoutput(kk)%head_routing_list, &
                  1,nsize,(/"m"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/500000.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "RiverVelocity:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_routing_list,&
               "RiverFlowVelocity",&
               "River_Flow_Velocity",&
               "River Flow Velocity","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_ROUTING_COUNT(kk), LVT_LIS_MOC_RIVVEL(kk), &
                  LVT_LISoutput(kk)%head_routing_list, &
                  1,nsize,(/"m/s"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/500000.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "FloodQ:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_routing_list,&
               "FloodQ",&
               "Floodplain_Water_Discharge",&
               "Floodplain Water Discharge","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_ROUTING_COUNT(kk), LVT_LIS_MOC_fldout(kk), &
                  LVT_LISoutput(kk)%head_routing_list, &
                  1,nsize,(/"m3/s"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/500000.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "FloodEvap:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_routing_list,&
               "FloodEvap",&
               "Floodplain_Evaporation",&
               "Floodplain Evaporation","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_ROUTING_COUNT(kk), LVT_LIS_MOC_fldevap(kk), &
                  LVT_LISoutput(kk)%head_routing_list, &
                  1,nsize,(/"m3"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/500000.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "FloodStor:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_routing_list,&
               "FloodStor",&
               "Floodplain_Water_Storage",&
               "Floodplain Water Storage","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_ROUTING_COUNT(kk), LVT_LIS_MOC_FLDSTO(kk), &
                  LVT_LISoutput(kk)%head_routing_list, &
                  1,nsize,(/"m3"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/500000.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "FloodDepth:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_routing_list,&
               "FloodDepth",&
               "Floodplain_Depth",&
               "Floodplain Depth","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_ROUTING_COUNT(kk), LVT_LIS_MOC_FLDDPH(kk), &
                  LVT_LISoutput(kk)%head_routing_list, &
                  1,nsize,(/"m"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/500000.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "FloodVelocity:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_routing_list,&
               "FloodVelocity",&
               "Floodplain Flow Velocity",&
               "Floodplain_Flow_Velocity","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_ROUTING_COUNT(kk), LVT_LIS_MOC_FLDVEL(kk), &
                  LVT_LISoutput(kk)%head_routing_list, &
                  1,nsize,(/"m/s"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/500000.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "FloodedFrac:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_routing_list,&
               "FloodedFrac",&
               "Flooded_Fraction",&
               "Flooded Fraction","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_ROUTING_COUNT(kk), LVT_LIS_MOC_FLDFRC(kk), &
                  LVT_LISoutput(kk)%head_routing_list, &
                  1,nsize,(/"-"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/500000.0/))

          endif
          call ESMF_ConfigFindLabel(modelSpecConfig, "FloodedArea:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_routing_list,&
               "FloodedArea",&
               "Flooded Area",&
               "Flooded_Area","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_ROUTING_COUNT(kk), LVT_LIS_MOC_FLDARE(kk), &
                  LVT_LISoutput(kk)%head_routing_list, &
                  1,nsize,(/"m2"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/500000.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "SurfElev:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_routing_list,&
               "SurfElev",&
               "Surface Water Elevation",&
               "Surface_Water_Elevation","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_ROUTING_COUNT(kk), LVT_LIS_MOC_SFCELV(kk), &
                  LVT_LISoutput(kk)%head_routing_list, &
                  1,nsize,(/"m"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/500000.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "RunoffStor:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_routing_list,&
               "RunoffStor",&
               "Runoff Reservoir Storage",&
               "Runoff_Reservoir_Storage","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_ROUTING_COUNT(kk), LVT_LIS_MOC_RNFSTO(kk), &
                  LVT_LISoutput(kk)%head_routing_list, &
                  1,nsize,(/"m3"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/500000.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "BaseflowStor:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_routing_list,&
               "BaseflowStor",&
               "Baseflow Reservoir Storage",&
               "Baseflow_Reservoir_Storage","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_ROUTING_COUNT(kk), LVT_LIS_MOC_BSFSTO(kk), &
                  LVT_LISoutput(kk)%head_routing_list, &
                  1,nsize,(/"m3"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/500000.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "RTM emissivity:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_rtm_list,&
               "RTM_emissivity","rtm_emissivity","rtm emissivity","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_RTM_COUNT(kk),LVT_LIS_MOC_RTM_EMISSIVITY(kk), &
                  LVT_LISoutput(kk)%head_rtm_list, &
                  1,nsize,(/"-"/),1,("-"),&
                  valid_min=(/0.7/),valid_max=(/1.0/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "RTM Tb:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_rtm_list,&
               "RTM_Tb","rtm_brightness_temperature",&
               "rtm brightness temperature","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_RTM_COUNT(kk),LVT_LIS_MOC_RTM_TB(kk), &
                  LVT_LISoutput(kk)%head_rtm_list, &
                  1,nsize,(/"K"/),1,("-"),&
                  valid_min=(/213.0/),valid_max=(/353.0/))

          endif


          call ESMF_ConfigFindLabel(modelSpecConfig, "RTM SoilMoist:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_rtm_list,&
               "RTM_SoilMoist","RTM_SoilMoist",&
               "RTM SoilMoist","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_RTM_COUNT(kk),LVT_LIS_MOC_RTM_SM(kk), &
                  LVT_LISoutput(kk)%head_rtm_list, &
                  1,nsize,(/"m3/m3"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/0.50/))

          endif

          call ESMF_ConfigFindLabel(modelSpecConfig, "Irrigated water:",rc=rc)
          call get_moc_attributes(modelspecConfig,LVT_LISoutput(kk)%head_irrig_list,&
               "IrrigatedWater","total_irrigated_water_amount",&
               "irrigated water amount","F",rc)
          if(rc.eq.1) then 
             call register_dataEntry(LVT_LIS_MOC_IRRIG_COUNT(kk),LVT_LIS_MOC_IRRIGATEDWATER(kk), &
                  LVT_LISoutput(kk)%head_irrig_list, &
                  1,nsize,(/"kg/m2s"/),1,("-"),&
                  valid_min=(/0.0/),valid_max=(/-9999.0/))

          endif

          call ESMF_ConfigDestroy(modelSpecConfig,rc=rc)

          allocate(LVT_LISoutput(kk)%ptr_into_lsm_list(LVT_LIS_MOC_LSM_COUNT(kk)))
          allocate(LVT_LISoutput(kk)%ptr_into_routing_list(LVT_LIS_MOC_ROUTING_COUNT(kk)))
          allocate(LVT_LISoutput(kk)%ptr_into_irrig_list(LVT_LIS_MOC_IRRIG_COUNT(kk)))
          allocate(LVT_LISoutput(kk)%ptr_into_rtm_list(LVT_LIS_MOC_RTM_COUNT(kk)))


          call set_ptr_into_list(LVT_LIS_MOC_LSM_COUNT(kk), &
               LVT_LISoutput(kk)%head_lsm_list, &
               LVT_LISoutput(kk)%ptr_into_lsm_list)

          call set_ptr_into_list(LVT_LIS_MOC_ROUTING_COUNT(kk), &
               LVT_LISoutput(kk)%head_routing_list, &
               LVT_LISoutput(kk)%ptr_into_routing_list)

          call set_ptr_into_list(LVT_LIS_MOC_RTM_COUNT(kk), &
               LVT_LISoutput(kk)%head_rtm_list, &
               LVT_LISoutput(kk)%ptr_into_rtm_list)

          call set_ptr_into_list(LVT_LIS_MOC_IRRIG_COUNT(kk), &
               LVT_LISoutput(kk)%head_irrig_list, &
               LVT_LISoutput(kk)%ptr_into_irrig_list)

          if(LVT_rc%var_based_strat .gt. 0.and.source.eq.1) then
             strat_var_found = .false. 
             if(LVT_LIS_rc(kk)%anlys_data_class.eq."LSM") then 
                lisdataEntry => LVT_LISoutput(source)%head_lsm_list
             elseif(LVT_LIS_rc(kk)%anlys_data_class.eq."Routing") then 
                lisdataEntry => LVT_LISoutput(source)%head_routing_list
             elseif(LVT_LIS_rc(kk)%anlys_data_class.eq."RTM") then 
                lisdataEntry => LVT_LISoutput(source)%head_rtm_list
             elseif(LVT_LIS_rc(kk)%anlys_data_class.eq."Irrigation") then 
                lisdataEntry => LVT_LISoutput(source)%head_irrig_list
             endif
             
             do while(associated(lisdataEntry))
                if(lisdataEntry%selectOpt.eq.1) then           
                   if(lisdataEntry%short_name.eq.LVT_rc%vname_strat) then 
                      LVT_histdata%strat_varEntry%short_name = &
                           trim(LVT_rc%vname_strat)
                      LVT_histdata%strat_varEntry%vlevels  = &
                           lisdataEntry%vlevels
                      allocate(LVT_histdata%strat_varEntry%value(nsize,&
                           LVT_histdata%strat_varEntry%vlevels))
                      allocate(LVT_histdata%strat_varEntry%count(nsize,&
                           LVT_histdata%strat_varEntry%vlevels))
                      
                      LVT_histdata%strat_varEntry%value = 0 
                      LVT_histdata%strat_varEntry%count = 0 
                      
                      allocate(LVT_stats%strat_var(nsize,&
                           LVT_histdata%strat_varEntry%vlevels))
                      
                      strat_var_found = .true. 
                      exit
                   endif
                endif
                
                lisdataEntry => lisdataEntry%next
             enddo
             if(.not.strat_var_found) then 
                write(LVT_logunit,*) '[ERR] '
                write(LVT_logunit,*) '[ERR] stratification variable '//trim(LVT_rc%vname_strat)
                write(LVT_logunit,*) '[ERR] not found among the LIS output variables'
                call LVT_endrun()
             endif
          endif
          
!Now map the data that was read into the LVT structures
          if(source.eq.1) then 
             dataEntry => LVT_histData%head_ds1_list
          elseif(source.eq.2) then 
             dataEntry => LVT_histData%head_ds2_list
          endif
       
          do while(associated(dataEntry))
!reset the pointers to the head of the linked list
             var_found = .false. 
             
             if(LVT_LIS_rc(source)%anlys_data_class.eq."LSM") then 
                lisdataEntry => LVT_LISoutput(source)%head_lsm_list
             elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Routing") then 
                lisdataEntry => LVT_LISoutput(source)%head_routing_list
             elseif(LVT_LIS_rc(source)%anlys_data_class.eq."RTM") then 
                lisdataEntry => LVT_LISoutput(source)%head_rtm_list
             elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Irrigation") then 
                lisdataEntry => LVT_LISoutput(source)%head_irrig_list
             endif
             do while(associated(lisdataEntry)) 
                if(lisdataEntry%short_name.eq.dataEntry%short_name) then 
                   var_found = .true. 
                   exit
                endif
                lisdataEntry => lisdataEntry%next
             enddo
             !if variable is enabled, but not found in the LIS output
             !then it must be computed
             if(dataEntry%selectNlevs.gt.0.and.(.not.var_found)) then 
                dataEntry%computeVar = 1
             endif
             dataEntry => dataEntry%next
          enddo

       enddo
    endif
  end subroutine LVT_LISoutputInit
     
  
!BOP
! 
! !ROUTINE: get_moc_attributes
!  \label{get_moc_attributes}
!
! !INTERFACE:
subroutine get_moc_attributes(modelSpecConfig, head_dataEntry, &
     short_name,standard_name, long_name, format_type, status)
! 
! !USES:

   implicit none
!
! !ARGUMENTS:
   type(ESMF_Config), intent(inout)        :: modelSpecConfig
   type(LVT_LISmetadataEntry), pointer     :: head_dataEntry
   character(len=*), intent(in)            :: short_name
   character(len=*), intent(in)            :: standard_name
   character(len=*), intent(in)            :: long_name
   character(len=*), intent(in)            :: format_type
   integer         , intent(inout)         :: status

!
! !DESCRIPTION:
! This routine reads the model output configuration attributes for the 
! given dataEntry output data structure. If the variable is selected
! in the attributes file, then the data entry object for that variable
! is initialized and added to the linked list representing all the 
! active variables. 
!
! The arguments are: 
! \begin{description}
!    \item[modelSpecConfig]
!       ESMF handle to the Model Output Configuration file
!    \item[head\_dataEntry]
!       pointer to the linked list containing the variables in the LIS
!       model output
!    \item[short\_name]
!       short name for the output variable
!    \item[standard\_name]
!       standard name for the output variable
!    \item[long\_name]
!       long name for the output variable
!    \item[format\_type]
!       format type for the output variable (scientific -'E', normal - 'F')
!    \item[status]
!       return code (=1 if variable is enabled in the 
!       output attributes table, =0 if it is not)
! \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

   integer :: rc

   integer                              :: selectOpt
   character*20          :: units  !unit of the the variable
   character*20          :: dir
   integer               :: vlevels
   integer               :: varid_def
   integer               :: gribSF          ! grib scale factor
   integer               :: timeAvgOpt
   integer               :: minMaxOpt
   integer               :: stdevOpt

   type(LVT_LISmetadataEntry),     pointer :: current, dataEntry

   if(status.eq.0) then 

      call ESMF_ConfigGetAttribute(modelSpecConfig,selectOpt,&
           default=0,rc=rc)

      call ESMF_ConfigGetAttribute(modelSpecConfig,units,rc=rc)
      
      call ESMF_ConfigGetAttribute(modelSpecConfig,dir,rc=rc)
      call ESMF_ConfigGetAttribute(modelSpecConfig,timeAvgOpt,rc=rc)
      call ESMF_ConfigGetAttribute(modelSpecConfig,minMaxOpt,rc=rc)
      call ESMF_ConfigGetAttribute(modelSpecConfig,stdevOpt,rc=rc)
      call ESMF_ConfigGetAttribute(modelSpecConfig,vlevels,rc=rc)
      call ESMF_ConfigGetAttribute(modelSpecConfig,varid_def,rc=rc)
      call ESMF_ConfigGetAttribute(modelSpecConfig,gribSF,rc=rc)         

      if(selectOpt.eq.1) then 
         status = 1

         allocate(dataEntry)
         if ( .not. associated(head_dataEntry) ) then
            head_dataEntry => dataEntry
         else
            current => head_dataEntry
            do while ( associated(current%next) )
               current => current%next
            enddo
            current%next => dataEntry
         endif
         dataEntry%next => null()
         
         dataEntry%selectOpt  = selectOpt            
         dataEntry%units      = units
         dataEntry%dir        = dir
         dataEntry%timeAvgOpt = timeAvgOpt
         dataEntry%minmaxOpt  = minMaxOpt
         dataEntry%stdevOpt     = stdevOpt
         dataEntry%vlevels    = vlevels
         dataEntry%varid_def     = varid_def
         dataEntry%gribSF     = gribSF
         
         dataEntry%short_name    = short_name
         dataEntry%standard_name = standard_name
         dataEntry%long_name     = long_name
         dataEntry%format        = format_type
         
         dataEntry%computeVar = 0 
         
      else
         status = 0 
      endif
   else
      status = 0 
   endif
 end subroutine get_moc_attributes

!BOP
! 
! !ROUTINE: register_dataEntry
! \label{register_dataEntry}
!
! !INTERFACE: 
  subroutine register_dataEntry(count, lvt_moc_index, &
       head_dataEntry, &
       nunits,nsize,unittypes,ndirs,&
       dirtypes, valid_min, valid_max)

    implicit none
! !ARGUMENTS: 
    integer                             :: count
    integer                             :: lvt_moc_index
    type(LVT_LISmetadataEntry), pointer :: head_dataEntry
    integer                             :: nunits
    integer                             :: nsize
    character(len=*)                    :: unittypes(nunits)
    integer                             :: ndirs
    character(len=*)                    :: dirtypes(ndirs)
    real                                :: valid_min(nunits)
    real                                :: valid_max(nunits)
!EOP

!
! !DESCRIPTION: 
!  This routine initializes the datastructures required for the 
!  specified output variable. 
! 
! 
!EOP

    
    type(LVT_LISmetadataEntry),  pointer :: dataEntry
    
    count = count + 1
    lvt_moc_index = count

    dataEntry => head_dataEntry
    do while(associated(dataEntry%next)) 
       dataEntry =>dataEntry%next
    enddo
    dataEntry%index = lvt_moc_index

    call allocate_dataEntry(dataEntry, nunits, nsize, unittypes, &
         ndirs, dirtypes, valid_min, valid_max)

    if(trim(LVT_rc%vname_strat).eq.trim(dataEntry%standard_name)) then 
       LVT_rc%var_strat_index = lvt_moc_index
    endif

    
  end subroutine register_dataEntry
  
!BOP
! 
! !ROUTINE: allocate_dataEntry
! \label{allocate_dataEntry}
!
! !INTERFACE: 
  subroutine allocate_dataEntry(dataEntry,  &
       nunits,nsize,unittypes, &
       ndirs, dirtypes,valid_min,valid_max)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine initializes the datastructures required for the 
!  specified output variable. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    type(LVT_LISmetadataEntry) :: dataEntry
    integer                 :: nunits
    integer                 :: nsize
    character(len=*)        :: unittypes(nunits)
    integer                 :: ndirs
    character(len=*)        :: dirtypes(ndirs)
    real                    :: valid_min(nunits)
    real                    :: valid_max(nunits)
!EOP
    integer                 :: i 

    if(dataEntry%selectOpt.ne.0.or.dataEntry%computeVar.eq.1) then 
       allocate(dataEntry%value(nsize,dataEntry%vlevels))
       allocate(dataEntry%count(nsize,dataEntry%vlevels))
       
       dataEntry%value = 0 
       dataEntry%count = 0 
       
       dataEntry%nunits = nunits
       allocate(dataEntry%unittypes(nunits))
       allocate(dataEntry%valid_min(nunits))
       allocate(dataEntry%valid_max(nunits))
       
       dataEntry%ndirs = ndirs
       allocate(dataEntry%dirtypes(ndirs))
       
       do i=1,nunits
          dataEntry%unittypes(i) = trim(unittypes(i))
       enddo
       do i=1, ndirs
          dataEntry%dirtypes(i) = trim(dirtypes(i))
       enddo
       do i=1,nunits
          dataEntry%valid_min(i) = valid_min(i)
          dataEntry%valid_max(i) = valid_max(i)
       enddo
       
       dataEntry%stdev_flag = .false.
             
    endif

  end subroutine allocate_dataEntry


!BOP
! !ROUTINE: set_ptr_into_list
! \label{set_ptr_into_list}
! 
! !INTERFACE: 
  subroutine set_ptr_into_list(count, head, array)
     implicit none
! !ARGUMENTS:
     integer                          :: count
     type(LVT_LISmetadataEntry), pointer :: head
     type(lisdep), dimension(count)      :: array
! 
! !DESCRIPTION: 
!  This routine takes an array of pointers and sets each element of the
!  array to point directly to its corresponding element in the given
!  history output linked list.  This allows for direct access to the
!  elements in the history output linked list.
!
!   The arguments are: 
!   \begin{description}
!   \item[count] number of elements in the given history output linked list
!   \item[head] head of the given history output linked list
!   \item[array] array of pointers to point directly to the elements of \newline
!                the given history output linked list
!   \end{description}
!EOP
!   

     type(LVT_LISmetadataEntry), pointer :: dataEntry
     integer :: i

     dataEntry    => head
     
     do i = 1, count
        array(i)%dataEntryPtr => dataEntry
        dataEntry => dataEntry%next
        
     enddo

  end subroutine set_ptr_into_list


!BOP
! 
! !ROUTINE: LVT_readLISModelOutput
! \label{LVT_readLISModelOutput}
!
! !INTERFACE: 
  subroutine LVT_readLISModelOutput(lsmoutfile, source, wout, wopt)
! 
! !USES:

!
! !INPUT PARAMETERS: 
    character(len=*),   intent(in)   :: lsmoutfile
    character(len=*),   intent(in)   :: wout
    integer,            intent(in)   :: source
    character(len=*),   intent(IN), optional :: wopt
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine reads the model output from LIS (binary, grib1, netcdf
!   formats). The routine also computes
!   derived variables (e.g. bowen ratio computed if the LIS output
!   contains latent and sensible heat flux variables).
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer :: ftn
    integer :: nunit
    integer :: iret
    real    :: rand
    integer :: ios

    type(LVT_LISmetadataEntry), pointer :: lisdataEntry

    integer      :: gid, nl, k
    type(LVT_metadataEntry),    pointer :: ebal,wbal, runoff, rnet
    type(LVT_LISmetadataEntry), pointer :: swnet,lwnet,qle,qh,qg
    type(LVT_LISmetadataEntry), pointer :: qf,qa,qv,delsurfheat
    type(LVT_metadataEntry),    pointer :: wrsi,br, ef
    type(LVT_LISmetadataEntry), pointer :: rainf,snowf,qs,qsb
    type(LVT_LISmetadataEntry), pointer :: delswe,delintercept,delsoilmoist
    type(LVT_LISmetadataEntry), pointer :: delcoldcont, delsurfstor
    type(LVT_LISmetadataEntry), pointer :: evap, potevap, et
    type(LVT_LISmetadataEntry), pointer :: swe, precip, soiltemp
    type(LVT_metadataEntry),    pointer :: sweoverp, etoverp,qsoverp, qsboverp
    type(LVT_metadataEntry),    pointer :: roottemp,rootmoist
    type(LVT_metadataEntry),    pointer :: watertabled
    type(LVT_LISmetadataEntry), pointer :: tws
    type(LVT_LISmetadataEntry), pointer :: soilmoist
    real         :: total_depth, value_temp

    if(PRESENT(wopt)) then 
       LVT_LIS_rc(source)%wopt = wopt
    endif

    ftn = 100

    if(trim(wout).eq."binary") then 
       if(LVT_masterproc) then 
          open(ftn,file=trim(lsmoutfile),form='unformatted')
       endif
       call readBinaryOutput(ftn,source)
       if(LVT_masterproc) then 
          close(ftn)
       endif
    elseif(trim(wout).eq."grib1") then 
       if(LVT_masterproc) then

          call grib_open_file(ftn,trim(lsmoutfile),'r',iret)
          if(iret.gt.0) then 
             write(LVT_logunit,*) '[ERR] Opening file ',trim(lsmoutfile), ' failed',iret
             call LVT_endrun()
          endif
          
       endif
       call readGrib1Output(ftn,source)
       if(LVT_masterproc) then 
          call grib_close_file(ftn,iret)
          if(iret.ne.0) then 
             write(LVT_logunit,*) '[ERR] Closing file ',trim(lsmoutfile), ' failed',iret
             call LVT_endrun()
          endif
       endif

    elseif(trim(wout).eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
       ios = nf90_open(path=trim(lsmoutfile),mode=NF90_NOWRITE,ncid=ftn)
       call LVT_verify(ios,'Error opening file'//trim(lsmoutfile))
       
       call readNETCDFOutput(ftn,source)

       ios = nf90_close(ftn)
       call LVT_verify(ios, 'Error in nf90_close')
#endif

    elseif(trim(wout).eq."grib2") then 

    endif

    if(LVT_LIS_rc(source)%anlys_data_class.eq."LSM") then 
       !equation ebal = swnet+lwnet-qle-qh-qg-qa-qv-qf-qa-delsurfheat-delcoldcont
       if(LVT_MOC_EBAL(source).gt.0) then

          if(source.eq.1) then 
             ebal => LVT_histData%ptr_into_ds1_list(&
                  LVT_MOC_EBAL(source))%dataEntryPtr
          else
             ebal => LVT_histData%ptr_into_ds2_list(&
                  LVT_MOC_EBAL(source))%dataEntryPtr
          endif

          if(ebal%selectNlevs.ge.1) then 
             if(LVT_LIS_MOC_SWNET(source).gt.0.and.&
                  LVT_LIS_MOC_LWNET(source).gt.0.and.&
                  LVT_LIS_MOC_QLE(source).gt.0.and.&
                  LVT_LIS_MOC_QH(source).gt.0.and.&
                  LVT_LIS_MOC_QG(source).gt.0) then 
                swnet => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_SWNET(source))%dataEntryPtr
                lwnet => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_LWNET(source))%dataEntryPtr
                qle => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QLE(source))%dataEntryPtr
                qh => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QH(source))%dataEntryPtr
                qg => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QG(source))%dataEntryPtr
             else
                write(LVT_logunit,*)&
                     '[ERR] Please make sure that all basic components of energy balance'
                write(LVT_logunit,*)&
                     '[ERR] (SWnet, LWnet, Qle, Qh, Qg) are enabled'
                call LVT_endrun()
             endif

             if(LVT_LIS_MOC_QF(source).gt.0) then 
                qf => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QF(source))%dataEntryPtr
             endif
             if(LVT_LIS_MOC_QV(source).gt.0) then 
                qv => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QV(source))%dataEntryPtr
             endif
             if(LVT_LIS_MOC_QA(source).gt.0) then 
                qa => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QA(source))%dataEntryPtr
             endif
             if(LVT_LIS_MOC_DELSURFHEAT(source).gt.0) then 
                delsurfheat => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_DELSURFHEAT(source))%dataEntryPtr
             endif
             if(LVT_LIS_MOC_DELCOLDCONT(source).gt.0) then 
                delcoldcont => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_DELCOLDCONT(source))%dataEntryPtr
             endif

             do gid=1,LVT_rc%npts
                ebal%value(gid,1) = &
                     (swnet%value(gid,1)+&
                     lwnet%value(gid,1)-&
                     qle%value(gid,1)-&
                     qh%value(gid,1)+&
                     qg%value(gid,1))

                if(LVT_LIS_MOC_QF(source).gt.0) then 
                   ebal%value(gid,1) = &
                        ebal%value(gid,1) -&
                        qf%value(gid,1)
                endif
                if(LVT_LIS_MOC_QV(source).gt.0) then 
                   ebal%value(gid,1) = &
                        ebal%value(gid,1) -&
                        qv%value(gid,1)
                endif
                if(LVT_LIS_MOC_QA(source).gt.0) then 
                   ebal%value(gid,1) = &
                        ebal%value(gid,1) -&
                        qa%value(gid,1)
                endif
                if(LVT_LIS_MOC_DELSURFHEAT(source).gt.0) then 
                   ebal%value(gid,1) = &
                        ebal%value(gid,1) -&
                        delsurfheat%value(gid,1)
                endif
                if(LVT_LIS_MOC_DELCOLDCONT(source).gt.0) then 
                   ebal%value(gid,1) = &
                        ebal%value(gid,1) -&
                        delcoldcont%value(gid,1)
                endif
             enddo
          endif
       endif
       !equation ebal = rainf+snowf-evap-qs-qsb -(delswe+delintercept+delsoilmoist+delsurfstor)
       if(LVT_MOC_WBAL(source).gt.0) then
          if(source.eq.1) then 
             wbal => LVT_histData%ptr_into_ds1_list(&
                  LVT_MOC_WBAL(source))%dataEntryPtr
          else
             wbal => LVT_histData%ptr_into_ds2_list(&
                  LVT_MOC_WBAL(source))%dataEntryPtr
          endif

          if(wbal%selectNlevs.ge.1) then 
             if (LVT_LIS_MOC_RAINF(source).gt.0.and.&
                  LVT_LIS_MOC_SNOWF(source).gt.0.and.&
                  LVT_LIS_MOC_EVAP(source).gt.0.and.&
                  LVT_LIS_MOC_QS(source).gt.0.and.&
                  LVT_LIS_MOC_QSB(source).gt.0) then 
                rainf => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_RAINF(source))%dataEntryPtr
                snowf => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_SNOWF(source))%dataEntryPtr
                evap => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_EVAP(source))%dataEntryPtr
                qs => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QS(source))%dataEntryPtr
                qsb => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QSB(source))%dataEntryPtr
             else
                write(LVT_logunit,*)&
                     '[ERR] Please make sure that all basic components of water balance'
                write(LVT_logunit,*)&
                     '[ERR] (Rainf, Snowf, Evap, Qs, Qsb) are enabled'
                call LVT_endrun()
             endif

             if(LVT_LIS_MOC_DELSWE(source).gt.0) then 
                delswe => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_DELSWE(source))%dataEntryPtr
             endif
             if(LVT_LIS_MOC_DELSWE(source).gt.0) then 
                delsoilmoist => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_DELSOILMOIST(source))%dataEntryPtr
             endif
             if(LVT_LIS_MOC_DELINTERCEPT(source).gt.0) then 
                delintercept => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_DELINTERCEPT(source))%dataEntryPtr
             endif
             if(LVT_LIS_MOC_DELSURFSTOR(source).gt.0) then 
                delsurfstor => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_DELSURFSTOR(source))%dataEntryPtr
             endif

             do gid=1,LVT_rc%npts
                wbal%value(gid,1) = &
                     (rainf%value(gid,1)+&
                     snowf%value(gid,1)-&
                     evap%value(gid,1)-&
                     qs%value(gid,1)+&
                     qsb%value(gid,1))
                if(LVT_LIS_MOC_DELSWE(source).gt.0) then
                   wbal%value(gid,1) = &
                        wbal%value(gid,1) -&
                        delswe%value(gid,1)
                endif
                if(LVT_LIS_MOC_DELSOILMOIST(source).gt.0) then 
                   wbal%value(gid,1) = &
                        wbal%value(gid,1) -&
                        delsoilmoist%value(gid,1)
                endif
                if(LVT_LIS_MOC_DELINTERCEPT(source).gt.0) then 
                   wbal%value(gid,1) = &
                        wbal%value(gid,1) -&
                        delintercept%value(gid,1)
                endif
                if(LVT_LIS_MOC_DELSURFSTOR(source).gt.0) then 
                   wbal%value(gid,1) = &
                        wbal%value(gid,1) -&
                        delsurfstor%value(gid,1)
                endif
             enddo
          endif
       endif
       
       if(LVT_MOC_RNET(source).gt.0) then 
          if(source.eq.1) then 
             rnet => LVT_histData%ptr_into_ds1_list(&
                  LVT_MOC_RNET(source))%dataEntryPtr
          else
             rnet => LVT_histData%ptr_into_ds2_list(&
                  LVT_MOC_RNET(source))%dataEntryPtr
          endif
          if(rnet%selectNlevs.ge.1) then 
             if(LVT_LIS_MOC_QLE(source).gt.0.and.&
                  LVT_LIS_MOC_QH(source).gt.0.and.&
                  LVT_LIS_MOC_QG(source).gt.0) then 
                qle => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QLE(source))%dataEntryPtr
                qh => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QH(source))%dataEntryPtr
                qg => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QG(source))%dataEntryPtr
             else
                write(LVT_logunit,*)& 
                     '[ERR] Please enable both Qle, Qh and Qg in the LIS output to compute Rnet'
                call LVT_endrun()
             endif
             do gid=1,LVT_rc%npts
                rnet%value(gid,1) =&
                     qle%value(gid,1) + qh%value(gid,1) + & 
                     qg%value(gid,1)
                rnet%count(gid,1) =  qle%count(gid,1)
             enddo
          endif
       endif

       if(LVT_MOC_BR(source).gt.0) then 
          if(source.eq.1) then 
             br => LVT_histData%ptr_into_ds1_list(&
                  LVT_MOC_BR(source))%dataEntryPtr
          else
             br => LVT_histData%ptr_into_ds2_list(&
                  LVT_MOC_BR(source))%dataEntryPtr             
          endif

          if(br%selectNlevs.ge.1) then 
             if(LVT_LIS_MOC_QLE(source).gt.0.and.LVT_LIS_MOC_QH(source).gt.0) then 
                qle => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QLE(source))%dataEntryPtr
                qh => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QH(source))%dataEntryPtr
             else
                write(LVT_logunit,*)& 
                     '[ERR] Please enable both Qle and Qh in the LIS output to compute bowen ratio'
                call LVT_endrun()
             endif
             do gid=1,LVT_rc%npts
                if(qle%value(gid,1).ne.0) then 
                   br%value(gid,1) =&
                        qh%value(gid,1)/&
                        qle%value(gid,1)
                   br%count(gid,1) = &
                        qle%count(gid,1)
                else
                   br%value(gid,1) = LVT_rc%udef
                endif
             enddo
          endif
       endif

       if(LVT_MOC_RUNOFF(source).gt.0) then 
          if(source.eq.1) then 
             runoff => LVT_histData%ptr_into_ds1_list(&
                  LVT_MOC_RUNOFF(source))%dataEntryPtr
          else
             runoff => LVT_histData%ptr_into_ds2_list(&
                  LVT_MOC_RUNOFF(source))%dataEntryPtr
          endif
          if(runoff%selectNlevs.ge.1) then 
             if(LVT_LIS_MOC_QS(source).gt.0.and.LVT_LIS_MOC_QSB(source).gt.0) then 
                qs => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QS(source))%dataEntryPtr
                qsb => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QSB(source))%dataEntryPtr
             else
                write(LVT_logunit,*)& 
                     '[ERR] Please enable both Qs and Qsb in the LIS output to compute total runoff'
                call LVT_endrun()
             endif
             do gid=1,LVT_rc%npts
                runoff%value(gid,1) =&
                     qs%value(gid,1) + & 
                     qsb%value(gid,1)
                runoff%count(gid,1) = &
                     qs%count(gid,1)
             enddo
          endif
       endif

       if(LVT_MOC_EF(source).gt.0) then 
          if(source.eq.1) then              
             ef => LVT_histData%ptr_into_ds1_list(&
                  LVT_MOC_EF(source))%dataEntryPtr
          else
             ef => LVT_histData%ptr_into_ds2_list(&
                  LVT_MOC_EF(source))%dataEntryPtr
          endif
          if(ef%selectNlevs.ge.1) then 
             if(LVT_LIS_MOC_QLE(source).gt.0.and.LVT_LIS_MOC_QH(source).gt.0) then 
                qle => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QLE(source))%dataEntryPtr
                qh => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QH(source))%dataEntryPtr
             else
                write(LVT_logunit,*)& 
                     '[ERR] Please enable both Qle and Qh in the LIS output to compute evaporative fraction'
                call LVT_endrun()
             endif

             do gid=1,LVT_rc%npts
                if((qle%value(gid,1)+&
                     qh%value(gid,1)).ne.0) then 
                   ef%value(gid,1) =&
                        qle%value(gid,1)/&
                        (qh%value(gid,1)+&
                        qle%value(gid,1))
                   ef%count(gid,1) = &
                        qle%count(gid,1)
                else
                   ef%value(gid,1) = LVT_rc%udef
                endif
             enddo
          endif
       endif
       if(LVT_MOC_SWEOVERP(source).gt.0) then 
          if(source.eq.1) then 
             sweoverp => LVT_histData%ptr_into_ds1_list(&
                  LVT_MOC_SWEOVERP(source))%dataEntryPtr
          else
             sweoverp => LVT_histData%ptr_into_ds2_list(&
                  LVT_MOC_SWEOVERP(source))%dataEntryPtr
          endif
          if(sweoverp%selectNlevs.ge.1) then 
             if(LVT_LIS_MOC_SWE(source).gt.0.and.LVT_LIS_MOC_TOTALPRECIP(source).gt.0) then 
                swe => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_SWE(source))%dataEntryPtr
                precip => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_TOTALPRECIP(source))%dataEntryPtr

                if(swe%units.ne.precip%units) then 
                   write(LVT_logunit,*) '[ERR] The units of SWE and P are inconsistent ..'
                   call LVT_endrun()
                endif

             else
                write(LVT_logunit,*)& 
                     '[ERR] Please enable both SWE and Totalprecip in the LIS output to compute SWE/P'
                call LVT_endrun()
             endif
             do gid=1,LVT_rc%npts
                if(precip%value(gid,1).ne.0) then 
                   sweoverp%value(gid,1) =&
                        swe%value(gid,1)/&
                        precip%value(gid,1)
                   sweoverp%count(gid,1) = &
                        swe%count(gid,1)
                else
                   sweoverp%value(gid,1) = LVT_rc%udef
                endif
             enddo
          endif
       endif

       if(LVT_MOC_ETOVERP(source).gt.0) then 
          if(source.eq.1) then 
             etoverp => LVT_histData%ptr_into_ds1_list(&
               LVT_MOC_ETOVERP(source))%dataEntryPtr
          else
             etoverp => LVT_histData%ptr_into_ds2_list(&
               LVT_MOC_ETOVERP(source))%dataEntryPtr
          endif
          if(etoverp%selectNlevs.ge.1) then 
             if(LVT_LIS_MOC_EVAP(source).gt.0.and.LVT_LIS_MOC_TOTALPRECIP(source).gt.0) then 
                et => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_EVAP(source))%dataEntryPtr
                precip => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_TOTALPRECIP(source))%dataEntryPtr

                if(et%units.ne.precip%units) then 
                   write(LVT_logunit,*) '[ERR] The units of ET and P are inconsistent ..'
                   call LVT_endrun()
                endif

             else
                write(LVT_logunit,*)& 
                     '[ERR] Please enable both ET and Totalprecip in the LIS output to compute ET/P'
                call LVT_endrun()
             endif
             do gid=1,LVT_rc%npts
                if(precip%value(gid,1).ne.0) then 
                   etoverp%value(gid,1) =&
                        et%value(gid,1)/&
                        precip%value(gid,1)
                   etoverp%count(gid,1) = &
                        et%count(gid,1)
                else
                   etoverp%value(gid,1) = LVT_rc%udef
                endif
             enddo
          endif
       endif

       if(LVT_MOC_QSOVERP(source).gt.0) then 
          if(source.eq.1) then 
             qsoverp => LVT_histData%ptr_into_ds1_list(&
                  LVT_MOC_QSOVERP(source))%dataEntryPtr
          else
             qsoverp => LVT_histData%ptr_into_ds2_list(&
                  LVT_MOC_QSOVERP(source))%dataEntryPtr             
          endif

          if(qsoverp%selectNlevs.ge.1) then 
             if(LVT_LIS_MOC_QS(source).gt.0.and.LVT_LIS_MOC_TOTALPRECIP(source).gt.0) then 
                qs => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QS(source))%dataEntryPtr
                precip => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_TOTALPRECIP(source))%dataEntryPtr

                if(qs%units.ne.precip%units) then 
                   write(LVT_logunit,*) '[ERR] The units of Qs and P are inconsistent ..'
                   call LVT_endrun()
                endif

             else
                write(LVT_logunit,*)& 
                     '[ERR] Please enable both Qs and Totalprecip in the LIS output to compute Qs/P'
                call LVT_endrun()
             endif
             do gid=1,LVT_rc%npts
                if(precip%value(gid,1).ne.0) then 
                   qsoverp%value(gid,1) =&
                        qs%value(gid,1)/&
                        precip%value(gid,1)
                   qsoverp%count(gid,1) = &
                        qs%count(gid,1)
                else
                   qsoverp%value(gid,1) = LVT_rc%udef
                endif
             enddo
          endif
       endif

       if(LVT_MOC_QSBOVERP(source).gt.0) then 
          if(source.eq.1) then 
             qsboverp => LVT_histData%ptr_into_ds1_list(&
                  LVT_MOC_QSBOVERP(source))%dataEntryPtr
          else
             qsboverp => LVT_histData%ptr_into_ds2_list(&
                  LVT_MOC_QSBOVERP(source))%dataEntryPtr
          endif
          if(qsboverp%selectNlevs.ge.1) then 
             if(LVT_LIS_MOC_QSB(source).gt.0.and.LVT_LIS_MOC_TOTALPRECIP(source).gt.0) then 
                qsb => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_QSB(source))%dataEntryPtr
                precip => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_TOTALPRECIP(source))%dataEntryPtr

                if(qsb%units.ne.precip%units) then 
                   write(LVT_logunit,*) '[ERR] The units of Qsb and P are inconsistent ..'
                   call LVT_endrun()
                endif

             else
                write(LVT_logunit,*)& 
                     '[ERR] Please enable both Qsb and Totalprecip in the LIS output to compute Qsb/P'
                call LVT_endrun()
             endif
             do gid=1,LVT_rc%npts
                if(precip%value(gid,1).ne.0) then 
                   qsboverp%value(gid,1) =&
                        qsb%value(gid,1)/&
                        precip%value(gid,1)
                   qsboverp%count(gid,1) = &
                        qsb%count(gid,1)
                else
                   qsboverp%value(gid,1) = LVT_rc%udef
                endif
             enddo
          endif
       endif

       !now compute the root zone variables, if necessary
       !weighted by thicknesses.. 
       if(LVT_MOC_ROOTMOIST(source).gt.0) then 
          if(source.eq.1) then 
             rootmoist => LVT_histData%ptr_into_ds1_list(&
                  LVT_MOC_ROOTMOIST(source))%dataEntryPtr
          else
             rootmoist => LVT_histData%ptr_into_ds2_list(&
                  LVT_MOC_ROOTMOIST(source))%dataEntryPtr
          endif

          if(rootmoist%selectNlevs.ge.1) then 
             soilmoist => LVT_LISoutput(source)%ptr_into_lsm_list(&
                  LVT_LIS_MOC_SOILMOIST(source))%dataEntryPtr

!------------------------------------------------------------------------------
!   If LIS output is involved, we can calculate root zone, else the plugin is
!   expected to compute this directly. 
!------------------------------------------------------------------------------
             if(LVT_LIS_rc(source)%model_name.eq."CLSM F2.5") then 
                ! for CLSM, simply pick out the second layer. 
                do gid=1,LVT_rc%npts
                   rootmoist%value(gid,1) = soilmoist%value(gid,2)
                   rootmoist%count(gid,1) = &
                        soilmoist%count(gid,2)
                enddo
             else
                total_depth =0 
                do k=1,LVT_LIS_rc(source)%nsmlayers
                   total_depth = total_depth + LVT_LIS_rc(source)%smthick(k)
                   if(total_depth.ge.LVT_rc%lis_rz_d) then 
                      nl = k
                      exit
                   endif
                enddo
                
                do gid=1,LVT_rc%npts
                   value_temp = 0 
                   total_depth = 0 
                   do k=1,nl
                      value_temp = value_temp + &
                           soilmoist%value(gid,k)*&
                           LVT_LIS_rc(source)%smthick(k)
                      total_depth = total_depth + LVT_LIS_rc(source)%smthick(k)
                   enddo
                   rootmoist%value(gid,1) = value_temp/total_depth
                   rootmoist%count(gid,1) = &
                        soilmoist%count(gid,1)
                enddo
             endif
          endif
       endif

       if(LVT_MOC_ROOTTEMP(source).gt.0) then 
          if(source.eq.1) then 
             roottemp => LVT_histData%ptr_into_ds1_list(&
                  LVT_MOC_ROOTTEMP(source))%dataEntryPtr
          else
             roottemp => LVT_histData%ptr_into_ds2_list(&
                  LVT_MOC_ROOTTEMP(source))%dataEntryPtr             
          endif
          if(roottemp%selectNlevs.ge.1) then 
             soiltemp => LVT_LISoutput(source)%ptr_into_lsm_list(&
                  LVT_LIS_MOC_SOILTEMP(source))%dataEntryPtr
!------------------------------------------------------------------------------
!   If LIS output is involved, we can calculate root zone, else the plugin is
!   expected to compute this directly. 
!------------------------------------------------------------------------------
             if(LVT_rc%lis_output_obs) then 
                total_depth =0 
                do k=1,LVT_LIS_rc(source)%nstlayers
                   total_depth = total_depth + LVT_LIS_rc(source)%stthick(k)
                   if(total_depth.ge.LVT_rc%lis_rz_d) then 
                      nl = k
                      exit
                   endif
                enddo

                do gid=1,LVT_rc%npts
                   value_temp = 0 
                   total_depth = 0 
                   do k=1,nl
                      value_temp = value_temp + soiltemp%value(gid,k)*&
                           LVT_LIS_rc(source)%stthick(k)
                      total_depth = total_depth + LVT_LIS_rc(source)%stthick(k)
                   enddo
                   roottemp%value(gid,1) = value_temp/total_depth
                   roottemp%count(gid,1) = &
                        soiltemp%count(gid,1)
                enddo
             endif
          endif
       endif

       if(LVT_MOC_WRSI(source).gt.0) then 
          if(source.eq.1) then 
             wrsi => LVT_histData%ptr_into_ds1_list(&
                  LVT_MOC_WRSI(source))%dataEntryPtr
          else
             wrsi => LVT_histData%ptr_into_ds2_list(&
                  LVT_MOC_WRSI(source))%dataEntryPtr
          endif
          if(wrsi%selectNlevs.ge.1) then 
             if(LVT_LIS_MOC_EVAP(source).gt.0.and.LVT_LIS_MOC_POTEVAP(source).gt.0) then 
                evap => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_EVAP(source))%dataEntryPtr
                potevap => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_POTEVAP(source))%dataEntryPtr
             else
                write(LVT_logunit,*) &
                     '[ERR] Please enable both EVAP and POTEVAP in the LIS output to compute WRSI'
                call LVT_endrun()
             endif
             do gid=1,LVT_rc%npts
                if(potevap%value(gid,1).ne.0) then 
                   wrsi%value(gid,1) =&
                        evap%value(gid,1)/&
                        potevap%value(gid,1)
                   wrsi%count(gid,1) = &
                        potevap%count(gid,1)
                else
                   wrsi%value(gid,1) = LVT_rc%udef
                endif
             enddo
          endif
       endif

       if(LVT_MOC_WATERTABLED(source).gt.0) then 

          if(source.eq.1) then 
             watertabled => LVT_histData%ptr_into_ds1_list(&
                  LVT_MOC_WATERTABLED(source))%dataEntryPtr
          else
             watertabled => LVT_histData%ptr_into_ds2_list(&
                  LVT_MOC_WATERTABLED(source))%dataEntryPtr
          endif
          
          if(watertabled%selectNlevs.gt.0) then 
             if(LVT_LIS_MOC_TWS(source).gt.0.and.LVT_LIS_MOC_SWE(source).gt.0.and.&
                  LVT_LIS_MOC_SOILMOIST(source).gt.0) then 
                tws => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_TWS(source))%dataEntryPtr
                swe => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_SWE(source))%dataEntryPtr
                soilmoist => LVT_LISoutput(source)%ptr_into_lsm_list(&
                     LVT_LIS_MOC_SOILMOIST(source))%dataEntryPtr                  
             else
                write(LVT_logunit,*) &
                     '[ERR] Please enable both TWS, SWE and SoilMoist in the LIS output to compute WaterTableD'
                call LVT_endrun()
             endif
             do gid=1,LVT_rc%npts
                ! for CLSM, simply pick out the second layer and here 
                ! we assume the units of meter for the watertable depth

                watertabled%value(gid,1) =&
                     (tws%value(gid,1) - &
                     soilmoist%value(gid,2)*1000.0 - &
                     swe%value(gid,1))/1000.0
                watertabled%count(gid,1) = &
                     tws%count(gid,1)
             enddo
          endif
       endif
    endif
  end subroutine LVT_readLISModelOutput


!BOP
! 
! !ROUTINE: readBinaryOutput
! \label{readBinaryOutput}
!
! !INTERFACE: 
  subroutine readBinaryOutput(ftn,source)
! 
! !USES: 

!
! !INPUT PARAMETERS: 
    integer,   intent(in)   :: ftn
    integer,   intent(in)   :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine reads a binary output file based on the list of selected 
!  output variables. 
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest \newline
!    \item[ftn] file unit for the output file \newline
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LVT\_readSingleBinaryVar](\ref{LVT_readSingleBinaryVar}) \newline
!     reads a single variable from a flat binary file. 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    type(LVT_metadataEntry), pointer :: dataEntry
    type(LVT_LISmetadataEntry), pointer :: lisdataEntry

    if(LVT_LIS_rc(source)%anlys_data_class.eq."LSM") then 
       lisdataEntry => LVT_LISoutput(source)%head_lsm_list
    elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Routing") then 
       lisdataEntry => LVT_LISoutput(source)%head_routing_list
    elseif(LVT_LIS_rc(source)%anlys_data_class.eq."RTM") then 
       lisdataEntry => LVT_LISoutput(source)%head_rtm_list
    elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Irrigation") then 
       lisdataEntry => LVT_LISoutput(source)%head_irrig_list
    endif
 
    do while(associated(lisdataEntry))
       call LVT_readSingleBinaryVar(ftn,source, lisdataEntry)
       lisdataEntry => lisdataEntry%next
    enddo

    if(source.eq.1) then 
       dataEntry => LVT_histData%head_ds1_list
    elseif(source.eq.2) then 
       dataEntry => LVT_histData%head_ds2_list
    endif
       
    do while(associated(dataEntry))
!reset the pointers to the head of the linked list
       if(LVT_LIS_rc(source)%anlys_data_class.eq."LSM") then 
          lisdataEntry => LVT_LISoutput(source)%head_lsm_list
       elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Routing") then 
          lisdataEntry => LVT_LISoutput(source)%head_routing_list
       elseif(LVT_LIS_rc(source)%anlys_data_class.eq."RTM") then 
          lisdataEntry => LVT_LISoutput(source)%head_rtm_list
       elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Irrigation") then 
          lisdataEntry => LVT_LISoutput(source)%head_irrig_list
       endif
       do while(associated(lisdataEntry)) 
          if(lisdataEntry%short_name.eq.dataEntry%short_name) then              
             call mapLISdataToLVT(source, dataEntry, lisdataEntry)
             exit
          endif
          lisdataEntry => lisdataEntry%next
       enddo
       dataEntry => dataEntry%next
    enddo

  end subroutine readBinaryOutput

!BOP
! 
! !ROUTINE: readGrib1Output
! \label{readGrib1Output}
!
! !INTERFACE: 
  subroutine readGrib1Output(ftn,source)
! 
! !USES: 

    implicit none 

    integer, intent(in)   :: source
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine reads a binary output file based on the list of selected 
!  output variables. 
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest \newline
!    \item[ftn] file unit for the output file \newline
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LVT\_readSingleGrib1Var](\ref{LVT_readSingleGrib1Var}) \newline
!     reads a single variable from a flat binary file. 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    integer   :: ftn
!EOP
    integer          :: nvars
    integer          :: index
    integer          :: igrib
    integer, allocatable :: pid(:)
    integer, allocatable :: tid(:)
    integer, allocatable :: lid(:)
    integer          :: iret
    real             :: undef_v
    real, allocatable    :: var(:,:) 
    type(LVT_metadataEntry), pointer :: dataEntry
    type(LVT_LISmetadataEntry), pointer :: lisdataEntry

    call grib_count_in_file(ftn,nvars)
    allocate(var(nvars,LVT_LIS_rc(source)%lnc*LVT_LIS_rc(source)%lnr))
    allocate(pid(nvars))
    allocate(tid(nvars))
    allocate(lid(nvars))

    do index = 1, nvars
       call grib_new_from_file(ftn,igrib,iret)
       call LVT_verify(iret,'grib_new_from_file failed in LVT_histDataMod')
       
       call grib_get(igrib,"indicatorOfParameter",pid(index),iret)
       call LVT_verify(iret,'grib_get failed for indicatorOfParameter')

       call grib_get(igrib, "timeRangeIndicator",tid(index),iret)
       call LVT_verify(iret,'grib_get failed for timeRangeIndicator')

       call grib_get(igrib, "bottomLevel",lid(index),iret)
       call LVT_verify(iret,'grib_get failed for bottomLevel')
       
       call grib_get(igrib,"values",var(index,:),iret)
       call LVT_verify(iret,'grib_get failed for values')

       call grib_get(igrib,"missingValue",undef_v,iret)
       call LVT_verify(iret,'grib_get failed for values')

       call grib_release(igrib,iret)
       call LVT_verify(iret,'grib_release failed in LVT_histDataMod')
    enddo

    if(LVT_LIS_rc(source)%anlys_data_class.eq."LSM") then 
       lisdataEntry => LVT_LISoutput(source)%head_lsm_list
    elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Routing") then 
       lisdataEntry => LVT_LISoutput(source)%head_routing_list
    elseif(LVT_LIS_rc(source)%anlys_data_class.eq."RTM") then 
       lisdataEntry => LVT_LISoutput(source)%head_rtm_list
    elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Irrigation") then 
       lisdataEntry => LVT_LISoutput(source)%head_irrig_list
    endif

    do while(associated(lisdataEntry))
       call LVT_readSingleGrib1Var(ftn,source,lisdataEntry, &
            nvars, var,undef_v,pid,tid,lid)
       lisdataEntry => lisdataEntry%next
    enddo

    if(source.eq.1) then 
       dataEntry => LVT_histData%head_ds1_list
    elseif(source.eq.2) then 
       dataEntry => LVT_histData%head_ds2_list
    endif
       
    do while(associated(dataEntry))
!reset the pointers to the head of the linked list
       if(LVT_LIS_rc(source)%anlys_data_class.eq."LSM") then 
          lisdataEntry => LVT_LISoutput(source)%head_lsm_list
       elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Routing") then 
          lisdataEntry => LVT_LISoutput(source)%head_routing_list
       elseif(LVT_LIS_rc(source)%anlys_data_class.eq."RTM") then 
          lisdataEntry => LVT_LISoutput(source)%head_rtm_list
       elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Irrigation") then 
          lisdataEntry => LVT_LISoutput(source)%head_irrig_list
       endif
       do while(associated(lisdataEntry)) 
          if(lisdataEntry%short_name.eq.dataEntry%short_name) then              
             call mapLISdataToLVT(source,dataEntry, lisdataEntry)
             exit
          endif
          lisdataEntry => lisdataEntry%next
       enddo
       dataEntry => dataEntry%next
    enddo

    deallocate(var)
    deallocate(pid)
    deallocate(tid)
    deallocate(lid)

  end subroutine readGrib1Output

!BOP
! 
! !ROUTINE: readNETCDFOutput
! \label{readNETCDFOutput}
!
! !INTERFACE: 
  subroutine readNETCDFOutput(ftn,source)
! 
! !USES: 
!
! !INPUT PARAMETERS: 
    integer,   intent(in)   :: ftn
    integer,   intent(in)   :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine reads a NETCDF output file based on the list of selected 
!  output variables. 
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest \newline
!    \item[ftn] file unit for the output file \newline
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LVT\_readSingleNETCDFVar](\ref{LVT_readSingleNETCDFVar}) \newline
!     reads a single variable from a flat binary file. 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    type(LVT_metadataEntry), pointer    :: dataEntry
    type(LVT_LISmetadataEntry), pointer :: lisdataEntry
! Read the LIS output
    if(LVT_LIS_rc(source)%anlys_data_class.eq."LSM") then 
       lisdataEntry => LVT_LISoutput(source)%head_lsm_list
    elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Routing") then 
       lisdataEntry => LVT_LISoutput(source)%head_routing_list
    elseif(LVT_LIS_rc(source)%anlys_data_class.eq."RTM") then 
       lisdataEntry => LVT_LISoutput(source)%head_rtm_list
    elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Irrigation") then 
       lisdataEntry => LVT_LISoutput(source)%head_irrig_list
    endif

    do while(associated(lisdataEntry))
       call LVT_readSingleNETCDFVar(ftn,source, lisdataEntry)
       lisdataEntry => lisdataEntry%next
    enddo

!Now map the data that was read into the LVT structures
    if(source.eq.1) then 
       dataEntry => LVT_histData%head_ds1_list
    elseif(source.eq.2) then 
       dataEntry => LVT_histData%head_ds2_list
    endif
       
    do while(associated(dataEntry))
!reset the pointers to the head of the linked list
       if(LVT_LIS_rc(source)%anlys_data_class.eq."LSM") then 
          lisdataEntry => LVT_LISoutput(source)%head_lsm_list
       elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Routing") then 
          lisdataEntry => LVT_LISoutput(source)%head_routing_list
       elseif(LVT_LIS_rc(source)%anlys_data_class.eq."RTM") then 
          lisdataEntry => LVT_LISoutput(source)%head_rtm_list
       elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Irrigation") then 
          lisdataEntry => LVT_LISoutput(source)%head_irrig_list
       endif
       do while(associated(lisdataEntry)) 
          if(lisdataEntry%short_name.eq.dataEntry%short_name) then              
             call mapLISdataToLVT(source,dataEntry, lisdataEntry)
             exit
          endif
          lisdataEntry => lisdataEntry%next
       enddo
       dataEntry => dataEntry%next
    enddo

  end subroutine readNETCDFOutput

!BOP
! 
! !ROUTINE: LVT_readSingleBinaryVar
! \label{LVT_readSingleBinaryVar}
!
! !INTERFACE:
  subroutine LVT_readSingleBinaryVar(ftn, source, dataEntry)
! 
! !USES: 

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine reads a single variable from the LIS output file in the 
!   binary format. Based on the external datamask, the routine also filters
!   each variable. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    integer :: ftn
    integer :: source
    type(LVT_LISmetadataEntry) :: dataEntry
!EOP
    real, allocatable :: value1d(:)
    real, allocatable :: value2d(:,:)
    real, allocatable :: value2d_ip(:,:)
    integer :: unit_id
    integer :: k,i,c,r,t,m,gid

    unit_id = -1
    if(dataEntry%selectOpt.eq.1) then 
       do i=1,dataEntry%nunits
          if(trim(dataEntry%units).eq.trim(dataEntry%unittypes(i))) then 
             unit_id = i
             exit
          endif
       enddo    
       if(unit_id.eq.-1) then 
          write(LVT_logunit,*) '[ERR] routine to diagnose ',trim(dataEntry%standard_name),&
               ' in units of ',trim(dataEntry%units),' is not defined'
          write(LVT_logunit,*) '[ERR] for diagnostic output...'
          call LVT_endrun()
       endif
       do k=1,dataEntry%vlevels
          if(trim(LVT_LIS_rc(source)%wopt).eq."1d tilespace") then 
             
             allocate(value1d(LVT_LIS_rc(source)%ntiles))
             allocate(value2d(LVT_LIS_rc(source)%lnc,LVT_LIS_rc(source)%lnr))
             allocate(value2d_ip(LVT_rc%lnc,LVT_rc%lnr))

             read(ftn) value1d

             if(LVT_rc%computeEnsMetrics.eq.1) then 
                do t=1,LVT_LIS_rc(source)%ntiles
                   c = LVT_LIS_domain(source)%tile(t)%col
                   r = LVT_LIS_domain(source)%tile(t)%row
                   if(LVT_stats%datamask(c,r).eq.1) then 
                      gid = LVT_domain%gindex(c,r)
                      if(gid.ne.-1) then 
                         if(value1d(t).ne.LVT_rc%udef) then 
                            dataEntry%value(t,k) = &
                                 dataEntry%value(t,k) + &
                                 value1d(t)
                            dataEntry%count(t,k) = &
                                 dataEntry%count(t,k) + 1
                         endif
                      endif
                   endif
                enddo
             else
                value2d = 0
                do i=1,LVT_LIS_rc(source)%ntiles,LVT_LIS_rc(source)%nensem
                   do m=1,LVT_LIS_rc(source)%nensem
                      t = i+m-1
                      c = LVT_LIS_domain(source)%tile(t)%col
                      r = LVT_LIS_domain(source)%tile(t)%row
                      value2d(c,r) = value2d(c,r) + value1d(t)*&
                           LVT_LIS_domain(source)%tile(t)%fgrd/&
                           LVT_LIS_rc(source)%nensem
                   enddo
                enddo
                
                call interp2lisgrid_2d(source, value2d,value2d_ip)
                
                do r=1,LVT_rc%lnr
                   do c=1,LVT_rc%lnc
                      if(LVT_stats%datamask(c,r).eq.1) then 
                         gid = LVT_domain%gindex(c,r)
                         if(gid.ne.-1) then 
                            if(value2d_ip(c,r).ne.LVT_rc%udef) then 
                               dataEntry%value(gid,k) = &
                                    dataEntry%value(gid,k) + &
                                    value2d_ip(c,r)
                               dataEntry%count(gid,k) = &
                                    dataEntry%count(gid,k) + 1
                            endif
                         endif
                      endif
                   enddo
                enddo
             endif

             deallocate(value1d)
             deallocate(value2d)
             deallocate(value2d_ip)

          elseif(trim(LVT_LIS_rc(source)%wopt).eq."2d gridspace") then 

             allocate(value2d(LVT_LIS_rc(source)%lnc, LVT_LIS_rc(source)%lnr))
             allocate(value2d_ip(LVT_rc%lnc, LVT_rc%lnr))

             value2d = LVT_rc%udef

             call LVT_readvar_gridded(ftn,source,value2d)
             call interp2lisgrid_2d(source,value2d,value2d_ip)
             
             do r=1,LVT_rc%lnr
                do c=1,LVT_rc%lnc
                   if(LVT_stats%datamask(c,r).eq.1) then 
                      gid = LVT_domain%gindex(c,r)
                      if(gid.ne.-1) then 
                         if(value2d_ip(c,r).ne.LVT_rc%udef) then 
                            dataEntry%value(gid,k) = dataEntry%value(gid,k) + &
                                 value2d_ip(c,r)
                            dataEntry%count(gid,k) = dataEntry%count(gid,k) + 1
                         endif
                      endif
                   endif
                enddo
             enddo

             deallocate(value2d)
             deallocate(value2d_ip)

          elseif(trim(LVT_LIS_rc(source)%wopt).eq."1d gridspace") then 
             
             write(LVT_logunit,*) '[ERR] Transformation from 1d gridspace is not fully tested..'
             write(LVT_logunit,*) '[ERR] Stopping in LVT_readSingleBinaryVar '
             call LVT_endrun()
             
             allocate(value1d(LVT_LIS_rc(source)%ngrid))
             allocate(value2d_ip(LVT_rc%lnc, LVT_rc%lnr))
             
             call LVT_readvar_gridded(ftn,source,value1d, 1)
             
             call interp2lisgrid_1d(source, value1d,value2d_ip)

             do r=1,LVT_rc%lnr
                do c=1,LVT_rc%lnc
                   if(LVT_stats%datamask(c,r).eq.1) then 
                      gid = LVT_domain%gindex(c,r)
                      if(gid.ne.-1) then 
                         if(value2d_ip(c,r).ne.LVT_rc%udef) then 
                            dataEntry%value(gid,k) = dataEntry%value(gid,k) + &
                                 value2d_ip(c,r)
                            dataEntry%count(gid,k) = dataEntry%count(gid,k) + 1
                         endif
                      endif
                   endif
                enddo
             enddo

             deallocate(value1d)
             deallocate(value2d_ip)

          endif
       enddo
    endif

  end subroutine LVT_readSingleBinaryVar

  subroutine interp2lisgrid_2d(source, var_in, var_out)
    
    implicit none

    integer        :: source
    real           :: var_in(LVT_LIS_rc(source)%lnc,LVT_LIS_rc(source)%lnr)
    real           :: var_out(LVT_rc%lnc,LVT_rc%lnr)

    real           :: gi(LVT_LIS_rc(source)%lnc*LVT_LIS_rc(source)%lnr)
    logical*1      :: li(LVT_LIS_rc(source)%lnc*LVT_LIS_rc(source)%lnr)
    real           :: go(LVT_rc%lnc*LVT_rc%lnr)
    logical*1      :: lo(LVT_rc%lnc*LVT_rc%lnr)

    integer        :: c,r,t
    integer        :: mi
    integer        :: mo
    integer        :: iret

    mi = LVT_LIS_rc(source)%lnc*LVT_LIS_rc(source)%lnr
    mo = LVT_rc%lnc*LVT_rc%lnr

    if(LVT_LIS_rc(source)%gridDesc(1).eq.0) then 

       do r=1,LVT_LIS_rc(source)%lnr
          do c=1,LVT_LIS_rc(source)%lnc
             t = c+(r-1)*LVT_LIS_rc(source)%lnc
             gi(t) = var_in(c,r)
             if(gi(t).ne.LVT_rc%udef.and. &
                  gi(t).lt.1E+10) then 
                li(t) = .true. 
             else
                li(t) = .false. 
             endif
          enddo
       enddo

       if(LVT_LIS_rc(source)%gridDesc(10).gt.LVT_rc%gridDesc(10)) then
          !interpolate
          call bilinear_interp(LVT_rc%gridDesc,li,gi,&
               lo,go,mi,mo,LVT_rc%rlat_dn,LVT_rc%rlon_dn,&
               LVT_rc%w11_dn,LVT_rc%w12_dn,LVT_rc%w21_dn,LVT_rc%w22_dn,&
               LVT_rc%n11_dn,LVT_rc%n12_dn,LVT_rc%n21_dn,LVT_rc%n22_dn,&
               LVT_rc%udef,iret)

       else !upscale/average
          call upscaleByAveraging(mi,mo,LVT_rc%udef,&
               LVT_LIS_rc(source)%n11_up,li,gi,lo,go)
       endif

       do r=1,LVT_rc%lnr
          do c=1,LVT_rc%lnc
             t = c+(r-1)*LVT_rc%lnc
             if(lo(t)) then 
                var_out(c,r) = go(t)                
             else
                var_out(c,r) = LVT_rc%udef
             endif
          enddo
       enddo
    else 
!interpolation in other projections not supported currently
       var_out = var_in
    endif
  end subroutine interp2lisgrid_2d

!BOP
!
! !ROUTINE: applySpatialAveragingMask
! \label{applySpatialAveragingMask}
!
! !INTERFACE: 
  subroutine applySpatialAveragingMask(value2d)
! !ARGUMENTS:        
    real               :: value2d(LVT_rc%lnc,LVT_rc%lnr)
!
! !DESCRIPTION:
!   This subroutine computes the spatial average based on 
!   externally specified region mask. The region mask is 
!   expected to be a categorical map. 
! 
!EOP 
    integer            :: c,r,i,ind_val
    real,    allocatable   :: binval(:)
    integer, allocatable   :: nbinval(:)
    
    if(LVT_rc%sp_avg_mode.eq."region-based") then 
       allocate(binval(nint(LVT_rc%regmask_max)))
       allocate(nbinval(nint(LVT_rc%regmask_max)))
       
       binval = 0
       nbinval = 0 
       do r=1,LVT_rc%lnr
          do c=1,LVT_rc%lnc
             if(LVT_rc%regmask(c,r).ne.LVT_rc%udef.and.&
                  value2d(c,r).ne.LVT_rc%udef) then 
                ind_val = nint(LVT_rc%regmask(c,r))
                if(value2d(c,r).gt.0) then 
                   binval(ind_val) = binval(ind_val) + value2d(c,r)
                   nbinval(ind_val) = nbinval(ind_val) + 1
                endif
             endif
          enddo
       enddo
       do i=1,nint(LVT_rc%regmask_max)
          if(nbinval(i).ne.0) then 
             binval(i) = binval(i)/nbinval(i)
          endif
       enddo

       do r=1,LVT_rc%lnr
          do c=1,LVT_rc%lnc
             ind_val = nint(LVT_rc%regmask(c,r))             
             if(value2d(c,r).ne.LVT_rc%udef.and.&
                  ind_val.ne.LVT_rc%udef) then 
                value2d(c,r) = binval(ind_val)
             else
                value2d(c,r) = LVT_rc%udef
             endif
          enddo
       enddo
       
       deallocate(binval)
       deallocate(nbinval)
    endif
    
  end subroutine applySpatialAveragingMask
!BOP
! 
! !ROUTINE: interp2lisgrid_1d
! \label(interp2lisgrid_1d)
!
! !INTERFACE:
  subroutine interp2lisgrid_1d(source,var_in, var_out)
! 
! !USES:   

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
    integer        :: source
    real           :: var_in(LVT_rc%ngrid)
    real           :: var_out(LVT_rc%lnc,LVT_rc%lnr)

    real           :: gi(LVT_LIS_rc(source)%lnc*LVT_LIS_rc(source)%lnr)
    logical*1      :: li(LVT_LIS_rc(source)%lnc*LVT_LIS_rc(source)%lnr)
    real           :: go(LVT_rc%lnc*LVT_rc%lnr)
    logical*1      :: lo(LVT_rc%lnc*LVT_rc%lnr)

    integer        :: c,r,t
    integer        :: mi
    integer        :: mo
    integer        :: iret
    integer        :: gid

    mi = LVT_LIS_rc(source)%lnc*LVT_LIS_rc(source)%lnr
    mo = LVT_rc%lnc*LVT_rc%lnr

    if(LVT_LIS_rc(source)%gridDesc(1).eq.0) then 

       do r=1,LVT_LIS_rc(source)%lnr
          do c=1,LVT_LIS_rc(source)%lnc
             t = c+(r-1)*LVT_LIS_rc(source)%lnc
             gid = LVT_domain%gindex(c,r)
             if(gid.ne.-1) then 
                gi(t) = var_in(gid)
             else
                gi(t) = LVT_rc%udef
             endif
             if(gi(t).ne.LVT_rc%udef) then 
                li(t) = .true. 
             else
                li(t) = .false. 
             endif
          enddo
       enddo

       if(LVT_LIS_rc(source)%gridDesc(10).ge.LVT_rc%gridDesc(10)) then !interpolate

          call bilinear_interp(LVT_rc%gridDesc,li,gi,&
               lo,go,mi,mo,LVT_rc%rlat_dn,LVT_rc%rlon_dn,&
               LVT_rc%w11_dn,LVT_rc%w12_dn,LVT_rc%w21_dn,LVT_rc%w22_dn,&
               LVT_rc%n11_dn,LVT_rc%n12_dn,LVT_rc%n21_dn,LVT_rc%n22_dn,&
               LVT_rc%udef,iret)

       else !upscale/average
          call upscaleByAveraging(mi,mo,LVT_rc%udef,&
               LVT_LIS_rc(source)%n11_up,li,gi,lo,go)
       endif

       do r=1,LVT_rc%lnr
          do c=1,LVT_rc%lnc
             t = c+(r-1)*LVT_rc%lnc
             if(lo(t)) then 
                var_out(c,r) = go(t)
             else
                var_out(c,r) = LVT_rc%udef
             endif
          enddo
       enddo          
    endif

  end subroutine interp2lisgrid_1d

!BOP
! 
! !ROUTINE: LVT_readSingleNETCDFVar
! \label{LVT_readSingleNETCDFVar}
!
! !INTERFACE:
  subroutine LVT_readSingleNETCDFVar(ftn, source, dataEntry)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine reads a single variable from the LIS output file in the 
!   NETCDF format and interpolates/aggregates the data to the LVT analysis
!   grid. Based on the external datamask, the routine also filters
!   each variable. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    integer :: ftn
    integer :: source
    type(LVT_LISmetadataEntry) :: dataEntry
!EOP
    real, allocatable :: value2d(:,:)
    real, allocatable :: value2d_ip(:,:)
    real, allocatable :: value3d(:,:,:)
    real, allocatable :: value1d(:)
    real, allocatable :: value1d_mvars(:,:)

    integer :: unit_id
    integer :: ios
    integer :: t,m,k,i,c,r,gid, varid
    character*100 :: short_name

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    unit_id = -1
    if(dataEntry%selectOpt.eq.1) then 
       do i=1,dataEntry%nunits
          if(trim(dataEntry%units).eq.trim(dataEntry%unittypes(i))) then 
             unit_id = i
             exit
          endif
       enddo    
       if(unit_id.eq.-1) then 
          write(LVT_logunit,*) '[ERR] routine to diagnose ',trim(dataEntry%standard_name),&
               ' in units of ',trim(dataEntry%units),' is not defined'
          write(LVT_logunit,*) '[ERR] for diagnostic output...'
          call LVT_endrun()
       endif
       if(trim(LVT_LIS_rc(source)%wopt).eq."1d tilespace") then 
          if(dataEntry%vlevels.eq.1) then           
             if(dataEntry%timeAvgOpt.eq.0) then 
                short_name = trim(dataEntry%short_name)//'_inst'
             elseif(dataEntry%timeAvgOpt.eq.1) then 
                short_name = trim(dataEntry%short_name)//'_tavg'
             elseif(dataEntry%timeAvgOpt.eq.3) then 
                short_name = trim(dataEntry%short_name)//'_acc'
             else
                write(LVT_logunit,*) '[ERR] Please specify if the time averaged/instantaneous/accumulated variable '
                write(LVT_logunit,*) '[ERR] is to be read. Change the OUTPUT config'
                write(LVT_logunit,*) '[ERR] to specify the time averaging option for'
                write(LVT_logunit,*) '[ERR] ',trim(dataEntry%short_name)
                call LVT_endrun()
             endif
             ios = nf90_inq_varid(ftn,trim(short_name),varid)
!This is for backward compatibility support - remove later
!If the suffixes do not work, try the variable name without the suffix
! and stop the program if that also fail.

             if(ios.ne.0) then 
                short_name = trim(dataEntry%short_name)
                ios = nf90_inq_varid(ftn,trim(short_name),varid)
                call LVT_verify(ios,'Error in nf90_inq_varid:'//&
                     trim(dataEntry%standard_name))
             endif
             
             allocate(value1d(LVT_LIS_rc(source)%ntiles))
             allocate(value2d(LVT_LIS_rc(source)%lnc,LVT_LIS_rc(source)%lnr))
             allocate(value2d_ip(LVT_rc%lnc,LVT_rc%lnr))
             
             ios = nf90_get_var(ftn,varid,value1d)
             call LVT_verify(ios,'Error in nf90_get_var for '&
                  //trim(dataEntry%short_name))
             
             if(LVT_rc%computeEnsMetrics.eq.1) then 
                do t=1,LVT_LIS_rc(source)%ntiles
                   c = LVT_LIS_domain(source)%tile(t)%col
                   r = LVT_LIS_domain(source)%tile(t)%row
                   if(LVT_stats%datamask(c,r).eq.1) then 
                      gid = LVT_domain%gindex(c,r)
                      if(gid.ne.-1) then 
                         if(value1d(t).ne.LVT_rc%udef) then 
                            dataEntry%value(t,1) = &
                                 dataEntry%value(t,1) + &
                                 value1d(t)
                            dataEntry%count(t,1) = &
                                 dataEntry%count(t,1) + 1
                         endif
                      endif
                   endif
                enddo
             else
                value2d = 0
                do i=1,LVT_LIS_rc(source)%ntiles,LVT_LIS_rc(source)%nensem
                   do m=1,LVT_LIS_rc(source)%nensem
                      t = i+m-1
                      c = LVT_LIS_domain(source)%tile(t)%col
                      r = LVT_LIS_domain(source)%tile(t)%row
                      value2d(c,r) = value2d(c,r) + value1d(t)*&
                           LVT_LIS_domain(source)%tile(t)%fgrd/&
                           LVT_LIS_rc(source)%nensem
                   enddo
                enddo
                
                call interp2lisgrid_2d(source, value2d,value2d_ip)

                do r=1,LVT_rc%lnr
                   do c=1,LVT_rc%lnc
                      if(LVT_stats%datamask(c,r).eq.1) then 
                         gid = LVT_domain%gindex(c,r)
                         if(gid.ne.-1) then 
                            if(value2d_ip(c,r).ne.LVT_rc%udef) then 
                               dataEntry%value(gid,1) = &
                                    dataEntry%value(gid,1) + &
                                    value2d_ip(c,r)
                               dataEntry%count(gid,1) = &
                                    dataEntry%count(gid,1) + 1
                            endif
                         endif
                      endif
                   enddo
                enddo
             endif

             deallocate(value1d)
             deallocate(value2d)
             deallocate(value2d_ip)
          else
             if(dataEntry%timeAvgOpt.eq.0) then 
                short_name = trim(dataEntry%short_name)//'_inst'
             elseif(dataEntry%timeAvgOpt.eq.1) then 
                short_name = trim(dataEntry%short_name)//'_tavg'
             elseif(dataEntry%timeAvgOpt.eq.3) then 
                short_name = trim(dataEntry%short_name)//'_acc'
             else
                write(LVT_logunit,*) '[ERR] Please specify if the time averaged/instantaneous/accumulated variable '
                write(LVT_logunit,*) '[ERR] is to be read. Change the OUTPUT config'
                write(LVT_logunit,*) '[ERR] to specify the time averaging option for'
                write(LVT_logunit,*) '[ERR] ',trim(dataEntry%short_name)
                call LVT_endrun()
             endif

             ios = nf90_inq_varid(ftn,trim(short_name),varid)
!This is for backward compatibility support - remove later
!If the suffixes do not work, try the variable name without the suffix
! and stop the program if that also fail.

             if(ios.ne.0) then 
                short_name = trim(dataEntry%short_name)
                ios = nf90_inq_varid(ftn,trim(short_name),varid)
                call LVT_verify(ios,'Error in nf90_inq_varid:'//&
                     trim(dataEntry%standard_name))
             endif
             
             allocate(value1d_mvars(LVT_LIS_rc(source)%ntiles,dataEntry%vlevels))
             allocate(value1d(LVT_LIS_rc(source)%ntiles))
             allocate(value2d(LVT_LIS_rc(source)%lnc,LVT_LIS_rc(source)%lnr))
             allocate(value2d_ip(LVT_rc%lnc,LVT_rc%lnr))
             
             ios = nf90_get_var(ftn,varid,value1d_mvars)
             call LVT_verify(ios,'Error in nf90_get_var for'//&
                  trim(dataEntry%short_name))
             
             do k=1, dataEntry%vlevels
                value1d = value1d_mvars(:,k)
                if(LVT_rc%computeEnsMetrics.eq.1) then 
                   do t=1,LVT_LIS_rc(source)%ntiles
                      c = LVT_LIS_domain(source)%tile(t)%col
                      r = LVT_LIS_domain(source)%tile(t)%row
                      if(LVT_stats%datamask(c,r).eq.1) then 
                         gid = LVT_domain%gindex(c,r)
                         if(gid.ne.-1) then 
                            if(value1d(t).ne.LVT_rc%udef) then 
                               dataEntry%value(t,k) = &
                                    dataEntry%value(t,k) + &
                                    value1d(t)
                               dataEntry%count(t,k) = &
                                    dataEntry%count(t,k) + 1
                            endif
                         endif
                      endif
                   enddo
                else
                   value2d = 0
                   do i=1,LVT_LIS_rc(source)%ntiles,LVT_LIS_rc(source)%nensem
                      do m=1,LVT_LIS_rc(source)%nensem
                         t = i+m-1
                         c = LVT_LIS_domain(source)%tile(t)%col
                         r = LVT_LIS_domain(source)%tile(t)%row
                         value2d(c,r) = value2d(c,r) + value1d(t)*&
                              LVT_LIS_domain(source)%tile(t)%fgrd/&
                              LVT_LIS_rc(source)%nensem
                      enddo
                   enddo
                   
                   call interp2lisgrid_2d(source,value2d,value2d_ip)
                
                   do r=1,LVT_rc%lnr
                      do c=1,LVT_rc%lnc
                         if(LVT_stats%datamask(c,r).eq.1) then 
                            gid = LVT_domain%gindex(c,r)
                            if(gid.ne.-1) then 
                               if(value2d_ip(c,r).ne.LVT_rc%udef) then 
                                  dataEntry%value(gid,k) = &
                                       dataEntry%value(gid,k) + &
                                       value2d_ip(c,r)
                                  dataEntry%count(gid,k) = &
                                       dataEntry%count(gid,k) + 1
                               endif
                            endif
                         endif
                      enddo
                   enddo
                   
                endif
             enddo
          
             deallocate(value1d_mvars)
             deallocate(value1d)
             deallocate(value2d)
             deallocate(value2d_ip)
          endif
          
       elseif(trim(LVT_LIS_rc(source)%wopt).eq."2d gridspace") then 
          if(dataEntry%vlevels.eq.1) then 
             if(dataEntry%timeAvgOpt.eq.0) then 
                short_name = trim(dataEntry%short_name)//'_inst'
             elseif(dataEntry%timeAvgOpt.eq.1) then 
                short_name = trim(dataEntry%short_name)//'_tavg'
             elseif(dataEntry%timeAvgOpt.eq.3) then 
                short_name = trim(dataEntry%short_name)//'_acc'
             else
                write(LVT_logunit,*) '[ERR] Please specify if the time averaged/instantaneous/accumulated variable '
                write(LVT_logunit,*) '[ERR] is to be read. Change the OUTPUT config'
                write(LVT_logunit,*) '[ERR] to specify the time averaging option for'
                write(LVT_logunit,*) '[ERR] ',trim(dataEntry%short_name)
                call LVT_endrun()
             endif
!to switch between longnames and shortnames
!             ios = nf90_inq_varid(ftn,trim(dataEntry%standard_name),varid)
             ios = nf90_inq_varid(ftn,trim(short_name),varid)
!This is for backward compatibility support - remove later
!If the suffixes do not work, try the variable name without the suffix
! and stop the program if that also fail.
             if(ios.ne.0) then 
                short_name = trim(dataEntry%short_name)
                ios = nf90_inq_varid(ftn,trim(short_name),varid)
                call LVT_verify(ios,'Error in nf90_inq_varid:'//&
                     trim(short_name))
             endif

             allocate(value2d(LVT_LIS_rc(source)%lnc, LVT_LIS_rc(source)%lnr))
             allocate(value2d_ip(LVT_rc%lnc, LVT_rc%lnr))

             value2d = LVT_rc%udef

             ios = nf90_get_var(ftn,varid,value2d)
             call LVT_verify(ios,'Error in nf90_get_var for'//&
                  trim(dataEntry%short_name))
             call interp2lisgrid_2d(source,value2d,value2d_ip)

             call applySpatialAveragingMask(value2d_ip)

             do r=1,LVT_rc%lnr
                do c=1,LVT_rc%lnc
                   if(LVT_stats%datamask(c,r).eq.1) then 
                      gid = LVT_domain%gindex(c,r)
                      if(gid.ne.-1) then 
                         if(value2d_ip(c,r).ne.LVT_rc%udef) then 
                            dataEntry%value(gid,1) = dataEntry%value(gid,1) + &
                                 value2d_ip(c,r)
                            dataEntry%count(gid,1) = dataEntry%count(gid,1) + 1
                         endif
                      endif
                   endif
                enddo
             enddo

             deallocate(value2d)
             deallocate(value2d_ip)
          else
!to switch between longnames and shortnames
!             ios = nf90_inq_varid(ftn,trim(dataEntry%standard_name),varid)
             if(dataEntry%timeAvgOpt.eq.0) then 
                short_name = trim(dataEntry%short_name)//'_inst'
             elseif(dataEntry%timeAvgOpt.eq.1) then 
                short_name = trim(dataEntry%short_name)//'_tavg'
             elseif(dataEntry%timeAvgOpt.eq.3) then 
                short_name = trim(dataEntry%short_name)//'_acc'
             else
                write(LVT_logunit,*) '[ERR] Please specify if the time averaged/instantaneous/accumulated variable '
                write(LVT_logunit,*) '[ERR] is to be read. Change the OUTPUT config'
                write(LVT_logunit,*) '[ERR] to specify the time averaging option for'
                write(LVT_logunit,*) '[ERR] ',trim(dataEntry%short_name)
                call LVT_endrun()
             endif
             ios = nf90_inq_varid(ftn,trim(short_name),varid)
!This is for backward compatibility support - remove later
!If the suffixes do not work, try the variable name without the suffix
! and stop the program if that also fail.
             if(ios.ne.0) then 
                short_name = trim(dataEntry%short_name)
                ios = nf90_inq_varid(ftn,trim(short_name),varid)
                call LVT_verify(ios,'Error in nf90_inq_varid:'//&
                     trim(dataEntry%short_name))
             endif

             allocate(value3d(LVT_LIS_rc(source)%lnc, &
                  LVT_LIS_rc(source)%lnr,dataEntry%vlevels))
             allocate(value2d_ip(LVT_rc%lnc, LVT_rc%lnr))

             ios = nf90_get_var(ftn,varid,value3d)
             call LVT_verify(ios,'Error in nf90_get_var for'//&
                  trim(dataEntry%short_name))

             do k=1, dataEntry%vlevels
                call interp2lisgrid_2d(source,value3d(:,:,k),value2d_ip)
               
                do r=1,LVT_rc%lnr
                   do c=1,LVT_rc%lnc
                      if(LVT_stats%datamask(c,r).eq.1) then 
                         gid = LVT_domain%gindex(c,r)
                         if(gid.ne.-1) then 
                            if(value2d_ip(c,r).ne.LVT_rc%udef) then 
                               dataEntry%value(gid,k) = &
                                    dataEntry%value(gid,k) + &
                                    value2d_ip(c,r)
                               dataEntry%count(gid,k) = &
                                    dataEntry%count(gid,k) + 1
                            endif
                         endif
                      endif
                   enddo
                enddo
             enddo

             deallocate(value3d)
             deallocate(value2d_ip)
          endif

       elseif(trim(LVT_LIS_rc(source)%wopt).eq."1d gridspace") then 
          
          print*, 'Reading in land vector format is not yet supported'
          print*, 'program stopping ...'
          stop

       endif
    endif
! The stratification variable is chosen from the first data stream
    if(LVT_rc%var_based_strat .gt. 0) then
       if(source.eq.1.and.dataEntry%selectOpt.eq.1) then 
          if(dataEntry%short_name.eq.LVT_rc%vname_strat) then 
             LVT_histData%strat_varEntry%value = dataEntry%value
             LVT_histData%strat_varEntry%count = dataEntry%count
          endif
       endif
    endif
    

#endif

  end subroutine LVT_readSingleNETCDFVar


!BOP
! !ROUTINE: mapLISdataToLVT
! \label{mapLISdatatoLVT}
!
! !INTERFACE: 
  subroutine mapLISdataToLVT(source,lvtdataEntry, lisdataEntry)
! !ARGUMENTS: 

    integer                    :: source
    type(LVT_metadataEntry)    :: lvtdataEntry
    type(LVT_LISmetadataEntry) :: lisdataEntry
    
!
! !DESCRIPTION: 
! 
! This routine maps the LIS output data to the LVT data
! structure for a particular variable
! 
!  The arguments are: 
!  \begin{description}
!   \item [lvtdataEntry]
!     the LVT data entry object
!   \item [lisdataEntry]
!     the LIS data entry object
!  \end{description}
!EOP
    integer             :: k 
    real                :: scale_f
!    do k=1, dataEntry%selectNlevs
    scale_f = 1.0

    do k=1, lvtdataEntry%vlevels
       if(lvtdataEntry%units.eq.lisdataEntry%units) then 
          scale_f = 1.0
       elseif(lvtdataEntry%units.eq."m3/m3".and.lisdataEntry%units.eq."kg/m2") then 
          if(LVT_LIS_rc(source)%model_name.eq."Noah.3.3") then 
             scale_f = 1/(LVT_CONST_RHOFW*LVT_LIS_rc(source)%smthick(k))
          endif
       elseif(lvtdataEntry%units.eq."kg/m2s".and.lisdataEntry%units.eq."W/m2") then 
          scale_f = 1/LVT_CONST_LATVAP
       elseif(lvtdataEntry%units.eq."W/m2".and.lisdataEntry%units.eq."kg/m2s") then 
          scale_f = LVT_CONST_LATVAP
       else
          write(LVT_logunit,*) '[ERR] '
          write(LVT_logunit,*) '[ERR] The units of the LIS output and the analysis variable'
          write(LVT_logunit,*) '[ERR] are mismatched'
          call LVT_endrun()
       endif
       lvtdataEntry%value(:,k) = lisdataEntry%value(:,k)*scale_f
       lvtdataEntry%count(:,k) = lisdataEntry%count(:,k)

    enddo
  end subroutine mapLISdataToLVT
!BOP
! 
! !ROUTINE: LVT_readSingleGrib1Var
! \label{LVT_readSingleGrib1Var}
!
! !INTERFACE:
  subroutine LVT_readSingleGrib1Var(ftn, source, dataEntry,nvars,gvar,&
       undef_v,pid,tid,lid)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine reads a single variable from the LIS output file in the 
!   grib1 format. Based on the external datamask, the routine also filters
!   each variable. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    integer :: ftn
    integer :: source
    type(LVT_LISmetadataEntry) :: dataEntry
    integer :: nvars
    real    :: undef_v
    integer :: pid(nvars)
    integer :: tid(nvars)
    integer :: lid(nvars)
    real    :: gvar(nvars,LVT_LIS_rc(source)%lnc*LVT_LIS_rc(source)%lnr)
    integer :: local_tid
! 
!
!EOP
    real, allocatable :: value1d(:)
    real, allocatable :: value2d(:,:)
    real, allocatable :: value2d_ip(:,:)
    logical*1, allocatable :: lo(:)
    integer :: unit_id
    integer :: k,i,c,r,t,gid,kk,index
    integer :: j,lubi, nsize, iret
    integer :: jpds(200), jgds(200)
    integer :: kf,kpds(200),gridDesc(200)
    
    unit_id = -1
    if(dataEntry%selectOpt.eq.1) then 
       do i=1,dataEntry%nunits
          if(trim(dataEntry%units).eq.trim(dataEntry%unittypes(i))) then 
             unit_id = i
             exit
          endif
       enddo    
       if(unit_id.eq.-1) then 
          write(LVT_logunit,*) '[ERR] routine to diagnose ',trim(dataEntry%standard_name),&
               ' in units of ',trim(dataEntry%units),' is not defined'
          write(LVT_logunit,*) '[ERR] for diagnostic output...'
          call LVT_endrun()
       endif

       nsize = LVT_LIS_rc(source)%lnc*LVT_LIS_rc(source)%lnr
       allocate(value1d(nsize))
       allocate(value2d(LVT_LIS_rc(source)%lnc, LVT_LIS_rc(source)%lnr))
       allocate(value2d_ip(LVT_rc%lnc, LVT_rc%lnr))
       allocate(lo(nsize))

       iret = 1

#if (defined AFWA_GRIB_CONFIGS)
       if(dataEntry%timeAvgOpt.eq.0) then !instantaneous only 
          local_tid = 1
       elseif(dataEntry%timeAvgOpt.eq.1) then !time averaged only 
          local_tid = 7
       elseif(dataEntry%timeAvgOpt.eq.3) then !accumulated
          local_tid = 133
       else
          local_tid = 1 !read only the time averaged field
       endif
#else
       if(dataEntry%timeAvgOpt.eq.0) then !instantaneous only 
          local_tid = 0
       elseif(dataEntry%timeAvgOpt.eq.1) then !time averaged only 
          local_tid = 7
       elseif(dataEntry%timeAvgOpt.eq.3) then !accumulated
          local_tid = 4
       else
          local_tid = 1 !read only the time averaged field
       endif
#endif

       do k=1,dataEntry%vlevels 
          do index=1,nvars
             if(pid(index).eq.dataEntry%varid_def.and.&
                  tid(index).eq.local_tid) then 
                if(dataEntry%vlevels.gt.1.and.&
                     trim(dataEntry%short_name).eq."SoilMoist") then 
                   if(lid(index).eq.LVT_LIS_rc(source)%smdepth(k)*100) then !covert to cm
                      value1d = gvar(index,:)
                   endif
                elseif(dataEntry%vlevels.gt.1.and.&
                     trim(dataEntry%short_name).eq."SoilTemp") then 
                   if(lid(index).eq.LVT_LIS_rc(source)%stdepth(k)*100) then !covert to cm
                      do t=1,nsize
                         if(gvar(index,t).lt.400) then
                            value1d(t) = gvar(index,t)
                         else
                            value1d(t) = LVT_rc%udef
                         endif
                      enddo
                   endif
                else
                   value1d = gvar(index,:)
                endif
                iret = 0
             endif
          enddo
          if(iret.ne.0) then 
             write(LVT_logunit,*) '[ERR] Reading variable ',trim(dataEntry%standard_name), ' failed'
             call LVT_endrun()
          endif

          do r=1,LVT_LIS_rc(source)%lnr
             do c=1,LVT_LIS_rc(source)%lnc
                if(value1d(c+(r-1)*LVT_LIS_rc(source)%lnc).ne.undef_v) then 
                   value2d(c,r) = value1d(c+(r-1)*LVT_LIS_rc(source)%lnc)
                else
                   value2d(c,r) = LVT_rc%udef
                endif
             enddo
          enddo
          
          call interp2lisgrid_2d(source,value2d, value2d_ip)

          do r=1,LVT_rc%lnr
             do c=1,LVT_rc%lnc
                if(LVT_stats%datamask(c,r).eq.1) then 
                   gid = LVT_domain%gindex(c,r)
                   if(gid.ne.-1) then 
                      if(value2d_ip(c,r).ne.LVT_rc%udef) then
                         dataEntry%value(gid,k) = &
                              dataEntry%value(gid,k) + &
                              value2d_ip(c,r)
                         dataEntry%count(gid,k) = &
                              dataEntry%count(gid,k) + 1
                      endif
                   endif
                endif
             enddo
          enddo          
       enddo

       deallocate(value1d)
       deallocate(value2d)
       deallocate(value2d_ip)
       deallocate(lo)
    endif
  end subroutine LVT_readSingleGrib1Var


!BOP
! 
! !ROUTINE: readvar_1dgridded_fromvector_real
! \label{readvar_1dgridded_fromvector_real}
!
! !INTERFACE:
! Private name: call LVT_readvar_gridded
  subroutine readvar_1dgridded_fromvector_real(ftn, source,var, oned)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in)   :: ftn
    integer, intent(in)   :: source
    real, intent(inout)   :: var(LVT_LIS_rc(source)%ngrid)
    integer, intent(in)   :: oned
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Reads a real variable from a binary 
!  sequential access, 1d gridded file. After reading
!  the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!   \item [oned]
!     dummy variable to distinguish the interface. 
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    real, allocatable :: gtmp(:)
    real, allocatable :: gtmp2d(:,:)
    integer :: i, gid
    integer :: c,r, nc, cnt
    integer :: c1,c2,r1,r2

    allocate(gtmp(LVT_LIS_rc(source)%glbngrid))
    read(ftn) gtmp

    allocate(gtmp2d(LVT_LIS_rc(source)%gnc, LVT_LIS_rc(source)%gnr))

    cnt = 1
    do r=1,LVT_LIS_rc(source)%gnr
       do c=1,LVT_LIS_rc(source)%gnc
          gid = c+(r-1)*LVT_LIS_rc(source)%gnc
          if(LVT_LIS_domain(source)%ntiles_pergrid(gid).gt.0) then 
             gtmp2d(c,r) = gtmp(cnt)
             cnt = cnt+1
          endif
       enddo
    enddo

    nc = (LVT_lis_ewe_halo_ind(source,LVT_localPet+1)-LVT_lis_ews_halo_ind(source,LVT_localPet+1))+1

    do r=LVT_lis_nss_halo_ind(source,LVT_localPet+1),LVT_lis_nse_halo_ind(source,LVT_localPet+1)
       do c=LVT_lis_ews_halo_ind(source,LVT_localPet+1),LVT_lis_ewe_halo_ind(source,LVT_localPet+1)
          c1 = c-LVT_lis_ews_halo_ind(source,LVT_localPet+1)+1
          r1 = r-LVT_lis_nss_halo_ind(source,LVT_localPet+1)+1
!          gid = r1+(c1-1)*nc
          gid = LVT_LIS_domain(source)%gindex(c1,r1)
          if(gid.ne.-1) then
             var(gid) = gtmp2d(c,r)
          endif
       enddo
    enddo

    deallocate(gtmp)
    deallocate(gtmp2d)


  end subroutine readvar_1dgridded_fromvector_real

!BOP
! 
! !ROUTINE: readvar_2dgridded_real
! \label{readvar_2dgridded_real}
!
! !INTERFACE:
! Private name: call LVT_readvar_gridded
  subroutine readvar_2dgridded_real(ftn,  source, var)
! 
! !USES:

    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in)   :: ftn
    integer, intent(in)   :: source
    real, intent(inout)   :: var(LVT_LIS_rc(source)%lnc, LVT_LIS_rc(source)%lnr)
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Reads a real variable from a binary 
!  sequential access, gridded file. After reading
!  the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    real, allocatable :: gtmp(:,:)
    integer :: c,r
    integer :: nc,c1,r1

    allocate(gtmp(LVT_LIS_rc(source)%gnc,LVT_LIS_rc(source)%gnr))
    read(ftn) gtmp
    
    nc = (LVT_lis_ewe_halo_ind(source,LVT_localPet+1)-LVT_lis_ews_halo_ind(source,LVT_localPet+1))+1

    do r=LVT_lis_nss_halo_ind(source,LVT_localPet+1),LVT_lis_nse_halo_ind(source,LVT_localPet+1)
       do c=LVT_lis_ews_halo_ind(source,LVT_localPet+1),LVT_lis_ewe_halo_ind(source,LVT_localPet+1)
          c1 = c-LVT_lis_ews_halo_ind(source,LVT_localPet+1)+1
          r1 = r-LVT_lis_nss_halo_ind(source,LVT_localPet+1)+1
          var(c1,r1) = gtmp(c,r)
       enddo
    enddo
    deallocate(gtmp)

  end subroutine readvar_2dgridded_real

!BOP
! 
! !ROUTINE: LVT_resetLISoutputContainers
! \label{LVT_resetLISoutputContainers}
!
! !INTERFACE: 
  subroutine LVT_resetLISoutputContainers(source)
! 
! !USES: 
!
! !INPUT PARAMETERS: 
    integer,   intent(in)   :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine reads a NETCDF output file based on the list of selected 
!  output variables. 
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest \newline
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LVT\_readSingleNETCDFVar](\ref{LVT_readSingleNETCDFVar}) \newline
!     reads a single variable from a flat binary file. 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    type(LVT_LISmetadataEntry), pointer :: lisdataEntry

    ! Read the LIS outputContainers
    if(LVT_LIS_rc(source)%anlys_data_class.eq."LSM") then 
       lisdataEntry => LVT_LISoutput(source)%head_lsm_list
    elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Routing") then 
       lisdataEntry => LVT_LISoutput(source)%head_routing_list
    elseif(LVT_LIS_rc(source)%anlys_data_class.eq."RTM") then 
       lisdataEntry => LVT_LISoutput(source)%head_rtm_list
    elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Irrigation") then 
       lisdataEntry => LVT_LISoutput(source)%head_irrig_list
    endif

    do while(associated(lisdataEntry))
       call resetSingleLISoutputContainer(source, lisdataEntry)
       lisdataEntry => lisdataEntry%next
    enddo
  end subroutine LVT_resetLISoutputContainers

!BOP
! 
! !ROUTINE: resetSingleLISoutputContainer
! \label{resetSingleLISoutputContainer}
!
! !INTERFACE:
  subroutine resetSingleLISoutputContainer(source, dataEntry)
! 
! !USES:   
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
!BOP
! !ARGUMENTS: 
    integer :: source
    type(LVT_LISmetadataEntry) :: dataEntry
    integer  :: k 
!EOP
    if(dataEntry%selectOpt.eq.1) then 
       do k=1,dataEntry%vlevels
          dataEntry%value(:,k) = 0.0
          dataEntry%count(:,k) = 0
       enddo
    endif

  end subroutine resetSingleLISoutputContainer


end module LVT_LISoutputHandlerMod
