!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: LVT_PRIV_rcMod
! \label(LVT_PRIV_rcMod)
!
! !INTERFACE:
module LVT_PRIV_rcMod 
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  
!  Module for specifying model independent variables in LIS 
!  verification toolkit. 
!!
! The variables specified in this module include: 
!
! \begin{description}
!  \item[wout]
!   output format option used in LIS (tiled/gridded/grib/netcdf/grib2) 
!  \item[pass]
!   Number of passes through the time series required for the computation
!   of the metric (e.g. std, anomaly correlations)
!  \item[wtsout]
!   Flag to check whether to output time series data
!  \item[obsCountThreshold]
!   Minimum number of observations (in time), used
!   as a threshold to compute statistics
!  \item[scCountThreshold]
!   Minimum number of points (in time), used
!   as a threshold to compute average seasonal cycles
!  \item[adcCountThreshold]
!   Minimum number of points (in time), used
!   as a threshold to compute average diurnal cycles.
! \end{description}
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY:
!  02 Oct 2008; Sujay Kumar; Initial Specification
!  03 Dec 2012; Shugong Wang; Add lis_version for backward support of LIS 6 
!
!EOP

  implicit none
  type lvtrcdec
     character*100          :: configfile
     integer                :: max_model_types

     character*50           :: runmode
     character*50           :: domain
     integer                :: nnest
     ! The following line is uncommented by Shugong Wang for backword support. 12/04/2012
     character*50           :: lsm
!     character*50           :: routing_model 
!     character*50           :: rtm_model 
     integer                :: nt

     integer                :: metric_sindex
     integer                :: metric_eindex

     integer                :: lsm_index 
     integer                :: lake_index 
     integer                :: glacier_index 

     integer                :: npts
     integer                :: ntiles
     integer                :: glbntiles
     integer                :: glbntiles_red
     integer                :: ngrid 
     integer                :: glbngrid
     integer                :: glbngrid_red


     integer                :: gnc
     integer                :: gnr
     integer                :: lnc
     integer                :: lnr
     integer                :: lnc_red
     integer                :: lnr_red
     integer,  allocatable      :: lnc_sc(:)
     integer,  allocatable      :: lnr_sc(:)

     integer                :: pnc  
     integer                :: pnr
     integer                :: input_lnc
     integer                :: input_lnr

     integer                :: nparam

     character*50           :: lvt_out_format
     character*50           :: lvt_wopt
     character*20           :: startmode
     character*100          :: odir
!     character*20           :: model_name
     character*50           :: diagfile

     integer                :: npesx
     integer                :: npesy
     integer                :: halox
     integer                :: haloy

     integer                :: vegsrc

     character*100          :: paramfile


     integer                :: vfile_form
     real                   :: gridDesc(50)
     real                   :: input_gridDesc(50)

     integer                :: sdoy        
     integer                :: sss         
     integer                :: smn         
     integer                :: shr         
     integer                :: sda         
     integer                :: smo         
     integer                :: syr         
     integer                :: endcode        
     integer                :: ess            
     integer                :: emn            
     integer                :: edoy           
     integer                :: ehr            
     integer                :: eda            
     integer                :: emo            
     integer                :: eyr            
     integer                :: endtime        
     real*8                 :: etime               
     real                   :: egmt
     integer                :: doy
     integer                :: yr
     integer                :: mo
     integer                :: prev_mo_rst
     integer                :: prev_mo_tavg
     integer                :: prev_yr_tavg
     integer                :: prev_yr_sout
     integer                :: prev_mo_sout
     integer                :: use_shift_mo
     integer                :: timeAvgOpt
     integer                :: da
     integer                :: hr
     integer                :: mn
     integer                :: ss 
     
     integer                :: nyr
     integer                :: nmo
     integer                :: nda
     integer                :: nhr
     integer                :: nmn
     integer                :: nss

     integer                :: dyr(2)
     integer                :: dmo(2)
     integer                :: dda(2)
     integer                :: dhr(2)
     integer                :: dmn(2)
     integer                :: dss(2)
     integer                :: ddoy(2)
     real*8                 :: dtime(2)      
     real                   :: dgmt(2)

     integer                :: d_nyr(2)
     integer                :: d_nmo(2)
     integer                :: d_nda(2)
     integer                :: d_nhr(2)
     integer                :: d_nmn(2)
     integer                :: d_nss(2)

     integer                :: lyr
     integer                :: lmo
     integer                :: lda
     integer                :: lhr
     integer                :: lmn
     integer                :: lss
     integer                :: ldoy
     real*8                 :: ltime
     real                   :: lgmt

     integer                :: l_nyr
     integer                :: l_nmo
     integer                :: l_nda
     integer                :: l_nhr
     integer                :: l_nmn
     integer                :: l_nss

     real*8                 :: time      
     real                   :: gmt
     
     real*8                 :: etime1
     real                   :: egmt1
     integer                :: edoy1
     integer                :: eyr1
     integer                :: emo1
     integer                :: eda1
     integer                :: ehr1
     integer                :: emn1
     integer                :: ess1
     integer                :: tscount
     integer                :: daycount
     integer                :: daycount_sout
     integer                :: monthCount
     integer                :: monthCount_sout
     character*100          :: rstfile

     integer                :: ts
     character*50           :: tsconv
     integer                :: nts

     integer                :: wtsout
     integer                :: extractts
     integer                :: tavgInterval
     integer                :: tlag

     character*200          :: outputSpecFile
     character*200          :: statsSpecFile
     character*50           :: statsodir
     integer                :: statswriteint
     integer                :: restartInterval
     integer                :: nstvars

     integer                :: nsmlayers
     integer                :: nstlayers
     real                   :: lis_sf_d
     real                   :: lis_rz_d
     
     character*20,  allocatable :: lisvarname(:)
     character*20,  allocatable :: lisvarunit(:)
     character*20,  allocatable :: obsvarname(:,:)
     character*20,  allocatable :: obsvarunit(:,:)

     character*50           :: obssource(2)

     integer                :: smoothObs
     integer                :: pass
     integer                :: curr_pass

     integer                :: dataMask
     integer                :: maskflag
     character*40           :: maskdir
     integer                :: monthly_mask(12)
     integer                :: ftn_summ_file
     integer                :: obsCountThreshold
     integer                :: computeInnovDist
     integer                :: computeGain
     integer                :: ComputeErrSC
     integer                :: computeADC
     integer                :: scCountThreshold
     integer                :: adcCountThreshold
     integer                :: scInterval
     integer                :: nasc
     integer                :: nadc
     character*10, allocatable  :: scname(:)
     character*10, allocatable  :: adcname(:)
     logical                :: obs_duplicate
     logical                :: computeFlag
     integer                :: computeEnsMetrics
     integer                :: computeICmetrics
     integer                :: ensLLType

     real                   :: pval_ci
     integer                :: var_based_strat
     integer                :: var_strat_index
     character*20           :: vname_strat
     real                   :: strat_var_threshold
     integer                :: strat_nlevels

     character*40           :: data_strat_attrib_file
     integer                :: data_based_strat
     character*40, allocatable  :: data_based_strat_file(:)
     character*40, allocatable  :: data_based_strat_var(:)
     integer                :: data_based_nstrats
     real,         allocatable  :: data_based_strat_max(:)
     real,         allocatable  :: data_based_strat_min(:)
     real,         allocatable  :: data_based_strat_delta(:)
     integer,      allocatable  :: data_based_strat_nbins(:)
     real,         allocatable  :: strat_data(:,:)

     integer                :: ntslocs
     character*200          :: tsspecfile
     integer                :: tsspecstyle

     integer                :: n_sc_locs
     character*200          :: sc_specfile
     integer                :: sc_specstyle

     integer                :: n_adc_locs
     character*200          :: adc_specfile
     integer                :: adc_specstyle

     real                   :: udef
     character*20           :: security_class
     character*20           :: distribution_class
     character*20           :: data_category
     character*20           :: area_of_data
     character*100          :: institution = 'NASA GSFC'

     integer                :: nscales
     logical                :: chkTS

     integer                :: anomalyTlength
     character*100          :: sp_avg_mode         
     character*100          :: reg_maskfile
     real,         allocatable  :: regmask(:,:)
     real                   :: regmask_max
     ! The following lines are added by Shugong Wang for backward support of LIS 6 
     integer                :: lis_version
     character*3            :: expcode
     logical                :: lis_output_obs
     
     real,     allocatable      :: rlat_dn(:)
     real,     allocatable      :: rlon_dn(:)
     real,     allocatable      :: w11_dn(:)
     real,     allocatable      :: w12_dn(:)
     real,     allocatable      :: w21_dn(:)
     real,     allocatable      :: w22_dn(:)
     integer,  allocatable      :: n11_dn(:)
     integer,  allocatable      :: n12_dn(:)
     integer,  allocatable      :: n21_dn(:)
     integer,  allocatable      :: n22_dn(:)

     integer,  allocatable      :: n11_up(:)

     logical                    :: ds1_dup
     logical                    :: ds2_dup
     
     character*100              :: trainingAlg

  end type lvtrcdec
  
  type lisrcdec
     integer                :: ts     
     character*50           :: anlys_data_class
     integer                :: nsf_model_types
     integer, allocatable       :: sf_model_type(:)
     character*50, allocatable  :: sf_model_type_name(:)
     integer, allocatable       :: sf_model_type_select(:)
     character*50, allocatable  :: sf_model_type_name_select(:)

     integer                :: nsurfacetypes
     character*100          :: domfile
     character*100          :: map_proj
     character*100          :: odir
     integer                :: nest
     character*100          :: style
     character*100          :: format

     integer                :: bareclass 
     integer                :: urbanclass
     integer                :: snowclass 
     integer                :: waterclass
     integer                :: wetlandclass
     integer                :: glacierclass
     integer                :: nvegtypes
     integer                :: nelevbands
     integer                :: nslopebands
     integer                :: naspectbands
     integer                :: nsoiltypes
     integer                :: nsoilfbands

     integer                :: ntiles
     integer                :: glbntiles
     integer                :: glbntiles_red
     integer                :: ngrid
     integer                :: glbngrid
     integer                :: glbngrid_red
     integer, allocatable   :: npatch(:)
     integer, allocatable   :: glbnpatch(:)
     integer, allocatable   :: glbnpatch_red(:)

     integer                :: lnc
     integer                :: lnr
     integer                :: lnc_red
     integer                :: lnr_red
     integer                :: gnc
     integer                :: gnr
     integer                :: nensem
     character*50           :: wopt
     character*100          :: useelevationmap
     character*100          :: useslopemap
     character*100          :: useaspectmap
     character*100          :: usetexturemap
     character*100          :: usesoilfractionmap
     real                   :: gridDesc(50)

     integer                :: surface_maxt
     real                   :: surface_minp    
     integer                :: soilt_maxt
     real                   :: soilt_minp    
     integer                :: soilf_maxt
     real                   :: soilf_minp    
     integer                :: elev_maxt
     real                   :: elev_minp    
     integer                :: slope_maxt
     real                   :: slope_minp    
     integer                :: aspect_maxt
     real                   :: aspect_minp    
     character*100          :: outputSpecFile

     real,     allocatable      :: rlat_dn(:)
     real,     allocatable      :: rlon_dn(:)
     real,     allocatable      :: w11_dn(:)
     real,     allocatable      :: w12_dn(:)
     real,     allocatable      :: w21_dn(:)
     real,     allocatable      :: w22_dn(:)
     integer,  allocatable      :: n11_dn(:)
     integer,  allocatable      :: n12_dn(:)
     integer,  allocatable      :: n21_dn(:)
     integer,  allocatable      :: n22_dn(:)

     integer,  allocatable      :: n11_up(:)
     character*50               :: model_name
     integer                :: nsmlayers
     integer                :: nstlayers
     real,          allocatable :: smthick(:)
     real,          allocatable :: stthick(:)
     real,          allocatable :: smdepth(:)
     real,          allocatable :: stdepth(:)

  end type lisrcdec
end module LVT_PRIV_rcMod
