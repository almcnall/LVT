!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readGIMMS_NDVIObs
! \label{readGIMMS_NDVIObs}
!
! !INTERFACE: 
subroutine readGIMMS_NDVIObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use GIMMS_NDVIobsMod
          
  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  NOTES: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  7 Mar 2015: Sujay Kumar, Initial Specification
! 
!EOP

  integer                :: c,r, tindex
  integer                :: flag
  integer                :: ftn
  character*100          :: fname_a, fname_b
  logical*1              :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: ndvi_out_a(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: ndvi_out_b(LVT_rc%lnc*LVT_rc%lnr)
  logical*1              :: input_bitmap(&
       GIMMSndviobs(source)%nc*GIMMSndviobs(source)%nr)
  real                   :: ndvi_inp(&
       GIMMSndviobs(source)%nc*GIMMSndviobs(source)%nr)
  real                   :: ndvi_inp1(&
       GIMMSndviobs(source)%nc, GIMMSndviobs(source)%nr)
  integer*2              :: ndvi_temp(&
       GIMMSndviobs(source)%nr,GIMMSndviobs(source)%nc)
  real                   :: varfield(LVT_rc%lnc,LVT_rc%lnr)
  logical                :: filea_exists, fileb_exists
  integer                :: iret
  integer                :: nc, nr
  
  varfield = LVT_rc%udef

  nc = GIMMSndviobs(source)%nc
  nr = GIMMSndviobs(source)%nr

  if(GIMMSndviobs(source)%mo.ne.LVT_rc%d_nmo(source)) then 

     if(GIMMSndviobs(source)%startFlag) then 
        GIMMSndviobs(source)%startFlag = .false. 
     endif
     GIMMSndviobs(source)%yr = LVT_rc%d_nyr(source)
     GIMMSndviobs(source)%mo = LVT_rc%d_nmo(source)

     call create_GIMMSndvi_filename(Gimmsndviobs(Source)%odir, &
          LVT_rc%dyr(source),LVT_rc%dmo(source),&
          fname_a, fname_b)
     
     inquire(file=trim(fname_a),exist=filea_exists) 
     inquire(file=trim(fname_b),exist=fileb_exists) 

     if(filea_exists.and.fileb_exists) then 
        write(LVT_logunit,*) '[INFO] Reading GIMMS NDVI file ',trim(fname_a)
            
        ftn = LVT_getNextUnitNumber()
        open(ftn,file=trim(fname_a),form='unformatted',access='direct',&
             recl= nc*nr*2)
        read(ftn,rec=1) ndvi_temp
        call LVT_releaseUnitNumber(ftn)
        
        input_bitmap = .false. 
        do r=1,nr
           do c=1,nc
              if(ndvi_temp(r,c).gt.0) then 
                 flag = ndvi_temp(r,c) - floor(ndvi_temp(r,c)/10.0)*10.0 + 1
                 if(flag.eq.1.or.flag.eq.2) then 
                    ndvi_inp1(c,nr-r+1) = & 
                         floor(ndvi_temp(r,c)/10.0)/1000.0
                 else
                    ndvi_inp1(c,nr-r+1) = -9999.0
                 endif
              else
                 ndvi_inp1(c,nr-r+1) = -9999.0
              endif
           enddo
        enddo

        do r=1,nr
           do c=1,nc
              ndvi_inp(c+(r-1)*nc) = ndvi_inp1(c,r)
              if(ndvi_inp(c+(r-1)*nc).gt.-9999.0) then 
                 input_bitmap(c+(r-1)*nc) = .true. 
              endif
           enddo
        enddo

        if(LVT_isAtAfinerResolution(gimmsndviobs(source)%datares)) then
           call neighbor_interp(LVT_rc%gridDesc, input_bitmap, &
                ndvi_inp, output_bitmap, ndvi_out_a, &
                nc*nr, &
                LVT_rc%lnc*LVT_rc%lnr, &
                gimmsndviobs(source)%rlat, & 
                gimmsndviobs(source)%rlon, &
                gimmsndviobs(source)%n11, &
                LVT_rc%udef, iret)
           
        else
           call upscaleByAveraging(&
                nc*nr, &
                LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
                gimmsndviobs(source)%n11, input_bitmap, &
                ndvi_inp, output_bitmap, ndvi_out_a)

        endif

        write(LVT_logunit,*) '[INFO] Reading GIMMS NDVI file ',trim(fname_b)
        ftn = LVT_getNextUnitNumber()
        open(ftn,file=trim(fname_b),form='unformatted',access='direct',&
             recl= nc*nr*2)
        read(ftn,rec=1) ndvi_temp
        call LVT_releaseUnitNumber(ftn)
        
        input_bitmap = .false. 
        do r=1,nr
           do c=1,nc
              if(ndvi_temp(r,c).gt.0) then 
                 flag = ndvi_temp(r,c) - floor(ndvi_temp(r,c)/10.0)*10.0 + 1
                 if(flag.eq.1.or.flag.eq.2) then 
                    ndvi_inp1(c,nr-r+1) = & 
                         floor(ndvi_temp(r,c)/10.0)/1000.0
                 else
                    ndvi_inp1(c,nr-r+1) = -9999.0
                 endif
              else
                 ndvi_inp1(c,nr-r+1) = -9999.0
              endif
           enddo
        enddo

        do r=1,nr
           do c=1,nc
              ndvi_inp(c+(r-1)*nc) = ndvi_inp1(c,r)
              if(ndvi_inp(c+(r-1)*nc).gt.-9999.0) then 
                 input_bitmap(c+(r-1)*nc) = .true. 
              endif
           enddo
        enddo

        if(LVT_isAtAfinerResolution(gimmsndviobs(source)%datares)) then
           call neighbor_interp(LVT_rc%gridDesc, input_bitmap, &
                ndvi_inp, output_bitmap, ndvi_out_b, &
                nc*nr, &
                LVT_rc%lnc*LVT_rc%lnr, &
                gimmsndviobs(source)%rlat, & 
                gimmsndviobs(source)%rlon, &
                gimmsndviobs(source)%n11, &
                LVT_rc%udef, iret)
           
        else
           call upscaleByAveraging(&
                nc*nr, &
                LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
                gimmsndviobs(source)%n11, input_bitmap, &
                ndvi_inp, output_bitmap, ndvi_out_b)

        endif
        
        do r=1, LVT_rc%lnr
           do c=1, LVT_rc%lnc
              if(ndvi_out_a(c+(r-1)*LVT_rc%lnc).gt.0.and.&
                   ndvi_out_b(c+(r-1)*LVT_rc%lnc).gt.0) then 
                 varfield(c,r) = (ndvi_out_a(c+(r-1)*LVT_rc%lnc)+& 
                      ndvi_out_b(c+(r-1)*LVT_rc%lnc))/2.0
              else
                 varfield(c,r) = LVT_rc%udef
              endif
           enddo
        enddo

     else
        varfield  = LVT_rc%udef
     endif
     
  endif
  
  call LVT_logSingleDataStreamVar(LVT_MOC_NDVI,source,varfield,&
       vlevel=1,units="-")
  
end subroutine readGIMMS_NDVIObs

!BOP
! 
! !ROUTINE: create_GIMMSndvi_filename
! \label{create_GIMMSndvi_filename}
!
! !INTERFACE: 
subroutine create_GIMMSndvi_filename(odir,yr,mo,filename_a, filename_b)
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
! This routine creates a timestamped filename for MODIS LST data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GlboSnow base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the GIMMSNDVI_LH file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr
  integer                      :: mo
  integer                      :: da
  character(len=*)             :: filename_a
  character(len=*)             :: filename_b
!
!EOP
  
  character*20            :: dir_name
  character*20            :: sat_name
  character*4             :: fyr
  character*2             :: fyr1, fyr2
  character*3             :: month(12), fmonth
  data month /'jan','feb','mar','apr','may','jun','jul','aug','sep',&
       'oct','nov','dec'/

  fmonth = month(mo)
  write(unit=fyr, fmt='(i4.4)') yr
  
  if(yr.eq.2014) then 
     dir_name = "2013s_new"     
     sat_name = "19"
  elseif(yr.eq.2013) then 
     dir_name = "2012s_new"     
     sat_name = "19"
  elseif(yr.eq.2012) then 
     dir_name = "2011s_new"     
     sat_name = "18"
  elseif(yr.eq.2011) then 
     dir_name = "2010s_new"     
     sat_name = "18"
  elseif(yr.eq.2000.and.mo.le.10) then 
     dir_name = "2000s_new"     
     sat_name = "14"
  elseif(yr.eq.2000.and.mo.gt.10) then 
     dir_name = "2000s_new"     
     sat_name = "15"
  elseif(yr.eq.2001.or.yr.eq.2002.or.yr.eq.2003) then 
     dir_name = "2000s_new"
     sat_name = "16"
  elseif(yr.eq.2004.or.yr.eq.2005.or.yr.eq.2006.or.&
       yr.eq.2007.or.yr.eq.2008) then 
     dir_name = "2000s_new"
     sat_name = "17"
  elseif(yr.eq.2009) then 
     dir_name = "2000s_new"
     sat_name = "18"
  elseif(yr.eq.1999.or.yr.eq.1998.or.yr.eq.1997.or.yr.eq.1996.or.&
       yr.eq.1995) then 
     dir_name = "1990s_new"
     sat_name = "14"
  elseif(yr.eq.1994.and.mo.le.8) then 
     dir_name = "1990s_new"
     sat_name = "11"
  elseif(yr.eq.1994.and.mo.gt.8) then 
     dir_name = "1990s_new"
     sat_name = "09"
  elseif(yr.eq.1993.or.yr.eq.1992.or.yr.eq.1991.or.yr.eq.1990) then 
     dir_name = "1990s_new"
     sat_name = "11"
  elseif(yr.eq.1989) then 
     dir_name = "1980s_new"
     sat_name = "11"
  elseif(yr.eq.1988.and.mo.le.10) then 
     dir_name = "1980s_new"
     sat_name = "09"
  elseif(yr.eq.1988.and.mo.gt.10) then 
     dir_name = "1980s_new"
     sat_name = "11"
  elseif(yr.eq.1987.or.yr.eq.1986) then 
     dir_name = "1980s_new"
     sat_name = "09"
  elseif(yr.eq.1985.and.mo.le.2) then 
     dir_name = "1980s_new"
     sat_name = "07"
  elseif(yr.eq.1985.and.mo.gt.2) then 
     dir_name = "1980s_new"
     sat_name = "09"
  elseif(yr.eq.1984.or.yr.eq.1983.or.yr.eq.1982) then 
     dir_name = "1980s_new"
     sat_name = "07"
  endif

  read(unit=fyr,fmt='(a2,a2)') fyr1, fyr2
  filename_a = trim(odir)//'/'//trim(dir_name)//&
       '/geo'//trim(fyr2)//trim(fmonth)//&
       '15a.n'//trim(sat_name)//'-VI3g'

  filename_b = trim(odir)//'/'//trim(dir_name)//&
       '/geo'//trim(fyr2)//trim(fmonth)//&
       '15b.n'//trim(sat_name)//'-VI3g'
  
end subroutine create_GIMMSndvi_filename


